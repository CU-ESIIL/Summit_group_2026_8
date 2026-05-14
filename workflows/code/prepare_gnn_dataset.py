#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import train_test_split
from torch_geometric.data import Data


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare microbial coalescence graph inputs for PyTorch Geometric."
    )
    parser.add_argument(
        "--graph-input-dir",
        default="workflows/output/graph_inputs",
        help="Directory containing CSV outputs from build_graph_inputs.R.",
    )
    parser.add_argument(
        "--out-dir",
        default="workflows/output/gnn_dataset",
        help="Directory to write PyG dataset outputs.",
    )
    parser.add_argument(
        "--edge-scope",
        choices=["resident", "final", "resident_final", "none"],
        default="resident",
        help="Which SpiecEasi taxon-taxon edges to use for each graph.",
    )
    parser.add_argument(
        "--split-seed",
        type=int,
        default=1,
        help="Random seed for train/validation/test splits.",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.20,
        help="Fraction of microcosms assigned to test set.",
    )
    parser.add_argument(
        "--val-size",
        type=float,
        default=0.20,
        help="Fraction of non-test microcosms assigned to validation set.",
    )
    parser.add_argument(
        "--scale-abundance-for-log",
        type=float,
        default=1e4,
        help="Scale abundance before log1p features.",
    )
    return parser.parse_args()


def load_inputs(graph_input_dir: Path) -> Dict[str, pd.DataFrame]:
    files = {
        "triplets": "coalescence_triplets.csv",
        "abundance_edges": "combined_sample_taxon_edges.csv",
        "taxa": "nodes_taxa.csv",
        "spieceasi_edges": "taxon_taxon_spieceasi_edges.csv",
        "run_summary": "spieceasi_run_summary.csv",
    }

    dfs = {}
    for key, filename in files.items():
        path = graph_input_dir / filename
        if not path.exists():
            raise FileNotFoundError(f"Required input file not found: {path}")
        dfs[key] = pd.read_csv(path)

    return dfs


def make_taxon_indices(taxa: pd.DataFrame) -> pd.DataFrame:
    taxa = taxa.copy()
    taxa = taxa.sort_values(["kingdom", "taxon_id"]).reset_index(drop=True)
    taxa["taxon_index_global"] = np.arange(len(taxa), dtype=np.int64)

    # Also make kingdom-specific indices, which are easier for per-kingdom graphs.
    taxa["taxon_index"] = (
        taxa.groupby("kingdom").cumcount().astype(np.int64)
    )

    return taxa


def build_abundance_lookup(abundance_edges: pd.DataFrame) -> Dict[Tuple[str, str, str, str], pd.DataFrame]:
    """
    Key:
      (kingdom, donor_id, community_type, sample_id)

    Value:
      DataFrame with taxon_id, abundance
    """
    lookup = {}
    grouped = abundance_edges.groupby(
        ["kingdom", "donor_id", "community_type", "sample_id"],
        sort=False,
    )

    for key, df in grouped:
        lookup[key] = df[["taxon_id", "abundance"]].copy()

    return lookup


def donor_source_profile(
    abundance_edges: pd.DataFrame,
    kingdom: str,
    donor_id: str,
) -> pd.Series:
    """
    Return mean donor abundance per taxon for a pooled donor source.

    In this workshop dataset, donor_sample_id is NA because donor communities
    are pooled source-level inputs. We average across donor rows for that
    kingdom and donor_id. If profiles are identical, this returns the same profile.
    """
    donor = abundance_edges[
        (abundance_edges["kingdom"] == kingdom)
        & (abundance_edges["donor_id"] == donor_id)
        & (abundance_edges["community_type"] == "donor")
    ]

    if donor.empty:
        return pd.Series(dtype=float)

    return donor.groupby("taxon_id")["abundance"].mean()


def sample_profile(
    lookup: Dict[Tuple[str, str, str, str], pd.DataFrame],
    kingdom: str,
    donor_id: str,
    community_type: str,
    sample_id: str,
) -> pd.Series:
    key = (kingdom, donor_id, community_type, sample_id)
    df = lookup.get(key)

    if df is None or df.empty:
        return pd.Series(dtype=float)

    return df.set_index("taxon_id")["abundance"]


def select_edges_for_graph(
    spieceasi_edges: pd.DataFrame,
    kingdom: str,
    donor_id: str,
    edge_scope: str,
) -> pd.DataFrame:
    if edge_scope == "none":
        return pd.DataFrame(columns=spieceasi_edges.columns)

    edges = spieceasi_edges[
        (spieceasi_edges["kingdom"] == kingdom)
        & (spieceasi_edges["donor_id"] == donor_id)
    ]

    if edge_scope == "resident":
        edges = edges[edges["community_type"] == "resident"]
    elif edge_scope == "final":
        edges = edges[edges["community_type"] == "final"]
    elif edge_scope == "resident_final":
        edges = edges[edges["community_type"].isin(["resident", "final"])]
    else:
        raise ValueError(f"Unsupported edge_scope: {edge_scope}")

    return edges.copy()


def make_edge_index_and_attr(
    edges: pd.DataFrame,
    taxon_to_index: Dict[str, int],
) -> Tuple[torch.Tensor, torch.Tensor]:
    if edges.empty:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0, 2), dtype=torch.float32)
        return edge_index, edge_attr

    # Keep edges where both taxa are in this kingdom's taxon index.
    mask = edges["taxon_a"].isin(taxon_to_index) & edges["taxon_b"].isin(taxon_to_index)
    edges = edges.loc[mask].copy()

    if edges.empty:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0, 2), dtype=torch.float32)
        return edge_index, edge_attr

    src = edges["taxon_a"].map(taxon_to_index).to_numpy(dtype=np.int64)
    dst = edges["taxon_b"].map(taxon_to_index).to_numpy(dtype=np.int64)

    # Undirected graph: add both directions.
    edge_index_np = np.vstack([
        np.concatenate([src, dst]),
        np.concatenate([dst, src]),
    ])

    weight = edges["weight"].to_numpy(dtype=np.float32)
    sign = edges["sign"].fillna(0).to_numpy(dtype=np.float32)

    edge_attr_np = np.vstack([
        np.concatenate([weight, weight]),
        np.concatenate([sign, sign]),
    ]).T

    return (
        torch.tensor(edge_index_np, dtype=torch.long),
        torch.tensor(edge_attr_np, dtype=torch.float32),
    )


def vectorize_profile(
    profile: pd.Series,
    taxon_ids: List[str],
) -> np.ndarray:
    return profile.reindex(taxon_ids).fillna(0.0).to_numpy(dtype=np.float32)


def build_graph_for_triplet(
    triplet: pd.Series,
    taxa_for_kingdom: pd.DataFrame,
    abundance_lookup: Dict[Tuple[str, str, str, str], pd.DataFrame],
    abundance_edges: pd.DataFrame,
    spieceasi_edges: pd.DataFrame,
    edge_scope: str,
    scale_for_log: float,
) -> Data:
    kingdom = triplet["kingdom"]
    donor_id = triplet["donor_id"]
    resident_sample_id = triplet["resident_sample_id"]
    final_sample_id = triplet["final_sample_id"]

    taxa_for_kingdom = taxa_for_kingdom.sort_values("taxon_index")
    taxon_ids = taxa_for_kingdom["taxon_id"].tolist()
    taxon_to_index = dict(zip(taxon_ids, taxa_for_kingdom["taxon_index"].tolist()))

    resident = sample_profile(
        abundance_lookup, kingdom, donor_id, "resident", resident_sample_id
    )
    final = sample_profile(
        abundance_lookup, kingdom, donor_id, "final", final_sample_id
    )
    donor = donor_source_profile(abundance_edges, kingdom, donor_id)

    resident_vec = vectorize_profile(resident, taxon_ids)
    donor_vec = vectorize_profile(donor, taxon_ids)
    final_vec = vectorize_profile(final, taxon_ids)

    resident_present = (resident_vec > 0).astype(np.float32)
    donor_present = (donor_vec > 0).astype(np.float32)
    final_present = (final_vec > 0).astype(np.float32)

    resident_log = np.log1p(resident_vec * scale_for_log).astype(np.float32)
    donor_log = np.log1p(donor_vec * scale_for_log).astype(np.float32)
    final_log = np.log1p(final_vec * scale_for_log).astype(np.float32)

    x_np = np.vstack([
        resident_vec,
        donor_vec,
        resident_present,
        donor_present,
        resident_log,
        donor_log,
    ]).T.astype(np.float32)

    selected_edges = select_edges_for_graph(
        spieceasi_edges,
        kingdom=kingdom,
        donor_id=donor_id,
        edge_scope=edge_scope,
    )
    edge_index, edge_attr = make_edge_index_and_attr(selected_edges, taxon_to_index)

    data = Data(
        x=torch.tensor(x_np, dtype=torch.float32),
        edge_index=edge_index,
        edge_attr=edge_attr,
        y_abundance=torch.tensor(final_vec, dtype=torch.float32),
        y_presence=torch.tensor(final_present, dtype=torch.float32),
        y_log_abundance=torch.tensor(final_log, dtype=torch.float32),
    )

    # Metadata stored as Python attributes.
    data.microcosm_id = triplet["microcosm_id"]
    data.kingdom = kingdom
    data.donor_id = donor_id
    data.resident_sample_id = resident_sample_id
    data.final_sample_id = final_sample_id
    data.donor_source_id = triplet["donor_source_id"]
    data.donor_is_pooled = bool(triplet["donor_is_pooled"])
    data.edge_scope = edge_scope
    data.taxon_ids = taxon_ids

    return data


def make_splits(metadata: pd.DataFrame, test_size: float, val_size: float, seed: int) -> Dict[str, List[int]]:
    """
    Stratify by kingdom and donor_id so every donor source is represented.
    """
    meta = metadata.copy()
    meta["stratum"] = meta["kingdom"].astype(str) + "::" + meta["donor_id"].astype(str)

    indices = np.arange(len(meta))

    train_val_idx, test_idx = train_test_split(
        indices,
        test_size=test_size,
        random_state=seed,
        stratify=meta["stratum"],
    )

    train_val_meta = meta.iloc[train_val_idx].reset_index(drop=False)

    train_local, val_local = train_test_split(
        np.arange(len(train_val_meta)),
        test_size=val_size,
        random_state=seed,
        stratify=train_val_meta["stratum"],
    )

    train_idx = train_val_meta.iloc[train_local]["index"].to_numpy(dtype=int)
    val_idx = train_val_meta.iloc[val_local]["index"].to_numpy(dtype=int)

    return {
        "train": train_idx.tolist(),
        "val": val_idx.tolist(),
        "test": test_idx.tolist(),
    }


def main():
    args = parse_args()

    graph_input_dir = Path(args.graph_input_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dfs = load_inputs(graph_input_dir)

    taxa = make_taxon_indices(dfs["taxa"])
    abundance_edges = dfs["abundance_edges"]
    triplets = dfs["triplets"]
    spieceasi_edges = dfs["spieceasi_edges"]

    abundance_lookup = build_abundance_lookup(abundance_edges)

    graphs: List[Data] = []
    metadata_rows = []

    for i, triplet in triplets.iterrows():
        kingdom = triplet["kingdom"]
        taxa_for_kingdom = taxa[taxa["kingdom"] == kingdom].copy()

        data = build_graph_for_triplet(
            triplet=triplet,
            taxa_for_kingdom=taxa_for_kingdom,
            abundance_lookup=abundance_lookup,
            abundance_edges=abundance_edges,
            spieceasi_edges=spieceasi_edges,
            edge_scope=args.edge_scope,
            scale_for_log=args.scale_abundance_for_log,
        )

        graphs.append(data)
        metadata_rows.append({
            "graph_index": i,
            "microcosm_id": triplet["microcosm_id"],
            "kingdom": triplet["kingdom"],
            "donor_id": triplet["donor_id"],
            "resident_sample_id": triplet["resident_sample_id"],
            "final_sample_id": triplet["final_sample_id"],
            "donor_source_id": triplet["donor_source_id"],
            "donor_is_pooled": triplet["donor_is_pooled"],
            "n_nodes": data.num_nodes,
            "n_edges": data.edge_index.shape[1],
            "edge_scope": args.edge_scope,
            "n_resident_present": int(data.x[:, 2].sum().item()),
            "n_donor_present": int(data.x[:, 3].sum().item()),
            "n_final_present": int(data.y_presence.sum().item()),
        })

    metadata = pd.DataFrame(metadata_rows)
    splits = make_splits(
        metadata,
        test_size=args.test_size,
        val_size=args.val_size,
        seed=args.split_seed,
    )

    torch.save(graphs, out_dir / "graphs.pt")
    taxa.to_csv(out_dir / "taxon_index.csv", index=False)
    metadata.to_csv(out_dir / "graph_metadata.csv", index=False)

    with open(out_dir / "split_manifest.json", "w") as f:
        json.dump(splits, f, indent=2)

    feature_manifest = {
        "node_features": [
            "resident_abundance",
            "donor_abundance",
            "resident_present",
            "donor_present",
            "resident_log1p_scaled_abundance",
            "donor_log1p_scaled_abundance",
        ],
        "edge_features": [
            "spieceasi_weight",
            "spieceasi_sign",
        ],
        "targets": [
            "final_abundance",
            "final_presence",
            "final_log1p_scaled_abundance",
        ],
        "edge_scope": args.edge_scope,
        "donor_note": "Donor abundance is pooled by kingdom and donor_id using the mean donor-source profile.",
    }

    with open(out_dir / "feature_manifest.json", "w") as f:
        json.dump(feature_manifest, f, indent=2)

    readme = f"""# GNN Dataset

Generated from `{graph_input_dir}`.

Each graph is one coalescence prediction unit from `coalescence_triplets.csv`.

## Graph definition

- Node = taxon
- Node features = resident and donor abundance/presence
- Edges = SpiecEasi taxon-taxon associations
- Target = final abundance and final presence

## Edge scope

`{args.edge_scope}`

## Files

- `graphs.pt`: list of PyTorch Geometric Data objects
- `graph_metadata.csv`: one row per graph
- `taxon_index.csv`: taxon IDs and node indices
- `split_manifest.json`: train/validation/test graph indices
- `feature_manifest.json`: feature and target definitions

## Important caveat

Donor input is pooled at the donor-source level for this workshop dataset.
"""
    (out_dir / "README.md").write_text(readme)

    print(f"Wrote {len(graphs)} graphs to {out_dir}")
    print(metadata.groupby(["kingdom", "donor_id"]).size())


if __name__ == "__main__":
    main()