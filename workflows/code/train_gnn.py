#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GATv2Conv, GINEConv, SAGEConv

try:
    from sklearn.metrics import average_precision_score, roc_auc_score
except ImportError:  # pragma: no cover - sklearn is expected in gnn-env.
    average_precision_score = None
    roc_auc_score = None


EXCLUDE_BATCH_KEYS = [
    "taxon_ids",
    "microcosm_id",
    "kingdom",
    "donor_id",
    "resident_sample_id",
    "final_sample_id",
    "donor_source_id",
    "donor_is_pooled",
    "edge_scope",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train a baseline PyTorch Geometric GNN on microbial coalescence graphs."
    )
    parser.add_argument(
        "--dataset-dir",
        default="workflows/output/gnn_dataset",
        help="Directory containing graphs.pt and dataset manifests.",
    )
    parser.add_argument(
        "--out-dir",
        default="workflows/output/gnn_training",
        help="Directory for checkpoints, metrics, and run summaries.",
    )
    parser.add_argument(
        "--kingdom",
        default="all",
        help="Subset to one kingdom, such as Fungi or Bacteria. Use all for both.",
    )
    parser.add_argument(
        "--model",
        choices=["sage", "gine", "gatv2"],
        default="sage",
        help="GNN backbone. sage ignores edge attributes; gine and gatv2 use them.",
    )
    parser.add_argument("--hidden-dim", type=int, default=64)
    parser.add_argument("--num-layers", type=int, default=2)
    parser.add_argument("--dropout", type=float, default=0.20)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--patience", type=int, default=15)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, cuda, or a specific device such as cuda:0.",
    )
    parser.add_argument(
        "--abundance-loss-scope",
        choices=["present", "all"],
        default="present",
        help="Compute log-abundance MSE on final-present taxa only or on all taxa.",
    )
    parser.add_argument("--abundance-loss-weight", type=float, default=1.0)
    parser.add_argument("--presence-loss-weight", type=float, default=1.0)
    parser.add_argument(
        "--presence-pos-weight",
        default="auto",
        help="BCE positive-class weight. Use auto or a positive number.",
    )
    parser.add_argument(
        "--max-pos-weight",
        type=float,
        default=100.0,
        help="Upper bound when --presence-pos-weight auto is used.",
    )
    parser.add_argument(
        "--no-standardize-inputs",
        action="store_true",
        help="Disable train-split feature standardization.",
    )
    parser.add_argument(
        "--max-graphs-per-split",
        type=int,
        default=None,
        help="Optional smoke-test limit applied separately to train/val/test splits.",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=0,
        help="PyG DataLoader worker count.",
    )
    parser.add_argument(
        "--save-test-predictions",
        action="store_true",
        help="Save test-set node predictions to test_predictions.pt.",
    )
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def choose_device(device_arg: str) -> torch.device:
    if device_arg == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(device_arg)


def load_json(path: Path) -> Dict:
    with path.open() as handle:
        return json.load(handle)


def load_graphs(path: Path) -> List:
    try:
        return torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        return torch.load(path, map_location="cpu")


def filter_split_indices(
    split_manifest: Dict[str, List[int]],
    metadata: pd.DataFrame,
    kingdom: str,
    max_graphs_per_split: Optional[int],
) -> Dict[str, List[int]]:
    if kingdom == "all":
        allowed = set(metadata["graph_index"].astype(int).tolist())
    else:
        allowed = set(
            metadata.loc[
                metadata["kingdom"].astype(str).str.lower() == kingdom.lower(),
                "graph_index",
            ]
            .astype(int)
            .tolist()
        )
        if not allowed:
            available = ", ".join(sorted(metadata["kingdom"].astype(str).unique()))
            raise ValueError(f"No graphs found for kingdom={kingdom!r}. Available: {available}")

    filtered = {}
    for split_name in ["train", "val", "test"]:
        indices = [int(i) for i in split_manifest.get(split_name, []) if int(i) in allowed]
        if max_graphs_per_split is not None:
            indices = indices[:max_graphs_per_split]
        filtered[split_name] = indices

    for split_name, indices in filtered.items():
        if not indices:
            raise ValueError(f"Split {split_name!r} is empty after filtering.")

    return filtered


def select_graphs(graphs: Sequence, indices: Sequence[int]) -> List:
    return [graphs[i] for i in indices]


def make_loader(
    graphs: Sequence,
    batch_size: int,
    shuffle: bool,
    num_workers: int,
) -> DataLoader:
    return DataLoader(
        list(graphs),
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        exclude_keys=EXCLUDE_BATCH_KEYS,
    )


def stream_feature_stats(graphs: Sequence) -> Tuple[torch.Tensor, torch.Tensor]:
    if not graphs:
        raise ValueError("Cannot compute feature statistics from an empty graph list.")

    n_features = int(graphs[0].x.shape[1])
    count = 0
    total = torch.zeros(n_features, dtype=torch.float64)
    total_sq = torch.zeros(n_features, dtype=torch.float64)

    for graph in graphs:
        x = graph.x.to(dtype=torch.float64)
        total += x.sum(dim=0)
        total_sq += (x * x).sum(dim=0)
        count += int(x.shape[0])

    mean = total / count
    var = (total_sq / count) - (mean * mean)
    std = torch.sqrt(torch.clamp(var, min=1e-12))
    std = torch.where(std == 0, torch.ones_like(std), std)
    return mean.to(dtype=torch.float32), std.to(dtype=torch.float32)


def count_presence_targets(graphs: Sequence) -> Tuple[float, float]:
    positives = 0.0
    total = 0.0
    for graph in graphs:
        y = graph.y_presence
        positives += float(y.sum().item())
        total += float(y.numel())
    return positives, total - positives


def parse_pos_weight(args: argparse.Namespace, train_graphs: Sequence) -> float:
    if args.presence_pos_weight != "auto":
        value = float(args.presence_pos_weight)
        if value <= 0:
            raise ValueError("--presence-pos-weight must be positive.")
        return value

    positives, negatives = count_presence_targets(train_graphs)
    if positives <= 0:
        return 1.0
    return min(negatives / positives, args.max_pos_weight)


class MicrocosmGNN(nn.Module):
    def __init__(
        self,
        in_channels: int,
        hidden_dim: int,
        num_layers: int,
        dropout: float,
        model_type: str,
        edge_dim: int,
        feature_mean: Optional[torch.Tensor],
        feature_std: Optional[torch.Tensor],
    ) -> None:
        super().__init__()
        if num_layers < 1:
            raise ValueError("num_layers must be at least 1.")

        self.model_type = model_type
        self.dropout = dropout

        if feature_mean is None or feature_std is None:
            self.register_buffer("feature_mean", None)
            self.register_buffer("feature_std", None)
        else:
            self.register_buffer("feature_mean", feature_mean.clone().detach())
            self.register_buffer("feature_std", feature_std.clone().detach())

        self.input_proj = nn.Linear(in_channels, hidden_dim)
        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()

        for _ in range(num_layers):
            if model_type == "sage":
                conv = SAGEConv(hidden_dim, hidden_dim)
            elif model_type == "gine":
                mlp = nn.Sequential(
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.ReLU(),
                    nn.Linear(hidden_dim, hidden_dim),
                )
                conv = GINEConv(mlp, edge_dim=edge_dim)
            elif model_type == "gatv2":
                conv = GATv2Conv(
                    hidden_dim,
                    hidden_dim,
                    heads=2,
                    concat=False,
                    edge_dim=edge_dim,
                    dropout=dropout,
                )
            else:
                raise ValueError(f"Unsupported model type: {model_type}")

            self.convs.append(conv)
            self.norms.append(nn.LayerNorm(hidden_dim))

        self.abundance_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )
        self.presence_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def standardize(self, x: torch.Tensor) -> torch.Tensor:
        if self.feature_mean is None or self.feature_std is None:
            return x
        return (x - self.feature_mean) / self.feature_std

    def forward(self, data) -> Tuple[torch.Tensor, torch.Tensor]:
        x = self.standardize(data.x)
        h = F.relu(self.input_proj(x))
        h = F.dropout(h, p=self.dropout, training=self.training)

        for conv, norm in zip(self.convs, self.norms):
            residual = h
            if self.model_type == "sage":
                h = conv(h, data.edge_index)
            else:
                h = conv(h, data.edge_index, data.edge_attr)
            h = F.relu(norm(h + residual))
            h = F.dropout(h, p=self.dropout, training=self.training)

        log_abundance = self.abundance_head(h).squeeze(-1)
        presence_logits = self.presence_head(h).squeeze(-1)
        return log_abundance, presence_logits


def compute_loss(
    log_abundance: torch.Tensor,
    presence_logits: torch.Tensor,
    batch,
    pos_weight: torch.Tensor,
    abundance_loss_scope: str,
    abundance_loss_weight: float,
    presence_loss_weight: float,
) -> Tuple[torch.Tensor, Dict[str, float]]:
    y_log = batch.y_log_abundance.view(-1).to(log_abundance.dtype)
    y_presence = batch.y_presence.view(-1).to(log_abundance.dtype)

    if abundance_loss_scope == "present":
        abundance_mask = y_presence > 0.5
    else:
        abundance_mask = torch.ones_like(y_presence, dtype=torch.bool)

    if abundance_mask.any():
        abundance_loss = F.mse_loss(log_abundance[abundance_mask], y_log[abundance_mask])
    else:
        abundance_loss = torch.zeros((), device=log_abundance.device)

    presence_loss = F.binary_cross_entropy_with_logits(
        presence_logits,
        y_presence,
        pos_weight=pos_weight,
    )
    loss = abundance_loss_weight * abundance_loss + presence_loss_weight * presence_loss
    parts = {
        "loss": float(loss.detach().cpu().item()),
        "abundance_loss": float(abundance_loss.detach().cpu().item()),
        "presence_loss": float(presence_loss.detach().cpu().item()),
    }
    return loss, parts


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
    pos_weight: torch.Tensor,
    args: argparse.Namespace,
) -> Dict[str, float]:
    model.train()
    totals = {"loss": 0.0, "abundance_loss": 0.0, "presence_loss": 0.0}
    batches = 0

    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad(set_to_none=True)
        log_abundance, presence_logits = model(batch)
        loss, parts = compute_loss(
            log_abundance,
            presence_logits,
            batch,
            pos_weight,
            args.abundance_loss_scope,
            args.abundance_loss_weight,
            args.presence_loss_weight,
        )
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()

        for key in totals:
            totals[key] += parts[key]
        batches += 1

    return {f"train_{key}": value / max(1, batches) for key, value in totals.items()}


def safe_auc(labels: np.ndarray, scores: np.ndarray) -> Optional[float]:
    if roc_auc_score is None or len(np.unique(labels)) < 2:
        return None
    return float(roc_auc_score(labels, scores))


def safe_average_precision(labels: np.ndarray, scores: np.ndarray) -> Optional[float]:
    if average_precision_score is None or len(np.unique(labels)) < 2:
        return None
    return float(average_precision_score(labels, scores))


@torch.no_grad()
def evaluate(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
    pos_weight: torch.Tensor,
    args: argparse.Namespace,
    collect_predictions: bool = False,
) -> Tuple[Dict[str, float], Optional[Dict[str, torch.Tensor]]]:
    model.eval()
    loss_parts = {"loss": 0.0, "abundance_loss": 0.0, "presence_loss": 0.0}
    batches = 0

    pred_log_chunks = []
    target_log_chunks = []
    pred_presence_chunks = []
    target_presence_chunks = []

    for batch in loader:
        batch = batch.to(device)
        log_abundance, presence_logits = model(batch)
        _, parts = compute_loss(
            log_abundance,
            presence_logits,
            batch,
            pos_weight,
            args.abundance_loss_scope,
            args.abundance_loss_weight,
            args.presence_loss_weight,
        )
        for key in loss_parts:
            loss_parts[key] += parts[key]
        batches += 1

        pred_log_chunks.append(log_abundance.detach().cpu())
        target_log_chunks.append(batch.y_log_abundance.view(-1).detach().cpu())
        pred_presence_chunks.append(torch.sigmoid(presence_logits).detach().cpu())
        target_presence_chunks.append(batch.y_presence.view(-1).detach().cpu())

    pred_log = torch.cat(pred_log_chunks)
    target_log = torch.cat(target_log_chunks)
    pred_presence = torch.cat(pred_presence_chunks)
    target_presence = torch.cat(target_presence_chunks)

    present_mask = target_presence > 0.5
    diff = pred_log - target_log
    rmse_all = torch.sqrt(torch.mean(diff * diff)).item()
    mae_all = torch.mean(torch.abs(diff)).item()

    if present_mask.any():
        present_diff = diff[present_mask]
        rmse_present = torch.sqrt(torch.mean(present_diff * present_diff)).item()
        mae_present = torch.mean(torch.abs(present_diff)).item()
    else:
        rmse_present = math.nan
        mae_present = math.nan

    predicted_class = pred_presence >= 0.5
    accuracy = (predicted_class == (target_presence >= 0.5)).to(torch.float32).mean().item()

    labels_np = target_presence.numpy()
    scores_np = pred_presence.numpy()
    metrics = {
        key: value / max(1, batches) for key, value in loss_parts.items()
    }
    metrics.update(
        {
            "log_rmse_all": float(rmse_all),
            "log_mae_all": float(mae_all),
            "log_rmse_present": float(rmse_present),
            "log_mae_present": float(mae_present),
            "presence_accuracy": float(accuracy),
            "presence_rate": float(target_presence.mean().item()),
            "presence_predicted_rate": float(predicted_class.to(torch.float32).mean().item()),
        }
    )

    auc = safe_auc(labels_np, scores_np)
    ap = safe_average_precision(labels_np, scores_np)
    if auc is not None:
        metrics["presence_auroc"] = auc
    if ap is not None:
        metrics["presence_average_precision"] = ap

    predictions = None
    if collect_predictions:
        predictions = {
            "pred_log_abundance": pred_log,
            "target_log_abundance": target_log,
            "pred_presence_probability": pred_presence,
            "target_presence": target_presence,
        }

    return metrics, predictions


def resident_baseline_metrics(graphs: Sequence) -> Dict[str, float]:
    pred_chunks = []
    target_chunks = []
    presence_chunks = []

    for graph in graphs:
        pred_chunks.append(graph.x[:, 4].detach().cpu())
        target_chunks.append(graph.y_log_abundance.detach().cpu())
        presence_chunks.append(graph.y_presence.detach().cpu())

    pred = torch.cat(pred_chunks)
    target = torch.cat(target_chunks)
    presence = torch.cat(presence_chunks) > 0.5

    diff = pred - target
    metrics = {
        "resident_log_rmse_all": float(torch.sqrt(torch.mean(diff * diff)).item()),
        "resident_log_mae_all": float(torch.mean(torch.abs(diff)).item()),
    }
    if presence.any():
        present_diff = diff[presence]
        metrics["resident_log_rmse_present"] = float(
            torch.sqrt(torch.mean(present_diff * present_diff)).item()
        )
        metrics["resident_log_mae_present"] = float(torch.mean(torch.abs(present_diff)).item())
    return metrics


def write_json(path: Path, payload: Dict) -> None:
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2)


def write_metrics_csv(path: Path, rows: List[Dict[str, float]]) -> None:
    if not rows:
        return
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_run_readme(
    out_dir: Path,
    args: argparse.Namespace,
    split_indices: Dict[str, List[int]],
    best_epoch: int,
    best_val_loss: float,
    test_metrics: Dict[str, float],
) -> None:
    readme = f"""# GNN Training Run

Generated by `workflows/code/train_gnn.py`.

## Inputs

- Dataset directory: `{args.dataset_dir}`
- Kingdom filter: `{args.kingdom}`
- Model: `{args.model}`
- Splits: train={len(split_indices["train"])}, val={len(split_indices["val"])}, test={len(split_indices["test"])}

## Outputs

- `best_model.pt`: checkpoint from the best validation-loss epoch
- `metrics.csv`: per-epoch train and validation metrics
- `run_config.json`: command-line configuration and split indices
- `test_metrics.json`: held-out test metrics from the best checkpoint

## Best Checkpoint

- Epoch: {best_epoch}
- Validation loss: {best_val_loss:.6g}
- Test loss: {test_metrics.get("loss", float("nan")):.6g}
- Test presence AUROC: {test_metrics.get("presence_auroc", float("nan")):.6g}
- Test log RMSE on present taxa: {test_metrics.get("log_rmse_present", float("nan")):.6g}
"""
    (out_dir / "README.md").write_text(readme)


def main() -> None:
    args = parse_args()
    set_seed(args.seed)

    dataset_dir = Path(args.dataset_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    graphs_path = dataset_dir / "graphs.pt"
    split_path = dataset_dir / "split_manifest.json"
    metadata_path = dataset_dir / "graph_metadata.csv"
    feature_manifest_path = dataset_dir / "feature_manifest.json"

    if not graphs_path.exists():
        raise FileNotFoundError(f"Missing graph file: {graphs_path}")
    if not split_path.exists():
        raise FileNotFoundError(f"Missing split manifest: {split_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing graph metadata: {metadata_path}")

    print(f"Loading graphs from {graphs_path}")
    graphs = load_graphs(graphs_path)
    metadata = pd.read_csv(metadata_path)
    split_manifest = load_json(split_path)
    feature_manifest = load_json(feature_manifest_path) if feature_manifest_path.exists() else {}

    split_indices = filter_split_indices(
        split_manifest,
        metadata,
        args.kingdom,
        args.max_graphs_per_split,
    )
    train_graphs = select_graphs(graphs, split_indices["train"])
    val_graphs = select_graphs(graphs, split_indices["val"])
    test_graphs = select_graphs(graphs, split_indices["test"])

    in_channels = int(train_graphs[0].x.shape[1])
    edge_dim = int(train_graphs[0].edge_attr.shape[1]) if train_graphs[0].edge_attr is not None else 0
    feature_mean = feature_std = None
    if not args.no_standardize_inputs:
        feature_mean, feature_std = stream_feature_stats(train_graphs)

    pos_weight_value = parse_pos_weight(args, train_graphs)
    device = choose_device(args.device)
    pos_weight = torch.tensor(pos_weight_value, dtype=torch.float32, device=device)

    train_loader = make_loader(train_graphs, args.batch_size, shuffle=True, num_workers=args.num_workers)
    val_loader = make_loader(val_graphs, args.batch_size, shuffle=False, num_workers=args.num_workers)
    test_loader = make_loader(test_graphs, args.batch_size, shuffle=False, num_workers=args.num_workers)

    model = MicrocosmGNN(
        in_channels=in_channels,
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers,
        dropout=args.dropout,
        model_type=args.model,
        edge_dim=edge_dim,
        feature_mean=feature_mean,
        feature_std=feature_std,
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.lr,
        weight_decay=args.weight_decay,
    )

    config = {
        "args": vars(args),
        "device": str(device),
        "torch_version": torch.__version__,
        "feature_manifest": feature_manifest,
        "split_indices": split_indices,
        "n_graphs": {key: len(value) for key, value in split_indices.items()},
        "input_channels": in_channels,
        "edge_dim": edge_dim,
        "presence_pos_weight": pos_weight_value,
        "standardized_inputs": not args.no_standardize_inputs,
    }
    write_json(out_dir / "run_config.json", config)

    print(
        "Training "
        f"{args.model} on {device} with splits "
        f"train={len(train_graphs)}, val={len(val_graphs)}, test={len(test_graphs)}"
    )
    print(f"Presence positive-class weight: {pos_weight_value:.4g}")

    best_val_loss = math.inf
    best_epoch = 0
    epochs_without_improvement = 0
    metric_rows: List[Dict[str, float]] = []
    checkpoint_path = out_dir / "best_model.pt"

    for epoch in range(1, args.epochs + 1):
        train_metrics = train_one_epoch(
            model,
            train_loader,
            optimizer,
            device,
            pos_weight,
            args,
        )
        val_metrics, _ = evaluate(
            model,
            val_loader,
            device,
            pos_weight,
            args,
        )
        row = {"epoch": epoch}
        row.update(train_metrics)
        row.update({f"val_{key}": value for key, value in val_metrics.items()})
        metric_rows.append(row)
        write_metrics_csv(out_dir / "metrics.csv", metric_rows)

        val_loss = val_metrics["loss"]
        print(
            f"epoch={epoch:03d} "
            f"train_loss={train_metrics['train_loss']:.5f} "
            f"val_loss={val_loss:.5f} "
            f"val_auc={val_metrics.get('presence_auroc', float('nan')):.4f} "
            f"val_rmse_present={val_metrics.get('log_rmse_present', float('nan')):.4f}"
        )

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_epoch = epoch
            epochs_without_improvement = 0
            torch.save(
                {
                    "epoch": epoch,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "val_metrics": val_metrics,
                    "config": config,
                },
                checkpoint_path,
            )
        else:
            epochs_without_improvement += 1

        if epochs_without_improvement >= args.patience:
            print(f"Stopping early after {args.patience} epochs without validation improvement.")
            break

    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint["model_state_dict"])
    test_metrics, test_predictions = evaluate(
        model,
        test_loader,
        device,
        pos_weight,
        args,
        collect_predictions=args.save_test_predictions,
    )
    test_metrics.update(resident_baseline_metrics(test_graphs))
    write_json(out_dir / "test_metrics.json", test_metrics)

    if args.save_test_predictions and test_predictions is not None:
        torch.save(test_predictions, out_dir / "test_predictions.pt")

    write_run_readme(
        out_dir=out_dir,
        args=args,
        split_indices=split_indices,
        best_epoch=best_epoch,
        best_val_loss=best_val_loss,
        test_metrics=test_metrics,
    )

    print(f"Best epoch: {best_epoch}")
    print(f"Test loss: {test_metrics['loss']:.5f}")
    if "presence_auroc" in test_metrics:
        print(f"Test presence AUROC: {test_metrics['presence_auroc']:.4f}")
    print(f"Wrote training artifacts to {out_dir}")


if __name__ == "__main__":
    main()
