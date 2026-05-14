# Graph Inputs Workflow

## What This Workflow Does

This workflow discovers microbial coalescence community matrices and taxonomy tables, validates their structure, builds a bipartite sample-taxon backbone, fits stratum-specific SpiecEasi taxon-taxon networks, and writes graph-ready CSV files for downstream graph neural network workflows.

## Inputs Detected

Community matrices were detected from `workflows/input/` first and then the repository root/immediate subdirectories as fallback.
- `Bacteria_inoculation_experiment_0Burn_W1_donor-community.csv`
- `Bacteria_inoculation_experiment_0Burn_W1_final-community.csv`
- `Bacteria_inoculation_experiment_0Burn_W1_resident-community.csv`
- `Bacteria_inoculation_experiment_0Burn_W2_donor-community.csv`
- `Bacteria_inoculation_experiment_0Burn_W2_final-community.csv`
- `Bacteria_inoculation_experiment_0Burn_W2_resident-community.csv`
- `Bacteria_inoculation_experiment_1Burn_W3_donor-community.csv`
- `Bacteria_inoculation_experiment_1Burn_W3_final-community.csv`
- `Bacteria_inoculation_experiment_1Burn_W3_resident-community.csv`
- `Bacteria_inoculation_experiment_1Burn_W4_donor-community.csv`
- `Bacteria_inoculation_experiment_1Burn_W4_final-community.csv`
- `Bacteria_inoculation_experiment_1Burn_W4_resident-community.csv`
- `Bacteria_inoculation_experiment_3Burn_W5_donor-community.csv`
- `Bacteria_inoculation_experiment_3Burn_W5_final-community.csv`
- `Bacteria_inoculation_experiment_3Burn_W5_resident-community.csv`
- `Bacteria_inoculation_experiment_3Burn_W6_donor-community.csv`
- `Bacteria_inoculation_experiment_3Burn_W6_final-community.csv`
- `Bacteria_inoculation_experiment_3Burn_W6_resident-community.csv`
- `Fungi_inoculation_experiment_0Burn_W1_donor-community.csv`
- `Fungi_inoculation_experiment_0Burn_W1_final-community.csv`
- `Fungi_inoculation_experiment_0Burn_W1_resident-community.csv`
- `Fungi_inoculation_experiment_0Burn_W2_donor-community.csv`
- `Fungi_inoculation_experiment_0Burn_W2_final-community.csv`
- `Fungi_inoculation_experiment_0Burn_W2_resident-community.csv`
- `Fungi_inoculation_experiment_1Burn_W3_donor-community.csv`
- `Fungi_inoculation_experiment_1Burn_W3_final-community.csv`
- `Fungi_inoculation_experiment_1Burn_W3_resident-community.csv`
- `Fungi_inoculation_experiment_1Burn_W4_donor-community.csv`
- `Fungi_inoculation_experiment_1Burn_W4_final-community.csv`
- `Fungi_inoculation_experiment_1Burn_W4_resident-community.csv`
- `Fungi_inoculation_experiment_3Burn_W5_donor-community.csv`
- `Fungi_inoculation_experiment_3Burn_W5_final-community.csv`
- `Fungi_inoculation_experiment_3Burn_W5_resident-community.csv`
- `Fungi_inoculation_experiment_3Burn_W6_donor-community.csv`
- `Fungi_inoculation_experiment_3Burn_W6_final-community.csv`
- `Fungi_inoculation_experiment_3Burn_W6_resident-community.csv`

Taxonomy tables used:
- `Bacteria`: `Bacteria_inoculation_experiment_taxonomy_table.csv`
- `Fungi`: `Fungi_inoculation_experiment_taxonomy_table.csv`

## Filename Metadata Parsing

Community filenames were parsed with the pattern `Kingdom_inoculation_experiment_<donor_id>_<community_type>-community.csv`.
This workflow uses:
- `kingdom`: `Bacteria` or `Fungi`
- `donor_id`: the middle filename token such as `0Burn_W1`
- `community_type`: `donor`, `resident`, or `final`

## Sample-Taxon Bipartite Backbone

Each community matrix is treated as samples in rows and taxa/features in columns.
The first column is used as `sample_id`; if it was not already named `sample_id`, the workflow renames it and records that warning in the summary outputs.
Only nonzero abundance entries are written to the bipartite edge table.

## Taxon Naming

Taxon node metadata comes from the kingdom-specific taxonomy tables when available.
Scientific names are constructed from `Genus` and `Species` only.
- `name = "Genus Species"` when both are present
- `name = "Genus"` when only genus is present
- `name = ""` when genus is missing
Higher taxonomy is retained only as optional metadata fields, not as the primary name.

## Coalescence Triplets

The workflow also writes `coalescence_triplets.csv`, which defines the supervised prediction units for this workshop dataset.
Resident and final communities are paired when they share the same `sample_id` within a given `kingdom` and `donor_id`.
Donor communities are represented at the pooled donor-source level: all resident/final samples with the same `donor_id` share the same donor-source input.
Accordingly, `donor_sample_id` is `NA`, `donor_source_id` stores the donor treatment/source, and `donor_is_pooled` is `TRUE`.

## SpiecEasi Networks

SpiecEasi is run separately for each `kingdom x donor_id x community_type` stratum so that bacteria and fungi are not pooled, donor/resident/final communities are not pooled, and donor IDs are not pooled.
The configured SpiecEasi method is `mb` with selection criterion `bstars`.
When a stratum retains more taxa than the workshop cap allows, taxa are ranked by prevalence and then mean abundance before taking the top subset for SpiecEasi. This cap affects only the inferred taxon-taxon layer; the full bipartite backbone is preserved.
For `method = mb`, SpiecEasi uses neighborhood selection. The resulting taxon-taxon edges should be interpreted as inferred association structure, not direct ecological interactions or strict partial correlations.
For `method = glasso`, the selected inverse covariance structure can be converted to partial-correlation-like associations. This script stores edge weights as absolute association strength and `sign` as the association direction from the selected SpiecEasi matrix.

Workflow parameters:
- `min_prevalence = 0.05`
- `min_taxa = 10`
- `scale_to_counts = TRUE`
- `scale_factor = 10000`
- `max_taxa_for_spieceasi = Inf`
- `spieceasi_time_limit_seconds = Inf`
- `spieceasi_method = "mb"`
- `spieceasi_sel_criterion = "bstars"`
- `spieceasi_ncores = 15`
- `nlambda = 30`
- `lambda_min_ratio = 0.01`
- `rep_num = 20`
- `random_seed = 1`

## Output Files

- `combined_sample_taxon_edges.csv`: nonzero sample-taxon bipartite edges with abundance and experimental context.
- `nodes_samples.csv`: unique sample nodes with kingdom, donor ID, and community type.
- `coalescence_triplets.csv`: supervised coalescence units linking resident and final samples by shared sample ID within kingdom and donor source; donor input is represented by pooled donor source ID.
- `nodes_taxa.csv`: unique taxon nodes with names and taxonomy metadata.
- `taxon_taxon_spieceasi_edges.csv`: undirected SpiecEasi taxon-taxon edges with inferred association weights and signs.
- `graph_edges_multirelational.csv`: flat edge file combining sample-taxon and taxon-taxon relations.
- `spieceasi_run_summary.csv`: one row per stratum describing network status and skip/failure messages.

## Rerun

Run the workflow from the repository root with:

```bash
Rscript workflows/code/build_graph_inputs.R
```

## GNN Use

The bipartite backbone captures observed experimental composition data, while the SpiecEasi layer adds inferred ecological association structure. The combined multirelational edge file can be used as a starting point for heterogeneous or relational GNN pipelines that link samples, taxa, and inferred taxon-taxon associations.

## Caveats

- Co-occurrence edges are inferred statistical associations, not measured direct interactions.
- Positive co-occurrence does not necessarily mean cooperation.
- Negative co-occurrence does not necessarily mean inhibition.
- Associations can reflect shared niches, environmental filtering, compositional effects, or indirect interactions.
- Prevalence filtering affects network density.
- Relative-abundance to pseudo-count conversion is a modeling choice.
- The bipartite backbone is the experimental data representation; the SpiecEasi graph is an inferred ecological association layer.
- With `method = mb`, edge weights are association-strength summaries from neighborhood selection and should not be described as strict partial correlations.
- Donor-community SpiecEasi strata may be intentionally skipped if repeated donor profiles have no sample-to-sample variation; donor composition is still represented in the sample-taxon backbone.

## Run Summary

- Sample-taxon edges: `943,126`
- Sample nodes: `288`
- Coalescence triplets: `288`
- Taxon nodes: `19,523`
- SpiecEasi edges: `144,436`
- Successful SpiecEasi strata: `24`
- Skipped SpiecEasi strata: `12`
- Failed SpiecEasi strata: `0`
