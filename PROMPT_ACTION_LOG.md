Prompt Action Log

YYYY-MM-DD

Prompt

You are working inside a workshop repository that includes a reporting website, documentation, and template materials. Do not reorganize the repo. Do not touch `docs/`, `site/`, `scripts/`, `code/`, or website files unless absolutely necessary.

All new workflow code must go under:

- `./workflows/code/`

All workflow inputs should be searched for under:

- `./workflows/input/`
- and, if not found there, the repository root and immediate subdirectories

All generated workflow outputs must go under:

- `./workflows/output/`

There is existing template material in `./workflows/`; it is okay to add new scripts and documentation there, but do not delete anything unless I explicitly ask.

Project context:
This is a microbial coalescence experiment being prepared for graph neural network modeling at a workshop. The goal is to generate graph-ready files from community composition data.

The biological setup:
- Donor microbial communities are introduced into resident communities.
- Final communities are observed after coalescence.
- The eventual GNN task is likely to predict final community states from donor and resident communities.
- We want two graph structures:
  1. A bipartite sample-taxon backbone, encoding abundance and experimental context.
  2. Taxon-taxon SpiecEasi co-occurrence networks, encoding inferred interspecies association structure.
- These should be merged into a flat multirelational edge file for downstream GNN use.

Your task:
Build a reproducible R-based workflow in `./workflows/code/` that discovers the input files, validates them, generates graph files, and writes a clear README/report under `./workflows/output/`.

Before coding:
1. Inspect the repository.
2. Inspect `./workflows/`, especially `./workflows/code/`, `./workflows/input/`, and `./workflows/output/`.
3. Search for microbial coalescence input files. Prefer `./workflows/input/`, but fall back to the repository root and nearby folders if needed.

Expected input files:
Community matrices may look like:

- `Fungi_inoculation_experiment_*_donor-community.csv`
- `Fungi_inoculation_experiment_*_resident-community.csv`
- `Fungi_inoculation_experiment_*_final-community.csv`
- `Bacteria_inoculation_experiment_*_donor-community.csv`
- `Bacteria_inoculation_experiment_*_resident-community.csv`
- `Bacteria_inoculation_experiment_*_final-community.csv`

Taxonomy tables may look like:

- `Fungi_inoculation_experiment_taxonomy_table.csv`
- `Bacteria_inoculation_experiment_taxonomy_table.csv`

Community matrix structure:
- Rows are samples.
- First column is usually `sample_id`, or should be treated as the sample identifier if unnamed or differently named.
- Remaining columns are taxon/ASV/feature IDs.
- Values are abundance or relative abundance.

Taxonomy table structure:
The taxonomy tables are NOT simple `id,name` tables. They likely have columns like:

- first column: feature/taxon ID, sometimes called `Unnamed: 0`
- `Kingdom`
- `Phylum`
- `Class`
- `Order`
- `Family`
- `Genus`
- `Species`

Important taxonomy rule:
Use `Genus` and `Species` to construct scientific names. Do NOT use `Kingdom` as the taxon name.

Construct:
- `name = "Genus Species"` when both Genus and Species are present.
- `name = "Genus"` when Species is missing but Genus is present.
- `name = ""` when Genus is missing.

No taxonomic assignments above genus should be used for traits or names.

Required scripts:
Create the following under `./workflows/code/`:

1. `build_graph_inputs.R`
Main workflow script. Running this script should generate all graph files.

2. Optionally, helper scripts such as:
- `graph_utils.R`
- `spieceasi_utils.R`

But keep the workflow easy to run from one top-level script.

The workflow should be runnable from the repository root with:


Rscript workflows/code/build_graph_inputs.R

Required outputs:
Write all outputs to ./workflows/output/graph_inputs/.

Generate these files:

combined_sample_taxon_edges.csv

A bipartite sample-taxon edge list with columns:

sample_id
taxon_id
abundance
kingdom
donor_id
community_type

Only include nonzero abundance edges.

nodes_samples.csv

Unique sample nodes with columns:

sample_id
kingdom
donor_id
community_type
nodes_taxa.csv

Unique taxon nodes with columns:

taxon_id
kingdom
name
genus
species
optionally: phylum, class, order, family

Use taxonomy tables to populate these fields.

taxon_taxon_spieceasi_edges.csv

Taxon-taxon co-occurrence networks built with SpiecEasi.

Run SpiecEasi separately for each stratum:

kingdom × donor_id × community_type

Do NOT pool bacteria with fungi.
Do NOT pool donor, resident, and final communities.
Do NOT pool different donor IDs.

The edge list should include:

kingdom
donor_id
community_type
taxon_a
taxon_b
weight
sign
method
n_samples
prevalence_a
prevalence_b

Use undirected edges only, with each taxon pair appearing once.

Default SpiecEasi method:
Use method = "glasso" by default because it gives a precision matrix that can be converted into partial-correlation-like edge weights. If the repo context or installed package constraints suggest otherwise, document the choice clearly.

For glasso:

Use the selected precision matrix from SpiecEasi.
Convert to partial correlations using:
P = -D^{-1/2} %*% Theta %*% D^{-1/2}

where D is the diagonal matrix from Theta.

Then:

weight = abs(partial correlation)
sign = sign(partial correlation)
method = "spieceasi_glasso"

Recommended configurable parameters:
Put these near the top of the script so workshop participants can change them:

min_prevalence <- 0.05
min_taxa <- 10
scale_to_counts <- TRUE
scale_factor <- 1e4
spieceasi_method <- "glasso"
nlambda <- 30
lambda_min_ratio <- 1e-2
rep_num <- 20
random_seed <- 1

Relative abundance handling:
If input values appear to be relative abundances and scale_to_counts = TRUE, multiply by scale_factor and round to pseudo-counts.

Zero handling:
Add a small pseudocount if needed for SpiecEasi stability.

Skipping:
If a stratum has too few taxa after prevalence filtering, skip it and record that skip in a summary file.

graph_edges_multirelational.csv

This is required, not optional.

Build this only after generating both:

combined_sample_taxon_edges.csv
taxon_taxon_spieceasi_edges.csv

It should combine bipartite sample-taxon edges and taxon-taxon co-occurrence edges into one flat edge list suitable for heterogeneous or multirelational GNN input.

Suggested columns:

kingdom
donor_id
community_type
from
to
weight
sign
relation
method

For bipartite edges:

from = sample_id
to = taxon_id
weight = abundance
sign = NA
relation = "sample_taxon"
method = "abundance"

For SpiecEasi co-occurrence edges:

from = taxon_a
to = taxon_b
weight = weight
sign = sign
relation = "taxon_taxon"
method = "spieceasi_glasso" or the selected method
spieceasi_run_summary.csv

A summary table with one row per stratum and columns like:

kingdom
donor_id
community_type
input_file
n_samples
n_taxa_original
n_taxa_after_filter
n_edges
status
message

Use status = "success", "skipped", or "failed".

README.md

Write this to ./workflows/output/graph_inputs/README.md.

The README should explain:

What the workflow does.
What input files were detected.
How filename metadata were parsed.
How sample-taxon bipartite edges were constructed.
How taxon names were constructed from taxonomy tables.
How SpiecEasi networks were constructed.
Why SpiecEasi was run separately by kingdom × donor_id × community_type.
What edge weights mean.
What sign means.
What each output file contains.
How to rerun the workflow.
How these files support GNN modeling.

Also include caveats:

Co-occurrence edges are inferred statistical associations, not measured direct interactions.
Positive co-occurrence does not necessarily mean cooperation.
Negative co-occurrence does not necessarily mean inhibition.
Associations can reflect shared niches, environmental filtering, compositional effects, or indirect interactions.
Prevalence filtering affects network density.
Relative-abundance to pseudo-count conversion is a modeling choice.
The bipartite backbone is the experimental data representation; the SpiecEasi graph is an inferred ecological association layer.

Implementation details:
Use R and tidyverse where helpful. Use SpiecEasi and Matrix. Add package checks at the top:

If SpiecEasi is missing, stop with installation guidance.
If tidyverse is missing, stop with installation guidance.
If Matrix is missing, stop with installation guidance.

Do not silently install packages.

Error handling:

If no community files are found, stop with a clear message.
If taxonomy tables are missing, still build sample-taxon edges but warn and leave taxon names/taxonomy blank.
If the first column of a community matrix is not named sample_id, rename it to sample_id and record a warning.
If SpiecEasi fails for a stratum, catch the error, record it in spieceasi_run_summary.csv, and continue with other strata.
If all SpiecEasi strata fail, still write the bipartite backbone and a multirelational file containing only sample-taxon edges, but make the README and summary clear.

Code quality:

Make the script modular with functions:
discover_input_files()
parse_community_filename()
read_community_matrix()
read_taxonomy_table()
build_taxon_nodes()
build_bipartite_edges()
prep_spieceasi_matrix()
run_spieceasi_for_stratum()
build_spieceasi_edges()
build_multirelational_edges()
write_readme()
Print useful progress messages.
Do not write outputs outside ./workflows/output/graph_inputs/.
Do not modify raw input files.
Do not modify website/reporting files.

After implementation:
Run or at least dry-run the workflow if possible. Then report:

files created
where outputs are written
any assumptions made
any skipped or failed strata
how to rerun the workflow
whether SpiecEasi was available in the environment

Make the pipeline workshop-friendly and reproducible. This is for a team that wants to use these graph files as inputs to a GNN project.

2026-05-13

Action

Added `workflows/code/build_graph_inputs.R`, a reproducible R workflow that discovers coalescence community/taxonomy inputs, validates structure, builds sample-taxon and multirelational graph edge outputs, runs stratum-specific SpiecEasi `glasso` networks with summary/error capture, and writes a workflow README under `workflows/output/graph_inputs/`. Performed syntax validation and a partial runtime check; full execution is expected to be slow in this environment because the community matrices are very wide.

2026-05-13

Handoff Note

Additional context from follow-up discussion:

- User has a previous SpiecEasi implementation using `method = "mb"`, `sel.criterion = "bstars"`, and `pulsar.params = list(rep.num = 20, ncores = parallel::detectCores() - 1)` on a relative-abundance `phyloseq` object.
- Current workflow in `workflows/code/build_graph_inputs.R` intentionally keeps `method = "glasso"` because the requested output is a weighted taxon-taxon edge file based on a precision matrix converted to partial-correlation-like weights.
- Comparison decision from this discussion: keep `glasso` for interpretability of weighted edges, but the user prefers to update the selection settings to use `sel.criterion = "bstars"` and `ncores = parallel::detectCores() - 1` in a later pass.
- User explicitly said not to worry much about runtime because they can run the workflow on HPC after the workshop.
- The workflow currently includes temporary runtime-guard choices added during local testing to keep the job bounded in this environment. These should be revisited in the next pass before HPC use, especially `max_taxa_for_spieceasi` and `spieceasi_time_limit_seconds`.
- Requested next step for the next Codex agent: preserve `glasso`, change SpiecEasi selection to `bstars`, add the user’s `ncores` setting, and reassess whether the temporary runtime guards should be relaxed or removed for HPC execution.

2026-05-13

Action

Updated `workflows/code/build_graph_inputs.R` for HPC-ready SpiecEasi defaults and lighter filtering pressure while preserving the existing graph-output structure. Changes: switched model selection to `sel.criterion = "bstars"`, added `pulsar.params$ncores = max(1, parallel::detectCores(logical = FALSE) - 1)`, set `max_taxa_for_spieceasi = Inf` so taxa are no longer capped by default, and set `spieceasi_time_limit_seconds = Inf` so runs are not artificially time-limited. Also updated the generated README parameter section to document the new `spieceasi_sel_criterion` and `spieceasi_ncores` settings.

2026-05-13

Action

Updated `workflows/code/build_graph_inputs.R` to use the user's original SpiecEasi defaults: `spieceasi_method = "mb"` and `spieceasi_sel_criterion = "bstars"`. Also widened method validation to allow both `mb` and `glasso` and made the taxon-taxon edge metadata method label dynamic (`spieceasi_<method>`).

2026-05-13

Action

Updated `.gitignore` to ignore generated workflow outputs under `workflows/output/**` to avoid pushing large local artifacts to GitHub. Added an exception for `workflows/output/.gitkeep` so the output directory can still be tracked when needed.

2026-05-14

Action

Recorded the active project goal in `workflows/README.md`: build a graph neural network workflow to predict community coalescence from donor and resident species co-occurrence/community matrices before communities are joined and compare predictions to observed final communities. Also cleaned up `workflows/README.md` formatting and expanded `workflows/code/README.md` with the current graph-input workflow and rerun command.

2026-05-14

Action

Added `workflows/code/train_coalescence_gnn.R`, an R `torch` workflow that reads `workflows/output/graph_inputs/` files, builds taxon-node prediction examples from `coalescence_triplets.csv` and `combined_sample_taxon_edges.csv`, trains separate bacteria/fungi multi-task models to predict final presence and abundance from donor/resident features, and writes model outputs under `workflows/output/gnn_model/`. The first-pass model intentionally ignores the empty SpiecEasi edge file and uses self-loop graph blocks plus learned taxon embeddings so it can run locally or on HPC before association edges are available.

2026-05-14

Action

Renamed the GNN training script to `workflows/code/train_coalescence_gnn_CMKR.R` so CMKR-created code files are clearly labeled. Updated the model to read `taxon_taxon_spieceasi_edges.csv` as a signed, weighted taxon-taxon association scaffold, build positive and negative normalized adjacency matrices, use those matrices in separate message-passing channels, and write `spieceasi_graph_summary.csv` plus `spieceasi_edges_used.csv` with model outputs.

2026-05-14

Action

Added explanatory comments to `workflows/code/train_coalescence_gnn_CMKR.R` without changing the executable code. Comments now describe package checks, configuration, input loading, taxon feature construction, SpiecEasi adjacency construction, data splitting, model architecture, training losses, evaluation metrics, and output files. Confirmed the script still parses with Rscript.

2026-05-14

Action

Added `workflows/code/visualize_gnn_performance_CMKR.R`, an R script that reads model outputs from `workflows/output/gnn_model`, visualizes kingdom-specific and combined GNN performance, writes plots and summary tables under `workflows/output/gnn_model/performance_visualizations_CMKR/`, and documents the generated visualization files. Confirmed the script parses with Rscript; local execution was not run because `data.table` and `ggplot2` are not installed in this Windows R environment.

2026-05-14

Action

Revised `workflows/code/visualize_gnn_performance_CMKR.R` to be fully separate from model training/compilation and fully tidyverse-based. The script now explicitly loads saved performance CSVs from `workflows/output/gnn_model`, uses `readr`, `dplyr`, `tidyr`, `purrr`, and `ggplot2` workflows rather than `data.table`, and still writes visualizations and tidy summaries under `workflows/output/gnn_model/performance_visualizations_CMKR/`. Confirmed the script parses with Rscript; local execution was not run because `tidyverse` is not installed in this Windows R environment.

2026-05-14

Action

Copied `workflows/code/train_coalescence_gnn_CMKR.R` to `workflows/code/train_coalescence_gnn_presence_absence_CMKR.R` for a presence/absence modeling variant. The copied script now binarizes all observed abundance values before feature and target construction, disables abundance loss by default, writes outputs to `workflows/output/gnn_model_presence_absence`, and replaces Bray-Curtis dissimilarity metrics with presence-set Jaccard and Sorensen similarity metrics for predicted versus observed final communities. Confirmed the copied script parses with Rscript.
