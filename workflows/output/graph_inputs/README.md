# Graph Inputs Workflow

# Graph inputs for microbial coalescence GNN modeling

This directory contains graph-ready data products generated from the microbial coalescence preprocessing workflow.

These files are **outputs of the preprocessing workflow**, but they are intended to serve as **inputs to downstream graph neural network (GNN) modeling**.

The main preprocessing script is:

```bash
Rscript workflows/code/build_graph_inputs.R

-----------------------------------


Generated files are written to:

workflows/output/graph_inputs/

Raw input community matrices and taxonomy tables are expected in:

workflows/input/coalescence_experiments/
Biological context

The workshop dataset represents microbial coalescence experiments in which donor microbial communities were introduced into resident communities, producing final communities after assembly.

Conceptually:

resident community + donor source -> final community

For this dataset, donor communities are represented at the donor-source level. Each donor source was homogenized and applied to multiple resident communities. Therefore, all resident/final samples associated with a given donor_id share the same measured donor-source composition.

This means the donor input is treatment-level, not replicate-specific:

resident sample i + pooled donor source j -> final sample i

Resident and final samples are paired by shared sample_id within each kingdom × donor_id combination. These supervised prediction units are recorded in:

coalescence_triplets.csv
What the preprocessing script does

build_graph_inputs.R performs the following steps:

Discovers bacterial and fungal community matrices.
Reads donor, resident, and final community files.
Reads bacterial and fungal taxonomy tables.
Converts wide community matrices into long sample-taxon abundance edges.
Builds sample node metadata.
Builds taxon node metadata from taxonomy tables.
Constructs supervised resident-final prediction triplets linked to donor sources.
Runs SpiecEasi separately within each kingdom × donor_id × community_type stratum when sufficient variation exists.
Exports inferred taxon-taxon association edges.
Combines observed sample-taxon edges and inferred taxon-taxon edges into one multirelational edge file.
Writes a SpiecEasi run summary.
Input files

The workflow expects files like:

Bacteria_inoculation_experiment_0Burn_W1_donor-community.csv
Bacteria_inoculation_experiment_0Burn_W1_resident-community.csv
Bacteria_inoculation_experiment_0Burn_W1_final-community.csv
Fungi_inoculation_experiment_0Burn_W1_donor-community.csv
Fungi_inoculation_experiment_0Burn_W1_resident-community.csv
Fungi_inoculation_experiment_0Burn_W1_final-community.csv

and taxonomy tables:

Bacteria_inoculation_experiment_taxonomy_table.csv
Fungi_inoculation_experiment_taxonomy_table.csv

Community matrices are expected to have:

samples as rows,
the first column as sample_id or another sample identifier,
taxa/ASV/feature IDs as the remaining columns,
abundance or relative-abundance values in the matrix.

Taxonomy tables are expected to include columns such as:

Kingdom
Phylum
Class
Order
Family
Genus
Species

Taxon names are built from Genus and Species:

Genus Species when both are available,
Genus when only genus is available,
blank when genus is unavailable.

Higher taxonomic ranks are not used as substitute genus names.

Output files

All generated files are written to:

workflows/output/graph_inputs/
combined_sample_taxon_edges.csv

Observed bipartite sample-taxon edge list.

Columns:

sample_id
taxon_id
abundance
kingdom
donor_id
community_type

Each row represents a nonzero observed abundance:

sample_id -> taxon_id
weight = abundance

This file encodes the observed composition of donor, resident, and final communities.

These are empirical edges from the raw community matrices. They are not inferred associations.

nodes_samples.csv

Sample node metadata.

Columns:

sample_id
kingdom
donor_id
community_type

Each row describes a sample node and its role in the experiment.

Example:

R-020,Bacteria,0Burn_W1,donor
R-020,Bacteria,0Burn_W1,resident
R-020,Bacteria,0Burn_W1,final
F-R-020,Fungi,0Burn_W1,donor

For downstream GNN loading, it is safest to construct composite sample node IDs such as:

kingdom::donor_id::community_type::sample_id

rather than relying on raw sample_id alone.

This prevents accidental node collisions when the same sample ID appears in multiple contexts.

nodes_taxa.csv

Taxon node metadata.

Columns include:

kingdom
taxon_id
taxonomy_kingdom
phylum
class
order
family
genus
species
name

The name column is constructed as:

Genus Species when both genus and species are available,
Genus when only genus is available,
blank when genus is unavailable.

Unresolved taxa are retained as taxon nodes but are not assigned artificial genus-level names.

These taxon nodes can later be joined with trait tables using taxon_id, genus, or species, depending on the trait source.

coalescence_triplets.csv

Supervised prediction units for the workshop dataset.

Columns:

microcosm_id
kingdom
donor_id
donor_source_id
donor_sample_id
resident_sample_id
final_sample_id
donor_is_pooled
n_donor_samples
pairing_basis

Each row represents one prediction unit:

resident sample + pooled donor source -> final sample

Example:

microcosm_id = Bacteria::0Burn_W1::R-020
resident input = Bacteria / 0Burn_W1 / resident / R-020
donor input = pooled donor source 0Burn_W1
target = Bacteria / 0Burn_W1 / final / R-020

Because donor communities are pooled donor-source observations, donor_sample_id is NA and donor_is_pooled is TRUE.

The donor source is represented by donor_source_id, which matches donor_id.

taxon_taxon_spieceasi_edges.csv

Taxon-taxon association edges inferred with SpiecEasi.

Columns include:

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

SpiecEasi is run separately for each:

kingdom × donor_id × community_type

This prevents bacteria and fungi from being mixed and avoids pooling donor, resident, and final communities into one association network.

The weight column represents inferred association strength.

The sign column indicates whether the inferred association is positive or negative.

When method = "mb", edges should be interpreted as inferred neighborhood-selection associations. They are not direct ecological interactions and should not be treated as strict partial correlations.

Some donor-community strata may be skipped because pooled donor profiles lack sample-to-sample variation. These skips are expected. Donor composition is still represented through sample-taxon abundance edges.

graph_edges_multirelational.csv

Combined edge file for heterogeneous or multirelational GNN loading.

Columns:

kingdom
donor_id
community_type
from
to
weight
sign
relation
method

This file combines two edge types.

sample_taxon

Observed abundance edges:

from = sample_id
to = taxon_id
weight = abundance
sign = NA
relation = sample_taxon
method = abundance

These edges encode observed microbial community composition.

taxon_taxon

Inferred SpiecEasi association edges:

from = taxon_a
to = taxon_b
weight = inferred association strength
sign = positive or negative association
relation = taxon_taxon
method = SpiecEasi method used

These edges encode inferred taxon-taxon association structure.

This is the primary graph-ready edge list for a heterogeneous or multirelational GNN.

spieceasi_run_summary.csv

Summary of SpiecEasi runs, one row per attempted stratum.

Columns include:

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

The status column records whether each stratum:

succeeded,
was skipped,
or failed.

Skipped donor strata are expected when donor profiles are repeated or otherwise lack variable taxa.

How these files support GNN modeling

The graph representation contains two major relation types:

sample_taxon
taxon_taxon
Sample-taxon edges

Sample-taxon edges are observed abundance edges from the original community matrices.

They encode:

which taxa occur in each pooled donor source,
which taxa occur in each resident community,
which taxa occur in each final community,
and the abundance of each taxon in each sample.

These edges are the empirical backbone of the dataset.

They are especially important because donor communities are represented at the pooled donor-source level. Even when donor-specific SpiecEasi networks cannot be inferred, donor composition is still present through sample_taxon abundance edges.

Taxon-taxon edges

Taxon-taxon edges are inferred association edges from SpiecEasi.

They provide an ecological association scaffold among taxa.

These associations may reflect:

direct interactions,
indirect interactions,
shared niches,
environmental filtering,
compositional effects,
or other sources of covariance.

They should not be interpreted as confirmed biological interactions.

Triplets

coalescence_triplets.csv defines the supervised prediction units:

resident_sample + pooled_donor_source -> final_sample

This lets downstream modelers build training examples where donor and resident communities are inputs and the final community is the prediction target.

For this workshop dataset, the donor component should be interpreted as a treatment-level donor-source profile, not a replicate-specific donor profile.

Suggested downstream modeling use
Option 1: Heterogeneous graph

Use:

nodes_samples.csv
nodes_taxa.csv
graph_edges_multirelational.csv
coalescence_triplets.csv

The model has sample nodes and taxon nodes.

Sample-taxon edges encode observed composition.
Taxon-taxon edges encode inferred association structure.
Triplets define which resident/final samples belong to the same prediction unit and which donor source they received.

This approach keeps the generated graph files close to their current edge-list form.

Option 2: Taxon-node graph per prediction unit

Use coalescence_triplets.csv to build one graph per microcosm_id.

For each taxon, construct node features such as:

resident_abundance
donor_abundance
resident_present
donor_present
taxonomy or trait features

Then predict:

final_present
final_abundance

In this setup:

nodes are taxa,
edges are taxon-taxon SpiecEasi associations,
donor and resident communities are node features,
the final community is the target.

This option is closer to the preferred design for future experiments where each biological unit has directly measured donor, resident, and final communities.

Caveats
Donor communities in this workshop dataset are pooled donor-source profiles, not replicate-specific donor observations.
Co-occurrence edges are inferred statistical associations, not direct interaction measurements.
Positive associations do not necessarily indicate mutualism or facilitation.
Negative associations do not necessarily indicate competition or inhibition.
SpiecEasi networks require sample-to-sample variation; strata with repeated identical profiles are intentionally skipped.
Relative-abundance data and sparse compositional data require careful interpretation.
For strict predictive evaluation, avoid using final-community-derived association networks from held-out test samples.
Re-running the workflow

From the repository root:

Rscript workflows/code/build_graph_inputs.R

Outputs will be regenerated in:

workflows/output/graph_inputs/
Quick interpretation summary

The workflow produces graph-ready data for predicting microbial coalescence outcomes.

The most important files are:

coalescence_triplets.csv
combined_sample_taxon_edges.csv
taxon_taxon_spieceasi_edges.csv
graph_edges_multirelational.csv
nodes_samples.csv
nodes_taxa.csv

Together, these files describe:

resident sample + pooled donor source -> final sample

using:

observed sample-taxon abundance edges
+
inferred taxon-taxon association edges

The donor data are imperfect because donor communities were measured at the pooled source level, but they are still usable as treatment-level donor inputs for the workshop GNN.
