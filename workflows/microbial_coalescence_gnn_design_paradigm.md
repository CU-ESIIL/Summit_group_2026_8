# Design paradigm for GNN modeling of microbial coalescence experiments

## Purpose

This document describes a graph-based modeling strategy for microbial coalescence experiments in which donor and resident microbial communities combine to produce a final community. It is written as guidance for a future Codex agent or data-science collaborator building a reproducible workflow.

For the current workshop dataset, donor communities are relatively blunt inputs because donor profiles are repeated across samples and donor-specific co-occurrence networks cannot be inferred. The workflow can still use donor composition through sample-taxon or taxon-level abundance features.

For the future seagrass inoculation experiment, each biological unit (one microcosm) will ideally have a clearly matched donor community, resident community, and final community:

```text
donor_i + resident_i -> final_i
```

That matched triplet design is the preferred structure for accurate predictive modeling.

---

## Core biological prediction task

The central supervised learning task is:

```text
Given donor_i and resident_i, predict final_i
```

where each community is represented as a sparse taxon abundance vector.

More formally:

```text
final_abundance_vector_i = f(resident_abundance_vector_i, donor_abundance_vector_i, taxon traits, taxon-taxon structure, covariates)
```

The model should learn which taxa persist, invade, disappear, or change abundance after coalescence.

---

## Recommended graph representation for high-accuracy future work

For the future seagrass inoculation experiment, use one graph instance per biological unit or microcosm.

### Graph instance

```text
Graph_i = one microcosm
Input = matched resident_i and donor_i communities
Target = final_i community
```

### Nodes

Nodes are taxa.

Use a consistent taxon universe across graph instances, at least within kingdom. Depending on the final modeling choice, bacteria and fungi can be modeled:

1. separately, with one model per kingdom, or
2. jointly, with kingdom as a taxon feature and possibly kingdom-specific edge types.

### Node features

For each taxon node in microcosm `i`, include features such as:

```text
resident_abundance_i
donor_abundance_i
resident_present_i
donor_present_i
resident_log_abundance_i
donor_log_abundance_i
resident_rank_i
donor_rank_i
taxon traits
taxonomic features
optional sequence or phylogenetic features
optional environmental or microcosm-level covariates broadcast to nodes
```

At minimum, include:

```text
resident_abundance_i
donor_abundance_i
resident_present_i
donor_present_i
```

The donor community does not require its own inferred donor-specific co-occurrence graph. The donor enters the model through taxon-level abundance and presence features.

### Edges

Edges are taxon-taxon relationships. They provide relational structure for message passing.

Potential edge types include:

```text
taxon_taxon_spieceasi
taxon_taxon_sequence_similarity
taxon_taxon_trait_similarity
taxon_taxon_known_interaction
```

For a first serious ecological GNN, prioritize:

```text
taxon_taxon_spieceasi
```

and optionally add:

```text
taxon_taxon_sequence_similarity
```

if sequence data are available and biologically meaningful.

### Edge attributes

For SpiecEasi-derived taxon-taxon edges, include:

```text
edge_weight
edge_sign
edge_method
edge_context
```

Example:

```text
from_taxon, to_taxon, edge_weight, edge_sign, edge_type, edge_context
ASV_1, ASV_7, 0.42, -1, spieceasi_mb, resident_training
```

Where:

- `edge_weight` is the magnitude of the inferred association.
- `edge_sign` is positive or negative.
- `edge_type` records the edge source.
- `edge_context` records where the association network came from, for example resident communities, donor communities, or training-only baseline samples.

---

## Important distinction: donor input vs donor network

A donor community can be highly informative even if no donor-specific co-occurrence network is inferred.

A single donor profile is an abundance vector. A co-occurrence network requires variation across multiple samples. Therefore, unless the experiment includes replicate subsamples within each donor community, a donor-specific SpiecEasi graph for each microcosm is not estimable.

This is not a problem.

The donor community should be encoded as node features:

```text
donor_abundance_i
donor_present_i
```

The co-occurrence graph acts as a background ecological association scaffold over taxa. The donor vector tells the model which taxa were introduced and at what abundances.

---

## Why this is better than separate donor and resident graphs

A tempting structure is:

```text
Graph(resident_i), Graph(donor_i) -> final_i
```

This is conceptually intuitive but not ideal unless each resident and donor graph has enough replicate variation to infer its own network.

A stronger structure is:

```text
One taxon graph per microcosm
Node features = donor and resident abundance information
Edges = inferred taxon-taxon association structure
Target = final abundance or presence
```

This lets the model learn rules like:

```text
If taxon A is abundant in the donor,
and taxon B is abundant in the resident,
and A-B have a negative inferred association,
and A and B differ in relevant traits,
then A may be less likely to establish in the final community.
```

That is the biological pattern the model should be able to discover.

---

## Recommended target structure

Predict two related outcomes per taxon:

```text
final_presence_i
final_abundance_i
```

A multi-task setup is recommended:

1. Classification head: predicts whether each taxon is present in the final community.
2. Regression head: predicts final abundance, ideally conditional on presence or with a zero-inflated/compositional loss.

This is usually better than forcing a single regression model to predict sparse compositional data.

Possible target columns per graph instance:

```text
taxon_id
final_present
final_abundance
final_log_abundance
```

---

## Data leakage caution

Be careful when using final-community data to infer SpiecEasi networks.

If the goal is strict prediction, do not infer taxon-taxon edges from all final communities and then evaluate on held-out final communities. That leaks information from the test set into the graph structure.

Safer options:

1. Infer association networks using only training samples within each cross-validation fold.
2. Use only pre-coalescence resident and donor communities to infer the association scaffold.
3. Use external or independently generated association/interaction networks.
4. If final-derived edges are used, document clearly that they are transductive and not a strict out-of-sample prediction setup.

For future high-accuracy work, the preferred workflow is:

```text
Split microcosms into train/validation/test.
Infer any learned association networks using training data only.
Apply the learned graph scaffold to validation/test graph instances.
```

---

## Current workshop dataset strategy

The current workshop dataset is less clean than the future seagrass design because donor communities may be repeated identical profiles. This means donor-specific SpiecEasi networks may be skipped due to zero sample-to-sample variation.

That is expected and acceptable.

For the workshop dataset:

1. Build the bipartite sample-taxon backbone:

```text
sample -> taxon, weight = abundance
```

2. Build SpiecEasi taxon-taxon networks only for strata with enough variation:

```text
kingdom x donor_id x community_type
```

3. Record skipped donor strata explicitly, not as failures.

4. Create a multirelational edge file:

```text
sample_taxon edges
taxon_taxon edges
```

5. Use donor information through observed abundance edges or through derived donor abundance node features.

The workshop graph files are useful for prototyping, but the future matched microcosm design should move toward graph-per-microcosm objects.

---

## Future seagrass inoculation experiment: preferred data structure

For each microcosm, collect and store:

```text
microcosm_id
treatment_id
block_id
site_id
date
donor_sample_id
resident_sample_id
final_sample_id
environmental_covariates
```

And three matched community profiles:

```text
donor community actually delivered
resident community at inoculation or transplant
final rhizoplane/rhizosphere community at harvest
```

This matched design is critical. Avoid relying on treatment-level averaged donors if possible.

---

## Recommended processed files for future GNN workflow

### 1. Taxon table

```text
taxon_id
kingdom
name
genus
species
phylum
class
order
family
trait columns...
```

### 2. Microcosm metadata

```text
microcosm_id
donor_sample_id
resident_sample_id
final_sample_id
treatment_id
block_id
site_id
environmental covariates...
```

### 3. Long abundance table

```text
sample_id
microcosm_id
sample_role    # donor, resident, final
taxon_id
abundance
presence
```

### 4. Graph node feature table

One row per microcosm x taxon:

```text
microcosm_id
taxon_id
resident_abundance
donor_abundance
resident_present
donor_present
trait features...
final_abundance
final_present
```

### 5. Taxon-taxon edge table

```text
from_taxon
to_taxon
edge_weight
edge_sign
edge_type
edge_context
```

### 6. Model-ready graph objects

For PyTorch Geometric or DGL, each microcosm should become a graph object:

```text
x = taxon node feature matrix
edge_index = taxon-taxon edges
edge_attr = edge weights, signs, edge types
y_presence = final presence vector
y_abundance = final abundance vector
metadata = microcosm_id, treatment, block, etc.
```

---

## Suggested model families

Start simple, then add complexity.

### Baselines

Always include non-GNN baselines:

```text
final = resident only
final = donor + resident linear model
random forest / XGBoost using donor and resident abundance features
MLP using donor and resident abundance features
```

If the GNN does not beat these, the graph is not helping.

### First GNN

Taxon-node graph per microcosm:

```text
GraphSAGE or GAT
node features = donor/resident abundances + traits
edges = SpiecEasi associations
targets = final presence and abundance
```

### More advanced

Heterogeneous/multirelational GNN with edge types:

```text
spieceasi_positive
spieceasi_negative
sequence_similarity
trait_similarity
```

Potentially separate bacteria and fungi first, then combine once the simpler setup behaves.

---

## Biological interpretation goal

The model should not only predict final communities. It should eventually help answer:

```text
Which donor taxa establish?
Which resident taxa resist invasion?
Which donor-resident trait pairings predict successful coalescence?
Which inferred taxon-taxon associations matter most for final community assembly?
```

Important interpretability targets:

```text
donor abundance effects
resident abundance effects
taxon traits
trait-pair interactions
edge sign and edge weight
community context
```

---

## Key caveats

- SpiecEasi edges are inferred statistical associations, not measured direct interactions.
- Positive associations may reflect shared habitat preferences, facilitation, or indirect effects.
- Negative associations may reflect competition, divergent niches, compositional constraints, or indirect effects.
- Donor communities can be highly useful even without donor-specific inferred networks.
- Do not infer graph structure from held-out final communities when evaluating strict prediction.
- Matched donor-resident-final sampling is the strongest experimental design for this modeling framework.

---

## One-sentence summary

Use donor and resident communities as taxon-level input features, use taxon-taxon edges as ecological association structure, and predict final community composition for each matched microcosm.
