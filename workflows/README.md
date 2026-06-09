# Workflows

This directory contains reproducible workflows for preparing microbial coalescence data for GNN modeling.

- Raw input data live in `workflows/input/`.
- Analysis scripts live in `workflows/code/`.
- Generated outputs live in `workflows/output/`.

## Main Graph Workflow

Use this to transform existing coalescence community matrices into graph-ready files:

```bash
Rscript workflows/code/build_graph_inputs.R
```

## Human Gut Fecal-Transplant Workflow

Use this to download SRA amplicon data, process it with DADA2, assign taxonomy, build a phyloseq object, and export coalescence-style starting files for human gut test cases:

```bash
Rscript workflows/code/run_fecal_transplant_pipeline.R
```

Editable study inputs for that pipeline live under:

- `workflows/input/fecal_transplant_data/`

Generated human-gut community files are written to:

- `workflows/input/human_gut_fecal_transplant/`

Intermediate workflow outputs are written to:

- `workflows/output/fecal_transplant_pipeline/`



