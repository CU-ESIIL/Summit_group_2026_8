# Workflows

# Workflows

This directory contains reproducible workflows for preparing microbial coalescence data for GNN modeling.

- Raw input data live in `workflows/input/`.
- Analysis scripts live in `workflows/code/`.
- Generated outputs live in `workflows/output/`.

The main graph-preparation workflow is:

```bash
Rscript workflows/code/build_graph_inputs.R




