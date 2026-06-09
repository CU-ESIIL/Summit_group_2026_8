# Workflow Code

This directory contains the graph-preparation workflow already used in this repository plus a new fecal-transplant preprocessing pipeline for building human gut test cases.

## Existing Graph Workflow

- `build_graph_inputs.R` builds graph-ready CSV files from coalescence community matrices.
- `prepare_gnn_dataset.py` converts those graph CSVs into a PyTorch Geometric dataset.

## New Fecal-Transplant Pipeline

The new pipeline mirrors the first five scripts in the reference `Sporophere_Bacteria/R` workflow:

1. `00_fecal_transplant_setup.R`
   Downloads taxonomy references, writes editable metadata templates, and downloads SRA FASTQ files.
2. `01_fecal_transplant_trim_primers.R`
   Removes primers if they have been configured, or stages raw reads unchanged until primer sequences are confirmed.
3. `02_fecal_transplant_build_asv_table.R`
   Runs DADA2 filtering, denoising, merging, chimera removal, and collapses technical runs into one ASV table per biological sample.
4. `03_fecal_transplant_assign_taxonomy.R`
   Assigns SILVA taxonomy to ASVs.
5. `04_fecal_transplant_build_phyloseq.R`
   Builds a phyloseq object from the ASV table, taxonomy, and editable sample inventory.
6. `05_fecal_transplant_export_coalescence_inputs.R`
   Exports coalescence-style `donor`, `resident`, and `final` CSVs plus a taxonomy table for downstream GNN testing.

You can run the full sequence with:

```bash
Rscript workflows/code/run_fecal_transplant_pipeline.R
```

## Important Editable Inputs

The setup script will create these if they do not already exist:

- `workflows/input/fecal_transplant_data/fecal_transplant_sample_inventory.csv`
- `workflows/input/fecal_transplant_data/fecal_transplant_triplets.csv`

The triplet file is where donor/resident/final pairings should be filled in after confirming them from the associated publication.

## IMC FMT Manifest Reconstruction

- `parse_imc_fmt_manifest.R`
  Reconstructs a conservative donor/resident/final manifest for the 2023 immune checkpoint colitis FMT study from the local SRA run table, supplement PDF, and workbook. It writes normalized inventories and a high-confidence GNN triple shortlist to `workflows/output/imc_fmt_2023/`.

## Meta-analysis Study Screening

- `screen_fmt_meta_analysis_studies.R`
  Reads the local supplementary workbook from the 2022 Nature Medicine FMT meta-analysis, applies a transparent screening rule for public-data usability, and writes ranked candidate-study tables to `workflows/output/fmt_meta_analysis_screen/`.

## Verma 2021 (PRJNA705895) Triads

- `parse_prjna705895_vermas2021_manifest.R`
  Parses `workflows/input/fecal_transplant_data/Verma_SRA_metadata.csv` plus the manuscript PDF text to build an unambiguous donor/pre/post triad manifest for `PRJNA705895`, writing outputs to `workflows/output/vermas_2021_prjna705895/`.

## Verma 2021 (Shotgun) Species Matrix

- `fecal_mg_processing_gzahn.sh`
  Metagenome processing pipeline for `PRJNA705895`: downloads SRR reads, runs QC, profiles species with MetaPhlAn4, and exports GNN-ready species abundance matrices plus coalescence-style donor/resident/final CSVs under `workflows/output/vermas_2021_prjna705895/mg_processing/`.
