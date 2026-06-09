# Fecal Transplant Input Data

This directory holds the source metadata and editable mapping files for the human gut fecal-transplant test case.

- `SRA_metadata.csv` is the Sequence Read Archive run table used to download raw FASTQ files.
- `fecal_transplant_sample_inventory.csv` will be created by `workflows/code/00_fecal_transplant_setup.R` if it does not already exist.
- `fecal_transplant_triplets.csv` will be created by `workflows/code/00_fecal_transplant_setup.R` if it does not already exist.

The editable files are intentionally kept in `workflows/input/` so sample pairing decisions remain easy to review and update in GitHub.

Important:

- The current SRA metadata does not fully specify donor-to-recipient pairings or which recipient timepoints should be treated as `resident` versus `final`.
- Confirm those relationships from the associated publication before exporting training/test inputs with `05_fecal_transplant_export_coalescence_inputs.R`.
