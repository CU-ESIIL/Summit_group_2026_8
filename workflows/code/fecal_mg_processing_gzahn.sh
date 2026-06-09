#!/usr/bin/env bash
set -euo pipefail

# fecal_mg_processing_gzahn.sh
#
# Repo-local metagenome processing pipeline for Verma 2021 (PRJNA705895).
# Goal: create species abundance matrices and coalescence-style donor/resident/final
# CSVs that can be ingested by the existing GNN workflows in this repository.
#
# This script DOES NOT install software. It only checks for required programs and
# runs them if present.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

TRIPLETS_CSV="${ROOT_DIR}/workflows/output/vermas_2021_prjna705895/prjna705895_gnn_triplets.csv"
SRA_METADATA_CSV="${ROOT_DIR}/workflows/input/fecal_transplant_data/Verma_SRA_metadata.csv"

OUT_DIR="${ROOT_DIR}/workflows/output/vermas_2021_prjna705895/mg_processing"
TMP_DIR="${OUT_DIR}/tmp"
FASTQ_DIR="${OUT_DIR}/fastq"
QC_DIR="${OUT_DIR}/qc_fastp"
LOG_DIR="${OUT_DIR}/logs"
DECONTAM_DIR="${OUT_DIR}/kneaddata"
METAPHLAN_DIR="${OUT_DIR}/metaphlan4"
MATRICES_DIR="${OUT_DIR}/matrices"
# Graph-input staging for the existing coalescence GNN pipeline.
GNN_INPUT_DIR="${ROOT_DIR}/workflows/input/verma2021_metagenome_coalescence"

THREADS="${THREADS:-8}"
# MetaPhlAn/Bowtie2 can be memory-hungry with the default SGB database.
# Use a smaller default nproc to reduce peak RAM; override as needed.
METAPHLAN_NPROC="${METAPHLAN_NPROC:-2}"

# Optional: set these in your shell if you want host-read removal.
KNEADDATA_DB="${KNEADDATA_DB:-${ROOT_DIR}/workflows/input/databases/hg38_bowtie2/hg38}"
METAPHLAN_DB="${METAPHLAN_DB:-}"

usage() {
  cat <<EOF
Usage:
  bash workflows/code/fecal_mg_processing_gzahn.sh [command]

Commands:
  check        Check required programs + inputs.
  manifest     Write a run-level manifest (samples + SRR accessions).
  download     Download SRR FASTQs for all runs in the triads.
  qc           Run fastp on downloaded FASTQs.
  decontam     (Optional) Run kneaddata host-read removal.
  metaphlan    Run MetaPhlAn4 species profiling.
  export       Build species abundance matrices + coalescence-style CSVs.
             Also stages per-recipient community matrices under workflows/input/
             so Rscript workflows/code/build_graph_inputs.R can discover them.
  all          Run: check -> manifest -> download -> qc -> metaphlan -> export

Environment:
  THREADS       Number of threads (default: 8)
  METAPHLAN_NPROC MetaPhlAn threads (default: 2)
  METAPHLAN_DB  MetaPhlAn4 database path (optional; uses tool defaults if empty)
  KNEADDATA_DB  KneadData Bowtie2 host db prefix (default: workflows/input/databases/hg38_bowtie2/hg38)

Notes:
  - This cohort is shotgun WGS (not 16S); DADA2 is not used.
  - Triads are derived from ${TRIPLETS_CSV}.
  - For Verma 2021, each recipient is treated as its own donor_id stratum
    (one sample per stratum). SpiecEasi will likely be skipped; use
    workflows/code/prepare_gnn_dataset.py --edge-scope none (or accept empty edges).
EOF
}

need_cmd() {
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Missing required program: ${cmd}" >&2
    return 1
  fi
}

mkdirs() {
  mkdir -p "${OUT_DIR}" "${TMP_DIR}" "${FASTQ_DIR}" "${QC_DIR}" "${LOG_DIR}" "${DECONTAM_DIR}" "${METAPHLAN_DIR}" "${MATRICES_DIR}"
}

check_inputs() {
  mkdirs
  test -f "${TRIPLETS_CSV}" || { echo "Missing triads file: ${TRIPLETS_CSV}" >&2; exit 2; }
  test -f "${SRA_METADATA_CSV}" || { echo "Missing SRA metadata file: ${SRA_METADATA_CSV}" >&2; exit 2; }

  need_cmd Rscript
  need_cmd awk
  need_cmd sed

  # Download tooling: pick one path and fail fast if neither exists.
  if command -v fasterq-dump >/dev/null 2>&1; then
    :
  elif command -v prefetch >/dev/null 2>&1 && command -v fastq-dump >/dev/null 2>&1; then
    :
  else
    echo "Need either: fasterq-dump (recommended) OR (prefetch + fastq-dump) from sra-tools." >&2
    exit 2
  fi

  # QC and profiling.
  need_cmd fastp
  need_cmd metaphlan

  # Optional tools.
  if [[ -n "${KNEADDATA_DB}" ]]; then
    need_cmd kneaddata
    need_cmd bowtie2
  fi

  echo "OK: inputs and required programs are present."
}

write_manifest() {
  mkdirs
  local out_tsv="${OUT_DIR}/run_manifest.tsv"

  # Manifest columns:
  # sample_id, role, recipient_id, run, library_name, collection_date
  #
  # sample_id format:
  #   donor_<recipient_id>[_r1/_r2 for multi-run donors]
  #   resident_<recipient_id>
  #   final_<recipient_id>
  #
  # For donor R with two runs, we keep two sample_ids: donor_R_1 and donor_R_2
  Rscript -e "
  tri <- read.csv('${TRIPLETS_CSV}', stringsAsFactors=FALSE)
  sra <- read.csv('${SRA_METADATA_CSV}', check.names=FALSE, stringsAsFactors=FALSE)
  names(sra) <- gsub(' ', '_', names(sra))
  sra\$Library_Name <- trimws(sra\$Library_Name)
  sra\$Sample_Name <- trimws(sra\$Sample_Name)
  sra\$Collection_Date <- trimws(sra\$Collection_Date)

  # index by Run
  idx <- setNames(seq_len(nrow(sra)), sra\$Run)

  rows <- list()
  add_row <- function(sample_id, role, recipient_id, run, library_name) {
    i <- idx[[run]]
    cd <- if (!is.null(i)) sra\$Collection_Date[[i]] else NA
    rows[[length(rows)+1]] <<- data.frame(
      sample_id=sample_id,
      role=role,
      recipient_id=recipient_id,
      run=run,
      library_name=library_name,
      collection_date=cd,
      stringsAsFactors=FALSE
    )
  }

  for (k in seq_len(nrow(tri))) {
    rid <- tri\$recipient_id[[k]]
    # donors can be semicolon-delimited
    donor_runs <- strsplit(tri\$donor_runs[[k]], ';', fixed=TRUE)[[1]]
    donor_runs <- donor_runs[nzchar(donor_runs)]
    if (length(donor_runs) == 1) {
      add_row(paste0('donor_', rid), 'donor', rid, donor_runs[[1]], tri\$donor_library[[k]])
    } else {
      for (j in seq_along(donor_runs)) {
        add_row(paste0('donor_', rid, '_', j), 'donor', rid, donor_runs[[j]], tri\$donor_library[[k]])
      }
    }
    add_row(paste0('resident_', rid), 'resident', rid, tri\$resident_run[[k]], tri\$resident_library[[k]])
    add_row(paste0('final_', rid), 'final', rid, tri\$final_run[[k]], tri\$final_library[[k]])
  }

  out <- do.call(rbind, rows)
  out <- out[order(out\$role, out\$recipient_id, out\$sample_id), ]
  write.table(out, file='${out_tsv}', sep='\\t', row.names=FALSE, quote=FALSE)
  cat('Wrote', nrow(out), 'rows to ${out_tsv}\\n')
  "

  echo "Manifest: ${out_tsv}"
}

download_fastqs() {
  mkdirs
  local manifest="${OUT_DIR}/run_manifest.tsv"
  test -f "${manifest}" || { echo "Missing manifest; run 'manifest' first." >&2; exit 2; }

  while IFS=$'\t' read -r sample_id role recipient_id run library_name collection_date; do
    [[ "${sample_id}" == "sample_id" ]] && continue
    local out1="${FASTQ_DIR}/${run}_1.fastq.gz"
    local out2="${FASTQ_DIR}/${run}_2.fastq.gz"
    if [[ -s "${out1}" && -s "${out2}" ]]; then
      continue
    fi

    echo "Downloading ${run} -> ${out1}, ${out2}"
    if command -v fasterq-dump >/dev/null 2>&1; then
      # fasterq-dump writes uncompressed fastq; compress afterwards.
      local tmp_run_dir="${TMP_DIR}/fasterq_${run}"
      mkdir -p "${tmp_run_dir}"
      fasterq-dump --split-files --threads "${THREADS}" --outdir "${tmp_run_dir}" "${run}"
      gzip -c "${tmp_run_dir}/${run}_1.fastq" > "${out1}"
      gzip -c "${tmp_run_dir}/${run}_2.fastq" > "${out2}"
      rm -f "${tmp_run_dir}/${run}_1.fastq" "${tmp_run_dir}/${run}_2.fastq"
      rmdir "${tmp_run_dir}" 2>/dev/null || true
    else
      # Fallback path for older sra-tools installs.
      prefetch "${run}"
      fastq-dump --split-files --gzip --outdir "${FASTQ_DIR}" "${run}"
      # fastq-dump names might be ${run}_1.fastq.gz, ${run}_2.fastq.gz already.
    fi
  done < "${manifest}"
}

run_fastp() {
  mkdirs
  local manifest="${OUT_DIR}/run_manifest.tsv"
  test -f "${manifest}" || { echo "Missing manifest; run 'manifest' first." >&2; exit 2; }

  while IFS=$'\t' read -r sample_id role recipient_id run library_name collection_date; do
    [[ "${sample_id}" == "sample_id" ]] && continue
    local in1="${FASTQ_DIR}/${run}_1.fastq.gz"
    local in2="${FASTQ_DIR}/${run}_2.fastq.gz"
    local out1="${QC_DIR}/${run}_1.fastq.gz"
    local out2="${QC_DIR}/${run}_2.fastq.gz"
    local html="${QC_DIR}/${run}.fastp.html"
    local json="${QC_DIR}/${run}.fastp.json"
    local log="${LOG_DIR}/${run}.fastp.log"
    if [[ -s "${out1}" && -s "${out2}" ]]; then
      continue
    fi
    test -s "${in1}" && test -s "${in2}" || { echo "Missing FASTQs for ${run}; run 'download' first." >&2; exit 2; }

    # Quick gzip integrity check to avoid wasting time on truncated downloads.
    if ! gzip -t "${in1}" >/dev/null 2>&1; then
      echo "Corrupt gzip: ${in1}" >&2
      exit 2
    fi
    if ! gzip -t "${in2}" >/dev/null 2>&1; then
      echo "Corrupt gzip: ${in2}" >&2
      exit 2
    fi

    {
      echo "[fastp] run=${run} sample_id=${sample_id} role=${role} recipient_id=${recipient_id}"
      echo "[fastp] in1=${in1}"
      echo "[fastp] in2=${in2}"
      echo "[fastp] out1=${out1}"
      echo "[fastp] out2=${out2}"
      echo "[fastp] threads=${THREADS}"
    } > "${log}"

    set +e
    fastp \
      --thread "${THREADS}" \
      --in1 "${in1}" --in2 "${in2}" \
      --out1 "${out1}" --out2 "${out2}" \
      --detect_adapter_for_pe \
      --html "${html}" --json "${json}" \
      --qualified_quality_phred 20 \
      --length_required 50 >> "${log}" 2>&1
    status=$?
    set -e

    if [[ ${status} -ne 0 ]]; then
      # Remove partial outputs so reruns won't "skip" truncated files.
      rm -f "${out1}" "${out2}" "${html}" "${json}"
      echo "fastp failed for ${run}; see ${log}" >&2
      exit ${status}
    fi
  done < "${manifest}"
}

run_kneaddata() {
  mkdirs
  if [[ -z "${KNEADDATA_DB}" ]]; then
    echo "KNEADDATA_DB is empty; skipping host-read removal." >&2
    echo "Set KNEADDATA_DB to your Bowtie2 human database path to enable." >&2
    exit 2
  fi

  local manifest="${OUT_DIR}/run_manifest.tsv"
  test -f "${manifest}" || { echo "Missing manifest; run 'manifest' first." >&2; exit 2; }

  while IFS=$'\t' read -r sample_id role recipient_id run library_name collection_date; do
    [[ "${sample_id}" == "sample_id" ]] && continue
    local in1="${QC_DIR}/${run}_1.fastq.gz"
    local in2="${QC_DIR}/${run}_2.fastq.gz"
    test -s "${in1}" && test -s "${in2}" || { echo "Missing fastp outputs for ${run}; run 'qc' first." >&2; exit 2; }

    local out_prefix="${DECONTAM_DIR}/${run}"
    local out1="${out_prefix}_paired_1.fastq"
    local out2="${out_prefix}_paired_2.fastq"
    if [[ -s "${out1}" && -s "${out2}" ]]; then
      continue
    fi

    kneaddata \
      --input "${in1}" --input "${in2}" \
      --output "${DECONTAM_DIR}" \
      --output-prefix "${run}" \
      --reference-db "${KNEADDATA_DB}" \
      --threads "${THREADS}" \
      --remove-intermediate-output
  done < "${manifest}"
}

run_metaphlan() {
  mkdirs
  local manifest="${OUT_DIR}/run_manifest.tsv"
  test -f "${manifest}" || { echo "Missing manifest; run 'manifest' first." >&2; exit 2; }

  while IFS=$'\t' read -r sample_id role recipient_id run library_name collection_date; do
    [[ "${sample_id}" == "sample_id" ]] && continue

    local in1=""
    local in2=""
    if [[ -n "${KNEADDATA_DB}" ]]; then
      local kd1="${DECONTAM_DIR}/${run}_paired_1.fastq"
      local kd2="${DECONTAM_DIR}/${run}_paired_2.fastq"
      if [[ -s "${kd1}" && -s "${kd2}" ]]; then
        in1="${kd1}"
        in2="${kd2}"
      fi
    fi
    if [[ -z "${in1}" ]]; then
      in1="${QC_DIR}/${run}_1.fastq.gz"
      in2="${QC_DIR}/${run}_2.fastq.gz"
    fi
    test -s "${in1}" && test -s "${in2}" || { echo "Missing inputs for MetaPhlAn for ${run}; run 'qc' (and optionally 'decontam') first." >&2; exit 2; }

    local out_profile="${METAPHLAN_DIR}/${run}.metaphlan_profile.tsv"
    local mapout="${METAPHLAN_DIR}/${run}.mapout.bz2"
    local log="${LOG_DIR}/${run}.metaphlan.log"
    if [[ -s "${out_profile}" ]]; then
      continue
    fi

    local db_args=()
    if [[ -n "${METAPHLAN_DB}" ]]; then
      db_args+=(--bowtie2db "${METAPHLAN_DB}")
    fi

    {
      echo "[metaphlan] run=${run} sample_id=${sample_id} role=${role} recipient_id=${recipient_id}"
      echo "[metaphlan] in1=${in1}"
      echo "[metaphlan] in2=${in2}"
      echo "[metaphlan] out_profile=${out_profile}"
      echo "[metaphlan] nproc=${METAPHLAN_NPROC}"
    } > "${log}"

    set +e
    metaphlan \
      "${in1},${in2}" \
      --input_type fastq \
      --nproc "${METAPHLAN_NPROC}" \
      "${db_args[@]}" \
      --mapout "${mapout}" \
      --no_map \
      --tax_lev s \
      -o "${out_profile}" >> "${log}" 2>&1
    status=$?
    set -e

    if [[ ${status} -ne 0 ]]; then
      rm -f "${out_profile}" "${mapout}"
      if [[ ${status} -eq 137 ]]; then
        echo "metaphlan failed for ${run} (exit 137). This is commonly an OOM kill; try lowering METAPHLAN_NPROC or run on a high-memory node." >&2
      fi
      echo "metaphlan failed for ${run}; see ${log}" >&2
      exit ${status}
    fi
  done < "${manifest}"
}

export_matrices() {
  mkdirs
  local manifest="${OUT_DIR}/run_manifest.tsv"
  test -f "${manifest}" || { echo "Missing manifest; run 'manifest' first." >&2; exit 2; }

  # Collect list of per-run MetaPhlAn profiles.
  local profiles_list="${TMP_DIR}/metaphlan_profiles.txt"
  awk -F'\t' 'NR>1 {print $4}' "${manifest}" | sort -u | awk -v d="${METAPHLAN_DIR}" '{print d"/"$1".metaphlan_profile.tsv"}' > "${profiles_list}"

  # Build merged + pooled species abundance matrices (samples x species).
  #
  # - The merged matrix is run/sample-level output (includes donor replicate runs).
  # - The pooled matrices are "coalescence-style" and are shaped to be consumable
  #   by the existing graph/GNN tooling in this repo:
  #     * community_type matrices share the same sample_id per triad
  #     * donor_id is encoded in the filename (one donor per recipient)
  #     * donor replicates are averaged into one pooled donor profile
  local species_matrix="${MATRICES_DIR}/species_abundance_matrix.csv"
  local donor_mat="${MATRICES_DIR}/Verma_2021_donor-community.csv"
  local resident_mat="${MATRICES_DIR}/Verma_2021_resident-community.csv"
  local final_mat="${MATRICES_DIR}/Verma_2021_final-community.csv"
  local taxonomy="${MATRICES_DIR}/Verma_2021_taxonomy_table.csv"
  local gnn_taxonomy="${GNN_INPUT_DIR}/Bacteria_inoculation_experiment_taxonomy_table.csv"

  Rscript -e "
  tri <- read.csv('${TRIPLETS_CSV}', stringsAsFactors=FALSE)
  man <- read.delim('${manifest}', stringsAsFactors=FALSE)

  # Read MetaPhlAn outputs (species level already requested via --tax_lev s).
  read_profile <- function(path) {
    x <- read.delim(path, sep='\\t', header=TRUE, comment.char='', quote='', stringsAsFactors=FALSE)
    # MetaPhlAn outputs 'clade_name' + 'relative_abundance' (and potentially more columns).
    cn <- names(x)
    clade_col <- cn[1]
    ra_col <- if ('relative_abundance' %in% cn) 'relative_abundance' else cn[length(cn)]
    x <- x[, c(clade_col, ra_col)]
    names(x) <- c('taxon', 'abundance')
    x
  }

  all_runs <- sort(unique(man\$run))
  taxa_union <- character()
  profiles <- list()

  for (r in all_runs) {
    p <- file.path('${METAPHLAN_DIR}', paste0(r, '.metaphlan_profile.tsv'))
    if (!file.exists(p)) stop('Missing MetaPhlAn profile: ', p)
    prof <- read_profile(p)
    # Keep only species strings (MetaPhlAn s__ clades). In practice these should already be species.
    prof <- prof[grepl('^s__', prof\$taxon), , drop=FALSE]
    profiles[[r]] <- prof
    taxa_union <- union(taxa_union, prof\$taxon)
  }

  taxa_union <- sort(taxa_union)

  # Build sample x taxa matrix (run-level, duplicating rows for donor replicate sample_ids).
  run_to_sample_ids <- split(man\$sample_id, man\$run)
  mat_rows <- list()
  for (r in all_runs) {
    prof <- profiles[[r]]
    v <- setNames(rep(0, length(taxa_union)), taxa_union)
    v[prof\$taxon] <- prof\$abundance
    # Some runs map to multiple sample_ids (donor replicates); duplicate rows.
    sids <- run_to_sample_ids[[r]]
    for (sid in sids) {
      mat_rows[[length(mat_rows)+1]] <- c(sample_id=sid, v)
    }
  }

  mat <- do.call(rbind, mat_rows)
  mat <- as.data.frame(mat, stringsAsFactors=FALSE)
  # Coerce numeric columns.
  for (j in 2:ncol(mat)) mat[[j]] <- as.numeric(mat[[j]])

  # Write full matrix.
  write.csv(mat, file='${species_matrix}', row.names=FALSE, quote=FALSE)

  # Build pooled donor/resident/final matrices with shared sample_id per triad.
  # sample_id is recipient_id (A, B, ...), consistent across all three matrices.
  pool_rows <- function(run_ids, sample_id) {
    # run_ids may contain NA/empty.
    run_ids <- run_ids[!is.na(run_ids) & nzchar(run_ids)]
    if (length(run_ids) < 1) return(NULL)
    # Average across runs (needed for donor replicates).
    m <- mat[mat\$sample_id %in% unlist(lapply(run_ids, function(r) run_to_sample_ids[[r]])), , drop=FALSE]
    if (nrow(m) < 1) return(NULL)
    vals <- colMeans(m[, -1, drop=FALSE])
    c(sample_id=sample_id, vals)
  }

  donor_rows <- list()
  resident_rows <- list()
  final_rows <- list()

  for (k in seq_len(nrow(tri))) {
    rid <- tri\$recipient_id[[k]]
    donor_runs <- strsplit(tri\$donor_runs[[k]], ';', fixed=TRUE)[[1]]
    donor_runs <- donor_runs[nzchar(donor_runs)]
    resident_run <- tri\$resident_run[[k]]
    final_run <- tri\$final_run[[k]]

    dr <- pool_rows(donor_runs, rid)
    rr <- pool_rows(c(resident_run), rid)
    fr <- pool_rows(c(final_run), rid)
    if (!is.null(dr)) donor_rows[[length(donor_rows)+1]] <- dr
    if (!is.null(rr)) resident_rows[[length(resident_rows)+1]] <- rr
    if (!is.null(fr)) final_rows[[length(final_rows)+1]] <- fr
  }

  donor_mat_df <- as.data.frame(do.call(rbind, donor_rows), stringsAsFactors=FALSE)
  resident_mat_df <- as.data.frame(do.call(rbind, resident_rows), stringsAsFactors=FALSE)
  final_mat_df <- as.data.frame(do.call(rbind, final_rows), stringsAsFactors=FALSE)
  for (j in 2:ncol(donor_mat_df)) donor_mat_df[[j]] <- as.numeric(donor_mat_df[[j]])
  for (j in 2:ncol(resident_mat_df)) resident_mat_df[[j]] <- as.numeric(resident_mat_df[[j]])
  for (j in 2:ncol(final_mat_df)) final_mat_df[[j]] <- as.numeric(final_mat_df[[j]])

  write.csv(donor_mat_df, file='${donor_mat}', row.names=FALSE, quote=FALSE)
  write.csv(resident_mat_df, file='${resident_mat}', row.names=FALSE, quote=FALSE)
  write.csv(final_mat_df, file='${final_mat}', row.names=FALSE, quote=FALSE)

  # Minimal taxonomy table for the downstream workflow (Genus/Species extracted from s__ labels).
  # taxon_id is the MetaPhlAn species string (s__Genus_species).
  tax <- data.frame(taxon_id=taxa_union, stringsAsFactors=FALSE)
  # Strip s__ prefix and convert underscores to spaces.
  sci <- gsub('^s__', '', tax\$taxon_id)
  sci <- gsub('_', ' ', sci)
  parts <- strsplit(sci, ' ', fixed=TRUE)
  tax\$Genus <- vapply(parts, function(p) if (length(p) >= 1) p[[1]] else '', character(1))
  tax\$Species <- vapply(parts, function(p) if (length(p) >= 2) paste(p[-1], collapse=' ') else '', character(1))
  tax\$name <- ifelse(nchar(tax\$Species) > 0, paste(tax\$Genus, tax\$Species), tax\$Genus)
  write.csv(tax, file='${taxonomy}', row.names=FALSE, quote=FALSE)

  # Also write taxonomy in the format expected by build_graph_inputs.R (first col is id + Genus/Species columns).
  write.csv(tax[, c('taxon_id','Genus','Species','name')], file='${gnn_taxonomy}', row.names=FALSE, quote=FALSE)

  cat('Wrote species matrix: ${species_matrix}\\n')
  cat('Wrote donor/resident/final matrices under: ${MATRICES_DIR}\\n')
  "

  # Stage community matrices for build_graph_inputs.R discovery.
  # It searches for: Bacteria_inoculation_experiment_<donor_id>_(donor|resident|final)-community.csv
  mkdir -p "${GNN_INPUT_DIR}"
  # We create one donor_id stratum per recipient_id (A..V) so the existing triplet builder works.
  Rscript -e "
  tri <- read.csv('${TRIPLETS_CSV}', stringsAsFactors=FALSE)
  donor <- read.csv('${donor_mat}', stringsAsFactors=FALSE, check.names=FALSE)
  resident <- read.csv('${resident_mat}', stringsAsFactors=FALSE, check.names=FALSE)
  final <- read.csv('${final_mat}', stringsAsFactors=FALSE, check.names=FALSE)

  out_dir <- '${GNN_INPUT_DIR}'
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  write_one <- function(df, rid, community_type) {
    row <- df[df\$sample_id == rid, , drop=FALSE]
    if (nrow(row) != 1) stop('Expected 1 row for ', rid, ' in ', community_type)
    # Rename sample_id to match build_graph_inputs expectations.
    names(row)[1] <- 'sample_id'
    fn <- file.path(out_dir, paste0('Bacteria_inoculation_experiment_Verma2021_', rid, '_', community_type, '-community.csv'))
    write.csv(row, fn, row.names=FALSE, quote=FALSE)
  }

  for (rid in tri\$recipient_id) {
    write_one(donor, rid, 'donor')
    write_one(resident, rid, 'resident')
    write_one(final, rid, 'final')
  }
  cat('Staged coalescence community matrices under:', out_dir, '\\n')
  "

  echo "Exported matrices to: ${MATRICES_DIR}"
  echo "Staged GNN input community matrices to: ${GNN_INPUT_DIR}"
}

cmd="${1:-}"
case "${cmd}" in
  check) check_inputs ;;
  manifest) check_inputs; write_manifest ;;
  download) check_inputs; write_manifest; download_fastqs ;;
  qc) check_inputs; write_manifest; run_fastp ;;
  decontam) check_inputs; write_manifest; run_fastp; run_kneaddata ;;
  metaphlan) check_inputs; write_manifest; run_fastp; run_metaphlan ;;
  export) check_inputs; write_manifest; export_matrices ;;
  all) check_inputs; write_manifest; download_fastqs; run_fastp; run_metaphlan; export_matrices ;;
  ""|-h|--help|help) usage ;;
  *) echo "Unknown command: ${cmd}" >&2; usage; exit 2 ;;
esac
