#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
bioproject_acc <- "PRJNA705895"

sra_table_path <- file.path(
  repo_root, "workflows", "input", "fecal_transplant_data", "Verma_SRA_metadata.csv"
)
manuscript_txt_path <- "/tmp/journal.pone.0251590.txt"

output_dir <- file.path(repo_root, "workflows", "output", "vermas_2021_prjna705895")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(sra_table_path)) {
  stop("Missing SRA table: ", sra_table_path, call. = FALSE)
}

parse_months <- function(library_name) {
  m <- str_match(library_name, "\\.(\\d+(?:\\.\\d+)?)\\.month$")
  as.numeric(m[, 2])
}

read_sra_table <- function(path) {
  # Use base R read.csv to avoid extra deps and preserve headers exactly.
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  as_tibble(x)
}

df <- read_sra_table(sra_table_path) %>%
  mutate(
    library_name = str_trim(.data[["Library Name"]]),
    run = .data[["Run"]],
    bioproject = .data[["BioProject"]],
    sample_name = .data[["Sample Name"]],
    collection_date = suppressWarnings(as.Date(.data[["Collection_Date"]])),
    role = case_when(
      str_detect(library_name, "^Donor\\.") ~ "donor",
      str_detect(library_name, "^PreFMT\\.") ~ "resident_pre",
      str_detect(library_name, "^PostFMT\\.") ~ "recipient_post",
      TRUE ~ "other"
    ),
    recipient_id = case_when(
      role == "resident_pre" ~ str_match(library_name, "^PreFMT\\.([^\\.]+)$")[, 2],
      role == "recipient_post" ~ str_match(library_name, "^PostFMT\\.([^\\.]+)\\..*$")[, 2],
      TRUE ~ NA_character_
    ),
    donor_id = case_when(
      role == "donor" ~ str_match(library_name, "^Donor\\.([^\\.]+)")[, 2],
      TRUE ~ NA_character_
    ),
    months = parse_months(library_name)
  )

if (!all(unique(df$bioproject) == bioproject_acc)) {
  stop("BioProject mismatch in SRA table; expected ", bioproject_acc, call. = FALSE)
}

write_csv(df, file.path(output_dir, "prjna705895_sra_table_normalized.csv"))

pre_samples <- df %>% filter(role == "resident_pre")
post_samples <- df %>% filter(role == "recipient_post")
donor_samples <- df %>% filter(role == "donor")

recipient_ids <- pre_samples %>%
  distinct(recipient_id) %>%
  filter(!is.na(recipient_id)) %>%
  arrange(recipient_id) %>%
  pull(recipient_id)

if (length(recipient_ids) == 0) {
  stop("No PreFMT.* rows found; cannot construct triads.", call. = FALSE)
}

pick_final_post <- function(posts_for_recipient) {
  # Manuscript states post-FMT samples used for strain-engraftment analysis were
  # collected ~30 days or later. In the SRA table, those are the ".1.month" (or later)
  # PostFMT timepoints. When a >=1-month sample is present, use the earliest such
  # timepoint; otherwise fall back to the earliest observed post sample.
  posts_for_recipient <- posts_for_recipient %>% mutate(months = parse_months(library_name))

  has_ge_1m <- any(!is.na(posts_for_recipient$months) & posts_for_recipient$months >= 1)
  if (has_ge_1m) {
    posts_for_recipient <- posts_for_recipient %>% filter(!is.na(months) & months >= 1)
  }

  posts_for_recipient %>%
    arrange(months, collection_date, library_name, run) %>%
    slice(1)
}

pick_donor <- function(donors_for_id) {
  # Most donors appear once (Donor.X). For K and O, multiple longitudinal donor stools
  # are present; use the donor sample at 0 months when available (same-day as PreFMT).
  if (nrow(donors_for_id) == 0) return(NULL)

  donors_for_id <- donors_for_id %>% mutate(months = parse_months(library_name))

  if (any(!is.na(donors_for_id$months) & donors_for_id$months == 0)) {
    donors_for_id <- donors_for_id %>% filter(!is.na(months) & months == 0)
  }

  if (nrow(donors_for_id) == 1) return(donors_for_id)

  # Donor.R has replicate1/2; keep both runs.
  donors_for_id %>% arrange(library_name, run)
}

build_one_triads_row <- function(rid) {
  pre <- pre_samples %>% filter(recipient_id == rid)
  posts <- post_samples %>% filter(recipient_id == rid)

  donor_key <- rid
  if (rid %in% c("T", "U")) donor_key <- "T/U"
  donors <- donor_samples %>% filter(donor_id == donor_key)

  if (nrow(pre) != 1) {
    return(tibble(
      recipient_id = rid,
      donor_library = NA_character_,
      donor_runs = NA_character_,
      resident_library = NA_character_,
      resident_run = NA_character_,
      final_library = NA_character_,
      final_run = NA_character_,
      note = sprintf("Expected exactly 1 PreFMT sample; found %d", nrow(pre))
    ))
  }

  if (nrow(posts) < 1) {
    return(tibble(
      recipient_id = rid,
      donor_library = NA_character_,
      donor_runs = NA_character_,
      resident_library = pre$library_name,
      resident_run = pre$run,
      final_library = NA_character_,
      final_run = NA_character_,
      note = "No PostFMT samples found"
    ))
  }

  final <- pick_final_post(posts)
  donor_pick <- pick_donor(donors)

  if (is.null(donor_pick) || nrow(donor_pick) < 1) {
    return(tibble(
      recipient_id = rid,
      donor_library = NA_character_,
      donor_runs = NA_character_,
      resident_library = pre$library_name,
      resident_run = pre$run,
      final_library = final$library_name,
      final_run = final$run,
      note = sprintf("No Donor sample found for donor_id=%s", donor_key)
    ))
  }

  tibble(
    recipient_id = rid,
    donor_library = paste(unique(donor_pick$library_name), collapse = ";"),
    donor_runs = paste(unique(donor_pick$run), collapse = ";"),
    resident_library = pre$library_name,
    resident_run = pre$run,
    final_library = final$library_name,
    final_run = final$run,
    note = NA_character_
  )
}

triads <- lapply(recipient_ids, build_one_triads_row) %>%
  bind_rows() %>%
  arrange(recipient_id)

write_csv(triads, file.path(output_dir, "prjna705895_gnn_triplets.csv"))

issues <- triads %>% filter(!is.na(note))
write_csv(issues, file.path(output_dir, "prjna705895_gnn_triplets_issues.csv"))

if (file.exists(manuscript_txt_path)) {
  # Save a short evidence snippet used for the post-FMT timing rule.
  txt <- readLines(manuscript_txt_path, warn = FALSE)
  idx <- grep("30\\s+days\\s+or\\s+later", txt, ignore.case = TRUE)
  if (length(idx) > 0) {
    lo <- max(1, idx[1] - 2)
    hi <- min(length(txt), idx[1] + 6)
    writeLines(txt[lo:hi], file.path(output_dir, "manuscript_evidence_postfmt_timing.txt"))
  }
}

message("Wrote outputs to ", output_dir)
