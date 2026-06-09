#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
})

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
input_xlsx <- "/home/geoff/Downloads/41591_2022_1964_MOESM3_ESM.xlsx"
output_dir <- file.path(repo_root, "workflows", "output", "fmt_meta_analysis_screen")

if (!file.exists(input_xlsx)) {
  stop("Meta-analysis workbook not found: ", input_xlsx, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

parse_triads <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    str_match(x, "^(\\d+)\\s*/\\s*(\\d+)")[, 1]
  )
}

triads_numerator <- function(x) {
  m <- str_match(x, "^(\\d+)\\s*/\\s*(\\d+)")
  as.integer(m[, 2])
}

triads_denominator <- function(x) {
  m <- str_match(x, "^(\\d+)\\s*/\\s*(\\d+)")
  as.integer(m[, 3])
}

clean_accession <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    str_detect(x, "Correct ID:\\s*(PRJNA\\d+)") ~ str_match(x, "Correct ID:\\s*(PRJNA\\d+)")[, 2],
    str_detect(x, "PRJEB\\d+") ~ str_match(x, "(PRJEB\\d+)")[, 2],
    str_detect(x, "PRJNA\\d+") ~ str_match(x, "(PRJNA\\d+)")[, 2],
    TRUE ~ NA_character_
  )
}

classify_recommendation <- function(public_reads, metadata_private, triads_included, triad_fraction) {
  case_when(
    !public_reads ~ "avoid_for_now",
    is.na(triads_included) ~ "review_manually",
    triads_included < 5 ~ "avoid_for_now",
    metadata_private & triads_included < 8 ~ "possible_with_extra_curation",
    metadata_private ~ "possible_with_extra_curation",
    triads_included >= 8 & triad_fraction >= 0.7 ~ "recommended_now",
    triads_included >= 5 ~ "possible_with_extra_curation",
    TRUE ~ "avoid_for_now"
  )
}

make_reason <- function(public_reads, metadata_private, triads_included, triad_fraction, source_reads, source_metadata) {
  parts <- character()

  if (public_reads) {
    parts <- c(parts, "public raw reads")
  } else {
    parts <- c(parts, paste0("reads source=", source_reads))
  }

  if (metadata_private) {
    parts <- c(parts, paste0("metadata depends on ", source_metadata))
  } else {
    parts <- c(parts, "metadata derivable from public sources")
  }

  if (!is.na(triads_included) && !is.na(triad_fraction)) {
    parts <- c(parts, sprintf("%d included triads (%.0f%% retained)", triads_included, triad_fraction * 100))
  }

  paste(parts, collapse = "; ")
}

meta_table <- read_excel(input_xlsx, sheet = "Supplementary Table 1") %>%
  filter(!is.na(.data[["Dataset name"]])) %>%
  transmute(
    dataset_name = .data[["Dataset name"]],
    first_author = .data[["First Author"]],
    year = .data[["Year"]],
    location = .data[["Location (of population)"]],
    study_design = .data[["Study design"]],
    recipient_disease = .data[["Recipient disease"]],
    bioproject_accession_raw = .data[["BioProject accession"]],
    bioproject_accession = clean_accession(.data[["BioProject accession"]]),
    source_read_files = .data[["Source of read files"]],
    source_metadata = .data[["Source of metadata"]],
    included_triads_raw = .data[["Included FMT triads (reason for exclusion)"]],
    follow_up = .data[["Follow-up samples after FMT"]],
    clinical_outcome = .data[["Clinical outcome used for analysis"]],
    clinical_data_use = .data[["Use of clinical data for analysis"]]
  ) %>%
  mutate(
    included_triads = triads_numerator(included_triads_raw),
    available_triads = triads_denominator(included_triads_raw),
    triad_fraction = included_triads / available_triads,
    public_reads = str_detect(source_read_files, "NCBI"),
    public_repository = case_when(
      str_detect(bioproject_accession, "^PRJNA") ~ "SRA",
      str_detect(bioproject_accession, "^PRJEB") ~ "ENA",
      TRUE ~ NA_character_
    ),
    metadata_private = str_detect(str_to_lower(source_metadata), "private correspondence"),
    in_house_only = str_detect(str_to_lower(source_read_files), "in-house") |
      str_detect(str_to_lower(source_read_files), "private correspondence"),
    recommendation_tier = classify_recommendation(
      public_reads = public_reads,
      metadata_private = metadata_private,
      triads_included = included_triads,
      triad_fraction = triad_fraction
    ),
    recommendation_reason = mapply(
      FUN = make_reason,
      public_reads = public_reads,
      metadata_private = metadata_private,
      triads_included = included_triads,
      triad_fraction = triad_fraction,
      source_reads = source_read_files,
      source_metadata = source_metadata,
      USE.NAMES = FALSE
    )
  ) %>%
  arrange(
    factor(recommendation_tier, levels = c("recommended_now", "possible_with_extra_curation", "avoid_for_now", "review_manually")),
    desc(included_triads),
    dataset_name
  )

recommended_studies <- meta_table %>%
  filter(recommendation_tier == "recommended_now") %>%
  select(
    dataset_name,
    first_author,
    year,
    recipient_disease,
    bioproject_accession,
    public_repository,
    included_triads,
    available_triads,
    triad_fraction,
    follow_up,
    clinical_data_use,
    recommendation_reason
  )

possible_studies <- meta_table %>%
  filter(recommendation_tier == "possible_with_extra_curation") %>%
  select(
    dataset_name,
    first_author,
    year,
    recipient_disease,
    bioproject_accession,
    public_repository,
    included_triads,
    available_triads,
    triad_fraction,
    source_metadata,
    follow_up,
    recommendation_reason
  )

write_csv(meta_table, file.path(output_dir, "fmt_meta_analysis_study_screen.csv"))
write_csv(recommended_studies, file.path(output_dir, "fmt_meta_analysis_recommended_studies.csv"))
write_csv(possible_studies, file.path(output_dir, "fmt_meta_analysis_possible_studies.csv"))

message("Wrote screening outputs to ", output_dir)
