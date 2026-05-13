#!/usr/bin/env Rscript

required_packages <- c("tidyverse", "Matrix", "SpiecEasi")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running this workflow, e.g. install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(SpiecEasi)
})

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
workflow_root <- file.path(repo_root, "workflows")
workflow_input_dir <- file.path(workflow_root, "input")
workflow_output_dir <- file.path(workflow_root, "output", "graph_inputs")

min_prevalence <- 0.05
min_taxa <- 10
scale_to_counts <- TRUE
scale_factor <- 1e4
max_taxa_for_spieceasi <- Inf
spieceasi_method <- "mb"
nlambda <- 30
lambda_min_ratio <- 1e-2
rep_num <- 20
random_seed <- 1
spieceasi_time_limit_seconds <- Inf
spieceasi_sel_criterion <- "bstars"
spieceasi_ncores <- max(1, parallel::detectCores(logical = FALSE) - 1)

sample_taxon_edges_filename <- "combined_sample_taxon_edges.csv"
sample_nodes_filename <- "nodes_samples.csv"
triplets_filename <- "coalescence_triplets.csv"
taxon_nodes_filename <- "nodes_taxa.csv"
taxon_taxon_edges_filename <- "taxon_taxon_spieceasi_edges.csv"
multirelational_edges_filename <- "graph_edges_multirelational.csv"
summary_filename <- "spieceasi_run_summary.csv"
readme_filename <- "README.md"

log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, paste0(..., collapse = "")))
}

clean_taxonomy_value <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- str_replace(x, "^[a-z]__", "")
  x <- str_replace_all(x, "^\\s+|\\s+$", "")
  x
}

detect_relative_abundance <- function(matrix_df) {
  if (nrow(matrix_df) == 0 || ncol(matrix_df) == 0) {
    return(FALSE)
  }
  row_sums <- rowSums(matrix_df, na.rm = TRUE)
  finite_row_sums <- row_sums[is.finite(row_sums)]
  if (length(finite_row_sums) == 0) {
    return(FALSE)
  }
  max_value <- suppressWarnings(max(as.matrix(matrix_df), na.rm = TRUE))
  all(abs(finite_row_sums - 1) < 1e-6) && is.finite(max_value) && max_value <= 1.5
}

discover_input_files <- function() {
  community_pattern <- "(Fungi|Bacteria)_inoculation_experiment_.*_(donor|resident|final)-community\\.csv$"
  taxonomy_pattern <- "(Fungi|Bacteria)_inoculation_experiment_taxonomy_table\\.csv$"

  preferred_files <- list.files(
    workflow_input_dir,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  fallback_dirs <- c(
    repo_root,
    list.dirs(repo_root, recursive = FALSE, full.names = TRUE)
  )
  fallback_files <- unique(unlist(lapply(
    fallback_dirs,
    function(path) list.files(path, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  )))

  community_files <- unique(c(
    preferred_files[str_detect(basename(preferred_files), community_pattern)],
    fallback_files[str_detect(basename(fallback_files), community_pattern)]
  ))
  taxonomy_files <- unique(c(
    preferred_files[str_detect(basename(preferred_files), taxonomy_pattern)],
    fallback_files[str_detect(basename(fallback_files), taxonomy_pattern)]
  ))

  if (length(community_files) == 0) {
    stop(
      "No microbial coalescence community files were found under workflows/input or the repository root/immediate subdirectories.",
      call. = FALSE
    )
  }

  tibble(
    community_file = sort(community_files),
    community_basename = basename(community_file)
  ) |>
    mutate(parsed = purrr::map(community_basename, parse_community_filename)) |>
    tidyr::unnest_wider(parsed) |>
    mutate(
      taxonomy_file = purrr::map_chr(
        kingdom,
        function(k) {
          match_file <- taxonomy_files[str_detect(basename(taxonomy_files), paste0("^", k, "_inoculation_experiment_taxonomy_table\\.csv$"))]
          if (length(match_file) == 0) "" else normalizePath(match_file[[1]], winslash = "/", mustWork = TRUE)
        }
      )
    )
}

parse_community_filename <- function(filename) {
  pattern <- "^(Fungi|Bacteria)_inoculation_experiment_(.+)_(donor|resident|final)-community\\.csv$"
  matched <- str_match(filename, pattern)
  if (any(is.na(matched))) {
    stop(sprintf("Could not parse community filename metadata from %s", filename), call. = FALSE)
  }
  tibble(
    kingdom = matched[, 2],
    donor_id = matched[, 3],
    community_type = matched[, 4]
  )
}

read_community_matrix <- function(path) {
  warnings <- character()
  log_message("Reading community matrix: ", path)
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2) {
    stop(sprintf("Community matrix %s is empty or missing data rows.", path), call. = FALSE)
  }

  header <- strsplit(lines[1], ",", fixed = TRUE)[[1]]
  header <- gsub('^"|"$', "", header)
  data_rows <- strsplit(lines[-1], ",", fixed = TRUE)
  n_cols <- length(header)
  matrix_rows <- lapply(data_rows, function(parts) {
    length(parts) <- n_cols
    parts
  })
  df <- as_tibble(as.data.frame(do.call(rbind, matrix_rows), stringsAsFactors = FALSE, check.names = FALSE))
  names(df) <- header

  if (ncol(df) < 2) {
    stop(sprintf("Community matrix %s must contain a sample column and at least one taxon column.", path), call. = FALSE)
  }

  original_first_col <- names(df)[1]
  if (is.na(original_first_col) || original_first_col == "" || original_first_col != "sample_id") {
    names(df)[1] <- "sample_id"
    warnings <- c(
      warnings,
      sprintf("Renamed first column from '%s' to 'sample_id' in %s", original_first_col, basename(path))
    )
  }

  df <- df |>
    mutate(sample_id = gsub('^"|"$', "", as.character(sample_id)))

  original_numeric_block <- df[-1]
  numeric_block <- suppressWarnings(original_numeric_block |> mutate(across(everything(), as.numeric)))
  if (any(vapply(names(numeric_block), function(col_name) any(is.na(numeric_block[[col_name]]) & !is.na(original_numeric_block[[col_name]])), logical(1)))) {
    warnings <- c(warnings, sprintf("Some abundance values were coerced to NA in %s", basename(path)))
  }
  df <- bind_cols(df["sample_id"], numeric_block)

  list(
    data = df,
    warnings = warnings,
    is_relative_abundance = detect_relative_abundance(numeric_block)
  )
}

read_taxonomy_table <- function(path) {
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) {
    return(list(data = NULL, warning = "Taxonomy table not found; taxon metadata fields will be left blank."))
  }

  log_message("Reading taxonomy table: ", path)
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) == 0) {
    return(list(data = NULL, warning = sprintf("Taxonomy table %s was empty.", basename(path))))
  }

  names(df)[1] <- "taxon_id"
  missing_cols <- setdiff(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), names(df))
  for (col in missing_cols) {
    df[[col]] <- NA_character_
  }

  df <- df |>
    transmute(
      taxon_id = as.character(taxon_id),
      taxonomy_kingdom = clean_taxonomy_value(Kingdom),
      phylum = clean_taxonomy_value(Phylum),
      class = clean_taxonomy_value(Class),
      order = clean_taxonomy_value(Order),
      family = clean_taxonomy_value(Family),
      genus = clean_taxonomy_value(Genus),
      species = clean_taxonomy_value(Species)
    ) |>
    mutate(
      name = case_when(
        genus != "" & species != "" ~ str_trim(paste(genus, species)),
        genus != "" ~ genus,
        TRUE ~ ""
      )
    ) |>
    distinct(taxon_id, .keep_all = TRUE)

  list(data = df, warning = NULL)
}

build_taxon_nodes <- function(sample_taxon_edges, taxonomy_tables) {
  taxa_from_edges <- sample_taxon_edges |>
    distinct(kingdom, taxon_id)

  taxonomy_combined <- bind_rows(taxonomy_tables) |>
    distinct(kingdom, taxon_id, .keep_all = TRUE)

  taxa_from_edges |>
    left_join(taxonomy_combined, by = c("kingdom", "taxon_id")) |>
    mutate(
      name = replace_na(name, ""),
      genus = replace_na(genus, ""),
      species = replace_na(species, ""),
      phylum = replace_na(phylum, ""),
      class = replace_na(class, ""),
      order = replace_na(order, ""),
      family = replace_na(family, "")
    ) |>
    arrange(kingdom, taxon_id)
}

build_bipartite_edges <- function(metadata_row, community_df) {
  community_df |>
    pivot_longer(-sample_id, names_to = "taxon_id", values_to = "abundance") |>
    filter(!is.na(abundance), abundance != 0) |>
    mutate(
      kingdom = metadata_row$kingdom,
      donor_id = metadata_row$donor_id,
      community_type = metadata_row$community_type
    ) |>
    select(sample_id, taxon_id, abundance, kingdom, donor_id, community_type)
}

prep_spieceasi_matrix <- function(community_df, metadata_row) {
  abundance_df <- community_df |> select(-sample_id)
  n_samples <- nrow(abundance_df)
  prevalence <- colMeans(abundance_df > 0, na.rm = TRUE)
  prevalence <- prevalence[!is.na(prevalence)]

  filtered_taxa <- names(prevalence)[prevalence >= min_prevalence]
  filtered_df <- abundance_df |>
    select(any_of(filtered_taxa))

  prep <- list(
    kingdom = metadata_row$kingdom,
    donor_id = metadata_row$donor_id,
    community_type = metadata_row$community_type,
    n_samples = n_samples,
    n_taxa_original = ncol(abundance_df),
    n_taxa_after_filter = ncol(filtered_df),
    skipped = FALSE,
    message = "",
    matrix = NULL,
    prevalence = prevalence
  )

  if (n_samples < 3) {
    prep$skipped <- TRUE
    prep$message <- "Skipped: fewer than 3 samples in stratum."
    return(prep)
  }

  if (ncol(filtered_df) < min_taxa) {
    prep$skipped <- TRUE
    prep$message <- sprintf("Skipped: %s taxa after prevalence filtering; requires at least %s.", ncol(filtered_df), min_taxa)
    return(prep)
  }

  if (is.finite(max_taxa_for_spieceasi) && ncol(filtered_df) > max_taxa_for_spieceasi) {
    taxa_ranking <- tibble(
      taxon_id = colnames(filtered_df),
      prevalence = unname(prevalence[colnames(filtered_df)]),
      mean_abundance = colMeans(filtered_df, na.rm = TRUE)
    ) |>
      arrange(desc(prevalence), desc(mean_abundance), taxon_id)
    keep_taxa <- taxa_ranking |>
      slice_head(n = max_taxa_for_spieceasi) |>
      pull(taxon_id)
    filtered_df <- filtered_df |>
      select(all_of(keep_taxa))
    prep$message <- paste(
      c(
        prep$message,
        sprintf(
          "Capped taxa from %s to %s for SpiecEasi by ranking taxa on prevalence then mean abundance.",
          nrow(taxa_ranking),
          max_taxa_for_spieceasi
        )
      ),
      collapse = " "
    )
  }

  prep$n_taxa_after_filter <- ncol(filtered_df)

  if (detect_relative_abundance(filtered_df) && isTRUE(scale_to_counts)) {
    filtered_df <- round(filtered_df * scale_factor)
    prep$message <- str_trim(paste(
      prep$message,
      sprintf("Relative abundances were scaled to pseudo-counts with factor %s.", format(scale_factor, scientific = FALSE))
    ))
  }

  # SpiecEasi needs sample-to-sample variation. Pooled donor strata in this
  # workshop dataset often contain repeated identical profiles, so they are
  # intentionally skipped rather than recorded as model failures.
  variable_taxa <- vapply(filtered_df, function(x) {
    x <- as.numeric(x)
    length(unique(x[is.finite(x)])) > 1
  }, logical(1))

  if (sum(variable_taxa) < min_taxa) {
    prep$skipped <- TRUE
    prep$n_taxa_after_filter <- sum(variable_taxa)
    prep$message <- sprintf(
      "Skipped: %s variable taxa after prevalence filtering; requires at least %s. This usually indicates repeated identical sample profiles.",
      sum(variable_taxa),
      min_taxa
    )
    return(prep)
  }

  if (sum(variable_taxa) < ncol(filtered_df)) {
    filtered_df <- filtered_df |> select(all_of(names(variable_taxa)[variable_taxa]))
    prep$message <- str_trim(paste(
      prep$message,
      sprintf("Removed %s non-variable taxa before SpiecEasi.", sum(!variable_taxa))
    ))
  }

  if (any(rowSums(filtered_df, na.rm = TRUE) == 0)) {
    filtered_df <- filtered_df + 1
    prep$message <- str_trim(paste(prep$message, "Added pseudocount of 1 to avoid zero-sum samples."))
  }

  prep$n_taxa_after_filter <- ncol(filtered_df)
  prep$matrix <- as.matrix(filtered_df)
  prep$prevalence <- prevalence[colnames(filtered_df)]
  prep
}

safe_get_spieceasi_matrix <- function(fit, getter_name) {
  getter <- get0(getter_name, envir = asNamespace("SpiecEasi"), mode = "function")
  if (is.null(getter)) {
    return(NULL)
  }
  tryCatch(
    Matrix::as.matrix(getter(fit)),
    error = function(err) NULL
  )
}

build_spieceasi_edges <- function(fit, taxa_names, prepped, input_file) {
  adjacency <- safe_get_spieceasi_matrix(fit, "getRefit")
  if (is.null(adjacency) || is.null(dim(adjacency)) || any(dim(adjacency) == 0)) {
    return(tibble())
  }

  adjacency[is.na(adjacency)] <- 0
  diag(adjacency) <- 0
  edge_positions <- which(upper.tri(adjacency) & adjacency != 0, arr.ind = TRUE)
  if (nrow(edge_positions) == 0) {
    return(tibble())
  }

  weight_matrix <- abs(adjacency)
  sign_matrix <- matrix(NA_real_, nrow = nrow(adjacency), ncol = ncol(adjacency))
  method_label <- paste0("spieceasi_", spieceasi_method)

  if (spieceasi_method == "mb") {
    beta <- safe_get_spieceasi_matrix(fit, "getOptBeta")
    if (!is.null(beta) && all(dim(beta) == dim(adjacency))) {
      beta[is.na(beta)] <- 0
      weight_matrix <- pmax(abs(beta), abs(t(beta)))
      signed_matrix <- beta + t(beta)
      sign_matrix <- sign(signed_matrix)
      sign_matrix[sign_matrix == 0] <- NA_real_
    }
  } else if (spieceasi_method == "glasso") {
    theta <- safe_get_spieceasi_matrix(fit, "getOptTheta")
    if (is.null(theta)) {
      theta <- safe_get_spieceasi_matrix(fit, "getOptiCov")
    }
    if (is.null(theta)) {
      theta <- safe_get_spieceasi_matrix(fit, "getOptCov")
    }
    if (!is.null(theta) && all(dim(theta) == dim(adjacency))) {
      diag_theta <- diag(theta)
      if (all(is.finite(diag_theta)) && all(diag_theta > 0)) {
        inv_sqrt_diag <- diag(1 / sqrt(diag_theta))
        partial_corr <- -inv_sqrt_diag %*% theta %*% inv_sqrt_diag
        diag(partial_corr) <- 0
        weight_matrix <- abs(partial_corr)
        sign_matrix <- sign(partial_corr)
      }
    }
  }

  tibble(
    kingdom = prepped$kingdom,
    donor_id = prepped$donor_id,
    community_type = prepped$community_type,
    taxon_a = taxa_names[edge_positions[, "row"]],
    taxon_b = taxa_names[edge_positions[, "col"]],
    weight = weight_matrix[edge_positions],
    sign = sign_matrix[edge_positions],
    method = method_label,
    n_samples = prepped$n_samples,
    prevalence_a = unname(prepped$prevalence[taxa_names[edge_positions[, "row"]]]),
    prevalence_b = unname(prepped$prevalence[taxa_names[edge_positions[, "col"]]])
  ) |>
    arrange(kingdom, donor_id, community_type, taxon_a, taxon_b)
}

run_spieceasi_for_stratum <- function(prepped, input_file) {
  summary_row <- tibble(
    kingdom = prepped$kingdom,
    donor_id = prepped$donor_id,
    community_type = prepped$community_type,
    input_file = input_file,
    n_samples = prepped$n_samples,
    n_taxa_original = prepped$n_taxa_original,
    n_taxa_after_filter = prepped$n_taxa_after_filter,
    n_edges = 0L,
    status = "skipped",
    message = prepped$message
  )

  if (isTRUE(prepped$skipped)) {
    return(list(edges = tibble(), summary = summary_row))
  }

  if (!spieceasi_method %in% c("glasso", "mb")) {
    summary_row$status <- "failed"
    summary_row$message <- sprintf("Unsupported spieceasi_method '%s'. Supported methods are 'glasso' and 'mb'.", spieceasi_method)
    return(list(edges = tibble(), summary = summary_row))
  }

  log_message(
    "Running SpiecEasi for ",
    prepped$kingdom, " / ", prepped$donor_id, " / ", prepped$community_type,
    " with ", prepped$n_samples, " samples and ", prepped$n_taxa_after_filter, " taxa."
  )

  result <- tryCatch(
    {
      set.seed(random_seed)
      if (is.finite(spieceasi_time_limit_seconds)) {
        setTimeLimit(elapsed = spieceasi_time_limit_seconds, transient = TRUE)
        on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
      }

      captured_warnings <- character()
      fit <- withCallingHandlers(
        spiec.easi(
          prepped$matrix,
          method = spieceasi_method,
          nlambda = nlambda,
          lambda.min.ratio = lambda_min_ratio,
          sel.criterion = spieceasi_sel_criterion,
          pulsar.params = list(rep.num = rep_num, ncores = spieceasi_ncores),
          verbose = FALSE
        ),
        warning = function(w) {
          captured_warnings <<- c(captured_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      edges <- build_spieceasi_edges(fit, colnames(prepped$matrix), prepped, input_file)
      summary_row$n_edges <- nrow(edges)
      summary_row$status <- "success"
      warning_text <- if (length(captured_warnings) > 0) {
        paste("Warnings:", paste(unique(captured_warnings), collapse = " | "))
      } else {
        ""
      }
      summary_row$message <- str_trim(paste("SpiecEasi completed successfully.", prepped$message, warning_text))
      list(edges = edges, summary = summary_row)
    },
    error = function(err) {
      summary_row$status <- "failed"
      summary_row$message <- paste("SpiecEasi failed:", conditionMessage(err))
      list(edges = tibble(), summary = summary_row)
    }
  )

  result
}

build_multirelational_edges <- function(sample_taxon_edges, taxon_taxon_edges) {
  sample_edges <- sample_taxon_edges |>
    transmute(
      kingdom,
      donor_id,
      community_type,
      from = sample_id,
      to = taxon_id,
      weight = abundance,
      sign = NA_real_,
      relation = "sample_taxon",
      method = "abundance"
    )

  if (nrow(taxon_taxon_edges) == 0) {
    return(sample_edges)
  }

  taxon_edges <- taxon_taxon_edges |>
    transmute(
      kingdom,
      donor_id,
      community_type,
      from = taxon_a,
      to = taxon_b,
      weight,
      sign,
      relation = "taxon_taxon",
      method
    )

  bind_rows(sample_edges, taxon_edges)
}

build_coalescence_triplets <- function(sample_nodes) {
  required_cols <- c("sample_id", "kingdom", "donor_id", "community_type")
  missing_cols <- setdiff(required_cols, names(sample_nodes))
  if (length(missing_cols) > 0) {
    warning("Cannot build coalescence triplets because nodes_samples is missing columns: ", paste(missing_cols, collapse = ", "))
    return(tibble())
  }

  resident_ids <- sample_nodes |>
    filter(community_type == "resident") |>
    select(kingdom, donor_id, sample_id) |>
    rename(resident_sample_id = sample_id)

  final_ids <- sample_nodes |>
    filter(community_type == "final") |>
    select(kingdom, donor_id, sample_id) |>
    rename(final_sample_id = sample_id)

  donor_counts <- sample_nodes |>
    filter(community_type == "donor") |>
    count(kingdom, donor_id, name = "n_donor_samples")

  if (nrow(resident_ids) == 0 || nrow(final_ids) == 0) {
    warning("Cannot build coalescence triplets because resident or final samples were not found.")
    return(tibble())
  }

  triplets <- inner_join(
    resident_ids,
    final_ids,
    by = c("kingdom", "donor_id"),
    relationship = "many-to-many"
  ) |>
    filter(resident_sample_id == final_sample_id) |>
    mutate(
      microcosm_id = paste(kingdom, donor_id, resident_sample_id, sep = "::"),
      donor_source_id = donor_id,
      donor_sample_id = NA_character_,
      donor_is_pooled = TRUE,
      pairing_basis = "resident_final_same_sample_id; donor_source_by_donor_id"
    ) |>
    left_join(donor_counts, by = c("kingdom", "donor_id")) |>
    mutate(n_donor_samples = replace_na(n_donor_samples, 0L)) |>
    select(
      microcosm_id,
      kingdom,
      donor_id,
      donor_source_id,
      donor_sample_id,
      resident_sample_id,
      final_sample_id,
      donor_is_pooled,
      n_donor_samples,
      pairing_basis
    ) |>
    arrange(kingdom, donor_id, resident_sample_id)

  if (nrow(triplets) == 0) {
    warning("No resident-final sample pairs matched by identical sample_id within kingdom and donor_id.")
  }

  triplets
}

write_readme <- function(output_dir, discovered_files, sample_taxon_edges, taxon_nodes, spieceasi_edges, summary_tbl, triplets, workflow_notes) {
  detected_community_files <- paste0("- `", basename(discovered_files$community_file), "`")
  detected_taxonomy_files <- discovered_files |>
    distinct(kingdom, taxonomy_file) |>
    mutate(
      line = if_else(
        taxonomy_file == "",
        paste0("- `", kingdom, "`: taxonomy table not found"),
        paste0("- `", kingdom, "`: `", basename(taxonomy_file), "`")
      )
    ) |>
    pull(line)

  successful_runs <- summary_tbl |> filter(status == "success")
  skipped_runs <- summary_tbl |> filter(status == "skipped")
  failed_runs <- summary_tbl |> filter(status == "failed")

  readme_lines <- c(
    "# Graph Inputs Workflow",
    "",
    "## What This Workflow Does",
    "",
    "This workflow discovers microbial coalescence community matrices and taxonomy tables, validates their structure, builds a bipartite sample-taxon backbone, fits stratum-specific SpiecEasi taxon-taxon networks, and writes graph-ready CSV files for downstream graph neural network workflows.",
    "",
    "## Inputs Detected",
    "",
    "Community matrices were detected from `workflows/input/` first and then the repository root/immediate subdirectories as fallback.",
    detected_community_files,
    "",
    "Taxonomy tables used:",
    detected_taxonomy_files,
    "",
    "## Filename Metadata Parsing",
    "",
    "Community filenames were parsed with the pattern `Kingdom_inoculation_experiment_<donor_id>_<community_type>-community.csv`.",
    "This workflow uses:",
    "- `kingdom`: `Bacteria` or `Fungi`",
    "- `donor_id`: the middle filename token such as `0Burn_W1`",
    "- `community_type`: `donor`, `resident`, or `final`",
    "",
    "## Sample-Taxon Bipartite Backbone",
    "",
    "Each community matrix is treated as samples in rows and taxa/features in columns.",
    "The first column is used as `sample_id`; if it was not already named `sample_id`, the workflow renames it and records that warning in the summary outputs.",
    "Only nonzero abundance entries are written to the bipartite edge table.",
    "",
    "## Taxon Naming",
    "",
    "Taxon node metadata comes from the kingdom-specific taxonomy tables when available.",
    "Scientific names are constructed from `Genus` and `Species` only.",
    "- `name = \"Genus Species\"` when both are present",
    "- `name = \"Genus\"` when only genus is present",
    "- `name = \"\"` when genus is missing",
    "Higher taxonomy is retained only as optional metadata fields, not as the primary name.",
    "",
    "## Coalescence Triplets",
    "",
    "The workflow also writes `coalescence_triplets.csv`, which defines the supervised prediction units for this workshop dataset.",
    "Resident and final communities are paired when they share the same `sample_id` within a given `kingdom` and `donor_id`.",
    "Donor communities are represented at the pooled donor-source level: all resident/final samples with the same `donor_id` share the same donor-source input.",
    "Accordingly, `donor_sample_id` is `NA`, `donor_source_id` stores the donor treatment/source, and `donor_is_pooled` is `TRUE`.",
    "",
    "## SpiecEasi Networks",
    "",
    "SpiecEasi is run separately for each `kingdom x donor_id x community_type` stratum so that bacteria and fungi are not pooled, donor/resident/final communities are not pooled, and donor IDs are not pooled.",
    paste0("The configured SpiecEasi method is `", spieceasi_method, "` with selection criterion `", spieceasi_sel_criterion, "`."),
    "When a stratum retains more taxa than the workshop cap allows, taxa are ranked by prevalence and then mean abundance before taking the top subset for SpiecEasi. This cap affects only the inferred taxon-taxon layer; the full bipartite backbone is preserved.",
    "For `method = mb`, SpiecEasi uses neighborhood selection. The resulting taxon-taxon edges should be interpreted as inferred association structure, not direct ecological interactions or strict partial correlations.",
    "For `method = glasso`, the selected inverse covariance structure can be converted to partial-correlation-like associations. This script stores edge weights as absolute association strength and `sign` as the association direction from the selected SpiecEasi matrix.",
    "",
    "Workflow parameters:",
    paste0("- `min_prevalence = ", min_prevalence, "`"),
    paste0("- `min_taxa = ", min_taxa, "`"),
    paste0("- `scale_to_counts = ", scale_to_counts, "`"),
    paste0("- `scale_factor = ", format(scale_factor, scientific = FALSE), "`"),
    paste0("- `max_taxa_for_spieceasi = ", max_taxa_for_spieceasi, "`"),
    paste0("- `spieceasi_time_limit_seconds = ", spieceasi_time_limit_seconds, "`"),
    paste0("- `spieceasi_method = \"", spieceasi_method, "\"`"),
    paste0("- `spieceasi_sel_criterion = \"", spieceasi_sel_criterion, "\"`"),
    paste0("- `spieceasi_ncores = ", spieceasi_ncores, "`"),
    paste0("- `nlambda = ", nlambda, "`"),
    paste0("- `lambda_min_ratio = ", lambda_min_ratio, "`"),
    paste0("- `rep_num = ", rep_num, "`"),
    paste0("- `random_seed = ", random_seed, "`"),
    "",
    "## Output Files",
    "",
    paste0("- `", sample_taxon_edges_filename, "`: nonzero sample-taxon bipartite edges with abundance and experimental context."),
    paste0("- `", sample_nodes_filename, "`: unique sample nodes with kingdom, donor ID, and community type."),
    paste0("- `", triplets_filename, "`: supervised coalescence units linking resident and final samples by shared sample ID within kingdom and donor source; donor input is represented by pooled donor source ID."),
    paste0("- `", taxon_nodes_filename, "`: unique taxon nodes with names and taxonomy metadata."),
    paste0("- `", taxon_taxon_edges_filename, "`: undirected SpiecEasi taxon-taxon edges with inferred association weights and signs."),
    paste0("- `", multirelational_edges_filename, "`: flat edge file combining sample-taxon and taxon-taxon relations."),
    paste0("- `", summary_filename, "`: one row per stratum describing network status and skip/failure messages."),
    "",
    "## Rerun",
    "",
    "Run the workflow from the repository root with:",
    "",
    "```bash",
    "Rscript workflows/code/build_graph_inputs.R",
    "```",
    "",
    "## GNN Use",
    "",
    "The bipartite backbone captures observed experimental composition data, while the SpiecEasi layer adds inferred ecological association structure. The combined multirelational edge file can be used as a starting point for heterogeneous or relational GNN pipelines that link samples, taxa, and inferred taxon-taxon associations.",
    "",
    "## Caveats",
    "",
    "- Co-occurrence edges are inferred statistical associations, not measured direct interactions.",
    "- Positive co-occurrence does not necessarily mean cooperation.",
    "- Negative co-occurrence does not necessarily mean inhibition.",
    "- Associations can reflect shared niches, environmental filtering, compositional effects, or indirect interactions.",
    "- Prevalence filtering affects network density.",
    "- Relative-abundance to pseudo-count conversion is a modeling choice.",
    "- The bipartite backbone is the experimental data representation; the SpiecEasi graph is an inferred ecological association layer.",
    "- With `method = mb`, edge weights are association-strength summaries from neighborhood selection and should not be described as strict partial correlations.",
    "- Donor-community SpiecEasi strata may be intentionally skipped if repeated donor profiles have no sample-to-sample variation; donor composition is still represented in the sample-taxon backbone.",
    "",
    "## Run Summary",
    "",
    paste0("- Sample-taxon edges: `", format(nrow(sample_taxon_edges), big.mark = ","), "`"),
    paste0("- Sample nodes: `", format(n_distinct(sample_taxon_edges$sample_id), big.mark = ","), "`"),
    paste0("- Coalescence triplets: `", format(nrow(triplets), big.mark = ","), "`"),
    paste0("- Taxon nodes: `", format(nrow(taxon_nodes), big.mark = ","), "`"),
    paste0("- SpiecEasi edges: `", format(nrow(spieceasi_edges), big.mark = ","), "`"),
    paste0("- Successful SpiecEasi strata: `", nrow(successful_runs), "`"),
    paste0("- Skipped SpiecEasi strata: `", nrow(skipped_runs), "`"),
    paste0("- Failed SpiecEasi strata: `", nrow(failed_runs), "`")
  )

  if (length(workflow_notes) > 0) {
    readme_lines <- c(
      readme_lines,
      "",
      "## Notes",
      "",
      paste0("- ", workflow_notes)
    )
  }

  readr::write_lines(readme_lines, file.path(output_dir, readme_filename))
}

main <- function() {
  dir.create(workflow_output_dir, recursive = TRUE, showWarnings = FALSE)
  log_message("Discovering workflow inputs.")
  discovered_files <- discover_input_files()

  taxonomy_cache <- list()
  taxonomy_warnings <- character()
  for (taxonomy_path in unique(discovered_files$taxonomy_file)) {
    key <- if (taxonomy_path == "") "missing" else taxonomy_path
    if (!key %in% names(taxonomy_cache)) {
      taxonomy_result <- read_taxonomy_table(taxonomy_path)
      taxonomy_cache[[key]] <- taxonomy_result
      if (!is.null(taxonomy_result$warning)) {
        taxonomy_warnings <- c(taxonomy_warnings, taxonomy_result$warning)
      }
    }
  }

  community_results <- pmap(
    discovered_files,
    function(community_file, community_basename, kingdom, donor_id, community_type, taxonomy_file) {
      metadata_row <- tibble(
        community_file = community_file,
        community_basename = community_basename,
        kingdom = kingdom,
        donor_id = donor_id,
        community_type = community_type,
        taxonomy_file = taxonomy_file
      )
      community_result <- read_community_matrix(community_file)
      community_result$metadata <- metadata_row
      community_result
    }
  )

  sample_taxon_edges <- map2_dfr(
    community_results,
    seq_along(community_results),
    function(result, idx) build_bipartite_edges(result$metadata, result$data)
  )

  sample_nodes <- sample_taxon_edges |>
    distinct(sample_id, kingdom, donor_id, community_type) |>
    arrange(kingdom, donor_id, community_type, sample_id)

  coalescence_triplets <- build_coalescence_triplets(sample_nodes)

  taxonomy_tables <- map2_dfr(
    community_results,
    seq_along(community_results),
    function(result, idx) {
      taxonomy_key <- if (result$metadata$taxonomy_file == "") "missing" else result$metadata$taxonomy_file
      taxonomy_df <- taxonomy_cache[[taxonomy_key]]$data
      if (is.null(taxonomy_df)) {
        return(tibble())
      }
      taxonomy_df |>
        mutate(kingdom = result$metadata$kingdom) |>
        select(kingdom, everything())
    }
  ) |>
    distinct(kingdom, taxon_id, .keep_all = TRUE)

  taxon_nodes <- build_taxon_nodes(sample_taxon_edges, taxonomy_tables)

  readr::write_csv(sample_taxon_edges, file.path(workflow_output_dir, sample_taxon_edges_filename))
  readr::write_csv(sample_nodes, file.path(workflow_output_dir, sample_nodes_filename))
  readr::write_csv(coalescence_triplets, file.path(workflow_output_dir, triplets_filename))
  readr::write_csv(taxon_nodes, file.path(workflow_output_dir, taxon_nodes_filename))

  spieceasi_results <- map(
    community_results,
    function(result) {
      prepped <- prep_spieceasi_matrix(result$data, result$metadata)
      run_spieceasi_for_stratum(prepped, result$metadata$community_file)
    }
  )

  spieceasi_edges <- bind_rows(map(spieceasi_results, "edges"))
  summary_tbl <- bind_rows(map(spieceasi_results, "summary")) |>
    arrange(kingdom, donor_id, community_type)

  multirelational_edges <- build_multirelational_edges(sample_taxon_edges, spieceasi_edges)

  readr::write_csv(spieceasi_edges, file.path(workflow_output_dir, taxon_taxon_edges_filename))
  readr::write_csv(multirelational_edges, file.path(workflow_output_dir, multirelational_edges_filename))
  readr::write_csv(summary_tbl, file.path(workflow_output_dir, summary_filename))

  workflow_notes <- c(
    unique(unlist(map(community_results, "warnings"))),
    taxonomy_warnings
  )
  if (nrow(spieceasi_edges) == 0) {
    workflow_notes <- c(workflow_notes, "All SpiecEasi strata were skipped or failed, so the multirelational edge file contains sample-taxon edges only.")
  }
  if (nrow(coalescence_triplets) == 0) {
    workflow_notes <- c(workflow_notes, "No resident-final coalescence triplets were identified by matching sample_id within kingdom and donor_id.")
  }

  write_readme(
    output_dir = workflow_output_dir,
    discovered_files = discovered_files,
    sample_taxon_edges = sample_taxon_edges,
    taxon_nodes = taxon_nodes,
    spieceasi_edges = spieceasi_edges,
    summary_tbl = summary_tbl,
    triplets = coalescence_triplets,
    workflow_notes = unique(workflow_notes[nzchar(workflow_notes)])
  )

  log_message("Graph input workflow completed. Outputs written to ", workflow_output_dir)
}

main()
