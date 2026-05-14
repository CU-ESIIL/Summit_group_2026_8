#!/usr/bin/env Rscript

# Visualize saved GNN performance outputs.
#
# This script is intentionally separate from model training/compilation. It
# only reads saved CSV outputs from workflows/output/gnn_model and writes plots,
# tidy summary tables, and a short README into a visualization subdirectory.
#
# Run from the repository root:
#   Rscript workflows/code/visualize_gnn_performance_CMKR.R

required_packages <- c("tidyverse")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install before running this workflow, e.g. install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
})

parse_cli_args <- function(args) {
  parsed <- list()
  for (arg in args) {
    arg <- sub("^--", "", arg)
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
    if (length(pieces) == 1) {
      parsed[[pieces[[1]]]] <- TRUE
    } else {
      key <- pieces[[1]]
      value <- paste(pieces[-1], collapse = "=")
      parsed[[key]] <- value
    }
  }
  parsed
}

arg_value <- function(args, name, default) {
  value <- args[[name]]
  if (is.null(value) || length(value) == 0 || is.na(value) || identical(value, "")) {
    default
  } else {
    value
  }
}

log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, paste0(..., collapse = "")))
}

setwd("C:/Users/coope/OneDrive/Desktop/esiil/Summit_group_2026_8")

repo_root <- normalizePath("../..", winslash = "/", mustWork = TRUE)
cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

model_dir <- normalizePath(
  file.path(repo_root, "workflows", "output", "gnn_model"),
  winslash = "/",
  mustWork = TRUE
)

output_dir <- file.path(
  model_dir,
  arg_value(cli_args, "output_subdir", "performance_visualizations_CMKR")
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

log_message("Reading saved GNN metrics from: ", model_dir)
log_message("Writing visualizations to: ", output_dir)

# Explicitly load the model performance files produced by the training script.
metrics_path <- file.path(model_dir, "evaluation_metrics_all_kingdoms.csv")
if (!file.exists(metrics_path)) {
  stop(sprintf("Missing required metrics file: %s", metrics_path), call. = FALSE)
}

metrics <- read_csv("workflows/output/gnn_model/evaluation_metrics_all_kingdoms.csv", show_col_types = FALSE) |>
  mutate(
    split = factor(split, levels = c("train", "validation", "test")),
    model = kingdom
  )

history_paths <- list.files(
  model_dir,
  pattern = "^training_history\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
histories <- history_paths |>
  purrr::map(readr::read_csv, show_col_types = FALSE) |>
  purrr::list_rbind()

graph_summary_paths <- list.files(
  model_dir,
  pattern = "^spieceasi_graph_summary\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
graph_summaries <- graph_summary_paths |>
  purrr::map(readr::read_csv, show_col_types = FALSE) |>
  purrr::list_rbind()

required_metric_columns <- c(
  "kingdom",
  "split",
  "n_microcosms",
  "presence_accuracy",
  "presence_precision",
  "presence_recall",
  "presence_specificity",
  "presence_f1",
  "rmse_log_abundance_present_taxa",
  "mae_log_abundance_present_taxa",
  "mean_bray_curtis",
  "median_bray_curtis"
)
missing_metric_columns <- setdiff(required_metric_columns, names(metrics))
if (length(missing_metric_columns) > 0) {
  stop(
    paste0(
      "Evaluation metrics file is missing required column(s): ",
      paste(missing_metric_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

classification_labels <- c(
  presence_accuracy = "Accuracy",
  presence_precision = "Precision",
  presence_recall = "Recall",
  presence_specificity = "Specificity",
  presence_f1 = "F1"
)
error_labels <- c(
  rmse_log_abundance_present_taxa = "RMSE log abundance",
  mae_log_abundance_present_taxa = "MAE log abundance",
  mean_bray_curtis = "Mean Bray-Curtis",
  median_bray_curtis = "Median Bray-Curtis"
)

classification_long <- metrics |>
  select(kingdom, model, split, all_of(names(classification_labels))) |>
  pivot_longer(
    cols = all_of(names(classification_labels)),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(metric_label = recode(metric, !!!classification_labels))

error_long <- metrics |>
  select(kingdom, model, split, all_of(names(error_labels))) |>
  pivot_longer(
    cols = all_of(names(error_labels)),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(metric_label = recode(metric, !!!error_labels))

combined_summary <- metrics |>
  filter(split %in% c("validation", "test")) |>
  group_by(split) |>
  summarise(
    n_models = n(),
    total_microcosms = sum(n_microcosms, na.rm = TRUE),
    weighted_presence_accuracy = weighted.mean(presence_accuracy, n_microcosms, na.rm = TRUE),
    weighted_presence_f1 = weighted.mean(presence_f1, n_microcosms, na.rm = TRUE),
    weighted_mean_bray_curtis = weighted.mean(mean_bray_curtis, n_microcosms, na.rm = TRUE),
    weighted_rmse_log_abundance_present_taxa = weighted.mean(
      rmse_log_abundance_present_taxa,
      n_microcosms,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

readr::write_csv(combined_summary, file.path(output_dir, "combined_performance_summary.csv"))
readr::write_csv(classification_long, file.path(output_dir, "classification_metrics_long.csv"))
readr::write_csv(error_long, file.path(output_dir, "error_metrics_long.csv"))

theme_cmkr <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

save_plot <- function(plot, filename, width = 10, height = 6) {
  png_path <- file.path(output_dir, paste0(filename, ".png"))
  pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300)
  ggsave(pdf_path, plot, width = width, height = height)
  log_message("Wrote ", png_path)
}

classification_plot <- classification_long |>
  ggplot(aes(x = split, y = value, fill = kingdom)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  facet_wrap(~metric_label, ncol = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Final-community presence prediction performance",
    subtitle = "Higher values indicate better classification of taxon presence after coalescence",
    x = "Data split",
    y = "Metric value",
    fill = "Model"
  ) +
  theme_cmkr()

save_plot(classification_plot, "presence_classification_metrics_by_split", width = 11, height = 7)

error_plot <- error_long |>
  ggplot(aes(x = split, y = value, fill = kingdom)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  facet_wrap(~metric_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Final-community abundance and composition error",
    subtitle = "Lower values indicate better abundance or community-composition prediction",
    x = "Data split",
    y = "Error value",
    fill = "Model"
  ) +
  theme_cmkr()
save_plot(error_plot, "abundance_and_bray_curtis_errors_by_split", width = 10, height = 7)

test_focus <- metrics |>
  filter(split == "test")
if (nrow(test_focus) > 0) {
  test_tradeoff <- test_focus |>
    ggplot(aes(
      x = presence_f1,
      y = mean_bray_curtis,
      size = n_microcosms,
      label = kingdom,
      color = kingdom
    )) +
    geom_point(alpha = 0.85) +
    geom_text(vjust = -0.8, show.legend = FALSE) +
    scale_size_continuous(range = c(4, 10)) +
    labs(
      title = "Held-out model performance tradeoff",
      subtitle = "Better models move toward higher F1 and lower Bray-Curtis distance",
      x = "Test presence F1",
      y = "Test mean Bray-Curtis distance",
      color = "Model",
      size = "Test microcosms"
    ) +
    theme_cmkr()
  save_plot(test_tradeoff, "test_f1_vs_bray_curtis_tradeoff", width = 8, height = 6)
}

if (nrow(histories) > 0) {
  loss_plot <- histories |>
    ggplot(aes(x = epoch, y = train_loss, color = kingdom)) +
    geom_line(linewidth = 0.8) +
    labs(
      title = "Training loss by epoch",
      subtitle = "Useful for diagnosing convergence, instability, or under-training",
      x = "Epoch",
      y = "Training loss",
      color = "Model"
    ) +
    theme_cmkr() +
    theme(axis.text.x = element_text(angle = 0))
  save_plot(loss_plot, "training_loss_by_epoch", width = 9, height = 5.5)

  validation_long <- histories |>
    select(
      kingdom,
      epoch,
      validation_presence_f1,
      validation_mean_bray_curtis,
      validation_rmse_log_abundance_present_taxa
    ) |>
    pivot_longer(
      cols = starts_with("validation_"),
      names_to = "metric",
      values_to = "value"
    ) |>
    mutate(
      metric_label = recode(
        metric,
        validation_presence_f1 = "Validation presence F1",
        validation_mean_bray_curtis = "Validation mean Bray-Curtis",
        validation_rmse_log_abundance_present_taxa = "Validation RMSE log abundance"
      )
    )

  validation_plot <- validation_long |>
    ggplot(aes(x = epoch, y = value, color = kingdom)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~metric_label, scales = "free_y", ncol = 1) +
    labs(
      title = "Validation performance during training",
      subtitle = "F1 should increase; Bray-Curtis and RMSE should decrease when training is improving",
      x = "Epoch",
      y = "Metric value",
      color = "Model"
    ) +
    theme_cmkr() +
    theme(axis.text.x = element_text(angle = 0))
  save_plot(validation_plot, "validation_metrics_by_epoch", width = 9, height = 8)
}

if (nrow(graph_summaries) > 0) {
  graph_labels <- c(
    n_edges_after_taxon_filter = "Edges retained",
    n_positive_edges = "Positive edges",
    n_negative_edges = "Negative edges",
    n_taxa_with_positive_neighbors = "Taxa with positive neighbors",
    n_taxa_with_negative_neighbors = "Taxa with negative neighbors"
  )

  graph_long <- graph_summaries |>
    pivot_longer(
      cols = all_of(names(graph_labels)),
      names_to = "metric",
      values_to = "value"
    ) |>
    mutate(metric_label = recode(metric, !!!graph_labels))

  readr::write_csv(graph_long, file.path(output_dir, "spieceasi_graph_summary_long.csv"))

  graph_plot <- graph_long |>
    ggplot(aes(x = kingdom, y = value, fill = kingdom)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    facet_wrap(~metric_label, scales = "free_y", ncol = 3) +
    labs(
      title = "SpiecEasi graph scaffold used by each model",
      subtitle = "Counts are after kingdom, taxon, and community-type filtering",
      x = "Model",
      y = "Count"
    ) +
    theme_cmkr()
  save_plot(graph_plot, "spieceasi_graph_scaffold_summary", width = 11, height = 6)
}

readme_lines <- c(
  "# GNN performance visualizations",
  "",
  sprintf("Source model directory: `%s`", model_dir),
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Files",
  "",
  "- `presence_classification_metrics_by_split.png/.pdf`: final-community presence classification metrics by kingdom and split.",
  "- `abundance_and_bray_curtis_errors_by_split.png/.pdf`: abundance and community-composition error metrics by kingdom and split.",
  "- `test_f1_vs_bray_curtis_tradeoff.png/.pdf`: held-out F1 versus Bray-Curtis distance.",
  "- `training_loss_by_epoch.png/.pdf`: training loss curves, when training histories are available.",
  "- `validation_metrics_by_epoch.png/.pdf`: validation performance curves, when training histories are available.",
  "- `spieceasi_graph_scaffold_summary.png/.pdf`: graph scaffold counts, when SpiecEasi graph summaries are available.",
  "- `combined_performance_summary.csv`: microcosm-weighted validation/test performance across kingdom-specific models.",
  "",
  "## Interpretation",
  "",
  "- Higher presence accuracy, precision, recall, specificity, and F1 indicate better final-presence prediction.",
  "- Lower log-abundance RMSE/MAE indicates better abundance prediction for taxa observed in the final community.",
  "- Lower Bray-Curtis distance indicates predicted final communities are compositionally closer to observed final communities.",
  "- Test metrics are the clearest held-out performance summary; training curves are diagnostic."
)
writeLines(readme_lines, file.path(output_dir, "README.md"))

log_message("Done. Visualization outputs written to: ", output_dir)
