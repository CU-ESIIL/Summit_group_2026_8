#!/usr/bin/env Rscript

# Visualize performance outputs from train_coalescence_gnn_CMKR.R.
#
# The script reads model metrics from workflows/output/gnn_model, then writes
# publication/workshop-friendly plots and compact summary tables back into the
# model output directory.
#
# Run from the repository root:
#   Rscript workflows/code/visualize_gnn_performance_CMKR.R
#
required_packages <- c("data.table", "ggplot2")
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
  library(data.table)
  library(ggplot2)
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

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
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

log_message("Reading GNN model outputs from: ", model_dir)
log_message("Writing visualizations to: ", output_dir)

metrics_path <- file.path(model_dir, "evaluation_metrics_all_kingdoms.csv")
if (!file.exists(metrics_path)) {
  metric_candidates <- list.files(
    model_dir,
    pattern = "^evaluation_metrics\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(metric_candidates) == 0) {
    stop(
      sprintf("No evaluation metrics found in %s", model_dir),
      call. = FALSE
    )
  }
  metrics <- rbindlist(lapply(metric_candidates, fread), fill = TRUE)
} else {
  metrics <- fread(metrics_path)
}

history_paths <- list.files(
  model_dir,
  pattern = "^training_history\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
histories <- if (length(history_paths) > 0) {
  rbindlist(lapply(history_paths, fread), fill = TRUE)
} else {
  data.table()
}

graph_summary_paths <- list.files(
  model_dir,
  pattern = "^spieceasi_graph_summary\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
graph_summaries <- if (length(graph_summary_paths) > 0) {
  rbindlist(lapply(graph_summary_paths, fread), fill = TRUE)
} else {
  data.table()
}

required_metric_columns <- c(
  "kingdom",
  "split",
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

metrics[, split := factor(split, levels = c("train", "validation", "test"))]
metrics[, model := kingdom]

classification_metrics <- c(
  "presence_accuracy",
  "presence_precision",
  "presence_recall",
  "presence_specificity",
  "presence_f1"
)
error_metrics <- c(
  "rmse_log_abundance_present_taxa",
  "mae_log_abundance_present_taxa",
  "mean_bray_curtis",
  "median_bray_curtis"
)

classification_long <- melt(
  metrics,
  id.vars = c("kingdom", "model", "split"),
  measure.vars = classification_metrics,
  variable.name = "metric",
  value.name = "value"
)
classification_long[, metric_label := fifelse(
  metric == "presence_accuracy",
  "Accuracy",
  fifelse(
    metric == "presence_precision",
    "Precision",
    fifelse(
      metric == "presence_recall",
      "Recall",
      fifelse(
        metric == "presence_specificity",
        "Specificity",
        "F1"
      )
    )
  )
)]

error_long <- melt(
  metrics,
  id.vars = c("kingdom", "model", "split"),
  measure.vars = error_metrics,
  variable.name = "metric",
  value.name = "value"
)
error_long[, metric_label := fifelse(
  metric == "rmse_log_abundance_present_taxa",
  "RMSE log abundance",
  fifelse(
    metric == "mae_log_abundance_present_taxa",
    "MAE log abundance",
    fifelse(metric == "mean_bray_curtis", "Mean Bray-Curtis", "Median Bray-Curtis")
  )
)]

combined_summary <- metrics[
  split %in% c("validation", "test"),
  .(
    n_models = .N,
    total_microcosms = sum(n_microcosms, na.rm = TRUE),
    weighted_presence_accuracy = stats::weighted.mean(presence_accuracy, n_microcosms, na.rm = TRUE),
    weighted_presence_f1 = stats::weighted.mean(presence_f1, n_microcosms, na.rm = TRUE),
    weighted_mean_bray_curtis = stats::weighted.mean(mean_bray_curtis, n_microcosms, na.rm = TRUE),
    weighted_rmse_log_abundance_present_taxa = stats::weighted.mean(
      rmse_log_abundance_present_taxa,
      n_microcosms,
      na.rm = TRUE
    )
  ),
  by = split
]
fwrite(combined_summary, file.path(output_dir, "combined_performance_summary.csv"))
fwrite(classification_long, file.path(output_dir, "classification_metrics_long.csv"))
fwrite(error_long, file.path(output_dir, "error_metrics_long.csv"))

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

classification_plot <- ggplot(
  classification_long,
  aes(x = split, y = value, fill = kingdom)
) +
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

error_plot <- ggplot(
  error_long,
  aes(x = split, y = value, fill = kingdom)
) +
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

test_focus <- metrics[split == "test"]
if (nrow(test_focus) > 0) {
  test_tradeoff <- ggplot(
    test_focus,
    aes(
      x = presence_f1,
      y = mean_bray_curtis,
      size = n_microcosms,
      label = kingdom,
      color = kingdom
    )
  ) +
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
  histories[, kingdom := as.character(kingdom)]

  loss_plot <- ggplot(histories, aes(x = epoch, y = train_loss, color = kingdom)) +
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

  validation_long <- melt(
    histories,
    id.vars = c("kingdom", "epoch"),
    measure.vars = c(
      "validation_presence_f1",
      "validation_mean_bray_curtis",
      "validation_rmse_log_abundance_present_taxa"
    ),
    variable.name = "metric",
    value.name = "value"
  )
  validation_long[, metric_label := fifelse(
    metric == "validation_presence_f1",
    "Validation presence F1",
    fifelse(
      metric == "validation_mean_bray_curtis",
      "Validation mean Bray-Curtis",
      "Validation RMSE log abundance"
    )
  )]

  validation_plot <- ggplot(validation_long, aes(x = epoch, y = value, color = kingdom)) +
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
  graph_long <- melt(
    graph_summaries,
    id.vars = c("kingdom", "selected_edge_community_types", "edge_weight_transform", "edge_min_weight"),
    measure.vars = c(
      "n_edges_after_taxon_filter",
      "n_positive_edges",
      "n_negative_edges",
      "n_taxa_with_positive_neighbors",
      "n_taxa_with_negative_neighbors"
    ),
    variable.name = "metric",
    value.name = "value"
  )
  graph_long[, metric_label := fifelse(
    metric == "n_edges_after_taxon_filter",
    "Edges retained",
    fifelse(
      metric == "n_positive_edges",
      "Positive edges",
      fifelse(
        metric == "n_negative_edges",
        "Negative edges",
        fifelse(
          metric == "n_taxa_with_positive_neighbors",
          "Taxa with positive neighbors",
          "Taxa with negative neighbors"
        )
      )
    )
  )]
  fwrite(graph_long, file.path(output_dir, "spieceasi_graph_summary_long.csv"))

  graph_plot <- ggplot(graph_long, aes(x = kingdom, y = value, fill = kingdom)) +
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
  "- `test_f1_vs_bray_curtis_tradeoff.png/.pdf`: held-out F1 versus Bray-Curtis distance, when test metrics are available.",
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
  "- Training curves are diagnostic; test metrics are the clearest held-out performance summary."
)
writeLines(readme_lines, file.path(output_dir, "README.md"))

log_message("Done. Visualization outputs written to: ", output_dir)
