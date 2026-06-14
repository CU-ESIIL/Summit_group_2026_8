# =============================================================================
# coalescence_gnn_rstudio.R
# End-to-end GNN pipeline for microbial coalescence
# Runs in RStudio using reticulate 1.41+ (uv-managed Python)
#
# HOW TO USE:
#   1. Change GRAPH_INPUTS_DIR to your actual path
#   2. Make sure coalescence_gnn.py is in the same folder as this script
#   3. Run the entire script (Ctrl+Shift+Enter or Source button)
# =============================================================================

# ── 0. INSTALL R PACKAGES (run once, then comment out) ───────────────────────
# install.packages("reticulate")
# install.packages("ggplot2")
# install.packages("tidyr")

# ── 1. DECLARE PYTHON REQUIREMENTS ───────────────────────────────────────────
# Must come before any import() or py_module_available() calls
# reticulate 1.41+ uses uv to auto-manage the Python environment

library(reticulate)

py_require(python_version = "3.11")
py_require("tensorflow==2.15.*")
py_require("numpy")
py_require("pandas")
py_require("scikit-learn")

# ── 2. VERIFY TENSORFLOW ──────────────────────────────────────────────────────
cat("=============================================================\n")
cat("Verifying Python environment...\n")
cat("Python path:", py_config()$python, "\n")

if (!py_module_available("tensorflow")) {
  stop(paste(
    "TensorFlow not found.",
    "Run reticulate::py_last_error() for details.",
    "Try: py_require('tensorflow-cpu>=2.13,<2.16') as a fallback."
  ))
}

tf <- import("tensorflow")
cat("TensorFlow version:", tf$`__version__`, "\n")
cat("=============================================================\n\n")

# ── 3. SET FILE PATHS ─────────────────────────────────────────────────────────
# !! CHANGE THIS TO YOUR ACTUAL DIRECTORY !!
GRAPH_INPUTS_DIR <- "/home/jovyan/data-store/"

# Output directory for results
OUTPUT_DIR <- "gnn_outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Pass paths to Python namespace
py$NODES_SAMPLES_PATH   <- file.path(GRAPH_INPUTS_DIR, "nodes_samples.csv")
py$NODES_TAXA_PATH      <- file.path(GRAPH_INPUTS_DIR, "nodes_taxa.csv")
py$EDGES_SAMPLE_TAXON   <- file.path(GRAPH_INPUTS_DIR, "combined_sample_taxon_edges.csv")
py$EDGES_TAXON_TAXON    <- file.path(GRAPH_INPUTS_DIR, "taxon_taxon_spieceasi_edges.csv")
py$TRIPLETS_PATH        <- file.path(GRAPH_INPUTS_DIR, "coalescence_triplets.csv")
py$OUTPUT_DIR           <- OUTPUT_DIR

cat("Input paths set:\n")
cat("  nodes_samples:  ", py$NODES_SAMPLES_PATH, "\n")
cat("  nodes_taxa:     ", py$NODES_TAXA_PATH, "\n")
cat("  edges_st:       ", py$EDGES_SAMPLE_TAXON, "\n")
cat("  edges_tt:       ", py$EDGES_TAXON_TAXON, "\n")
cat("  triplets:       ", py$TRIPLETS_PATH, "\n")
cat("  output_dir:     ", py$OUTPUT_DIR, "\n\n")

# Verify input files exist
input_files <- c(
  py$NODES_SAMPLES_PATH,
  py$NODES_TAXA_PATH,
  py$EDGES_SAMPLE_TAXON,
  py$TRIPLETS_PATH
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(" -", missing_files, collapse = "\n"))
}
cat("All required input files found.\n\n")

# ── 4. RUN THE GNN PIPELINE ───────────────────────────────────────────────────
cat("=============================================================\n")
cat("Launching GNN pipeline...\n")
cat("=============================================================\n\n")

# Get the directory of this R script so we can find the .py file
script_dir <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) getwd()
)
py_script_path <- file.path(script_dir, "coalescence_gnn.py")

if (!file.exists(py_script_path)) {
  stop(
    "Cannot find coalescence_gnn.py at: ", py_script_path,
    "\nMake sure coalescence_gnn.py is in the same folder as this R script."
  )
}

source_python(py_script_path)

# ── 5. RETRIEVE AND DISPLAY RESULTS ──────────────────────────────────────────
cat("\n=============================================================\n")
cat("Results\n")
cat("=============================================================\n")
cat(sprintf("Test KL Divergence:             %.4f\n", py$test_kl))
cat(sprintf("Test Bray-Curtis Dissimilarity: %.4f\n", py$test_bc))
cat("  (Bray-Curtis: 0 = identical communities, 1 = completely different)\n\n")

# Load result CSVs into R
predictions_r <- read.csv(file.path(OUTPUT_DIR, "gnn_predictions_test.csv"))
history_r     <- read.csv(file.path(OUTPUT_DIR, "gnn_training_history.csv"))

cat("Prediction output shape:", nrow(predictions_r), "samples x",
    ncol(predictions_r) - 1, "taxa\n")
cat("Training ran for", nrow(history_r), "epochs\n\n")

# ── 6. PLOT TRAINING HISTORY ──────────────────────────────────────────────────
cat("Generating training plots...\n")

library(ggplot2)

history_r$epoch <- seq_len(nrow(history_r))

# KL Divergence loss plot
p_loss <- ggplot(history_r) +
  geom_line(aes(x = epoch, y = train_loss, colour = "Train"),
            linewidth = 1.0) +
  geom_line(aes(x = epoch, y = val_loss,   colour = "Validation"),
            linewidth = 1.0) +
  scale_colour_manual(values = c("Train" = "#2196F3", "Validation" = "#F44336")) +
  labs(
    title    = "CoalescenceGNN — Training History",
    subtitle = "KL Divergence loss (lower = better)",
    x        = "Epoch",
    y        = "KL Divergence",
    colour   = "Split"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# Bray-Curtis validation metric plot
p_bc <- ggplot(history_r, aes(x = epoch, y = val_bc)) +
  geom_line(colour = "#4CAF50", linewidth = 1.0) +
  geom_hline(yintercept = min(history_r$val_bc),
             linetype = "dashed", colour = "grey50") +
  annotate("text",
           x     = which.min(history_r$val_bc),
           y     = min(history_r$val_bc) + 0.005,
           label = sprintf("Best: %.4f", min(history_r$val_bc)),
           size  = 3.5, colour = "grey30") +
  labs(
    title    = "CoalescenceGNN — Validation Bray-Curtis",
    subtitle = "Community dissimilarity (lower = better)",
    x        = "Epoch",
    y        = "Bray-Curtis Dissimilarity"
  ) +
  theme_minimal(base_size = 13)

print(p_loss)
print(p_bc)

ggsave(file.path(OUTPUT_DIR, "training_loss.png"),    p_loss, width = 8, height = 5, dpi = 150)
ggsave(file.path(OUTPUT_DIR, "validation_bray_curtis.png"), p_bc,  width = 8, height = 5, dpi = 150)

# ── 7. SUMMARY TABLE ─────────────────────────────────────────────────────────
cat("\n=============================================================\n")
cat("Output files written to:", OUTPUT_DIR, "\n")
cat("=============================================================\n")
cat("  gnn_predictions_test.csv       — predicted final community abundances\n")
cat("  gnn_training_history.csv       — loss and metric per epoch\n")
cat("  best_coalescence_gnn.weights.h5 — saved model weights\n")
cat("  training_loss.png              — KL divergence training curve\n")
cat("  validation_bray_curtis.png     — Bray-Curtis validation curve\n")
cat("\nPipeline complete.\n")