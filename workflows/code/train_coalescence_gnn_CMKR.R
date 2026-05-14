#!/usr/bin/env Rscript

# Train a taxon-node graph neural network for microbial coalescence.
#
# The supervised unit is:
#   resident sample + pooled donor source -> observed final community
#
# SpiecEasi taxon-taxon edges are used as a signed, weighted association
# scaffold. Positive and negative associations are passed through separate graph
# channels, while donor/resident community composition remains the supervised
# input.
#
# Run from the repository root:
#   Rscript workflows/code/train_coalescence_gnn_CMKR.R
#
# Example local smoke run:
#   Rscript workflows/code/train_coalescence_gnn_CMKR.R --kingdoms=Bacteria --max_taxa_per_kingdom=500 --epochs=5
#
# Example larger HPC run:
#   Rscript workflows/code/train_coalescence_gnn_CMKR.R --max_taxa_per_kingdom=Inf --epochs=200 --batch_size=16

required_packages <- c("data.table", "torch")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running this workflow. ",
      "For example: install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(torch)
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

parse_bool <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

parse_csv <- function(x) {
  values <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

parse_number_or_inf <- function(x) {
  x <- as.character(x)
  if (tolower(x) %in% c("inf", "infinity", "all", "none")) {
    Inf
  } else {
    as.numeric(x)
  }
}

log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, paste0(..., collapse = "")))
}

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

config <- list(
  graph_input_dir = arg_value(
    cli_args,
    "graph_input_dir",
    file.path(repo_root, "workflows", "output", "graph_inputs")
  ),
  output_dir = arg_value(
    cli_args,
    "output_dir",
    file.path(repo_root, "workflows", "output", "gnn_model")
  ),
  kingdoms = parse_csv(arg_value(cli_args, "kingdoms", Sys.getenv("GNN_KINGDOMS", "Bacteria,Fungi"))),
  max_taxa_per_kingdom = parse_number_or_inf(arg_value(
    cli_args,
    "max_taxa_per_kingdom",
    Sys.getenv("GNN_MAX_TAXA_PER_KINGDOM", "2000")
  )),
  epochs = as.integer(arg_value(cli_args, "epochs", Sys.getenv("GNN_EPOCHS", "50"))),
  batch_size = as.integer(arg_value(cli_args, "batch_size", Sys.getenv("GNN_BATCH_SIZE", "8"))),
  hidden_dim = as.integer(arg_value(cli_args, "hidden_dim", Sys.getenv("GNN_HIDDEN_DIM", "128"))),
  taxon_embedding_dim = as.integer(arg_value(
    cli_args,
    "taxon_embedding_dim",
    Sys.getenv("GNN_TAXON_EMBEDDING_DIM", "16")
  )),
  dropout = as.numeric(arg_value(cli_args, "dropout", Sys.getenv("GNN_DROPOUT", "0.20"))),
  learning_rate = as.numeric(arg_value(cli_args, "learning_rate", Sys.getenv("GNN_LEARNING_RATE", "0.001"))),
  weight_decay = as.numeric(arg_value(cli_args, "weight_decay", Sys.getenv("GNN_WEIGHT_DECAY", "0.0001"))),
  train_fraction = as.numeric(arg_value(cli_args, "train_fraction", Sys.getenv("GNN_TRAIN_FRACTION", "0.70"))),
  validation_fraction = as.numeric(arg_value(
    cli_args,
    "validation_fraction",
    Sys.getenv("GNN_VALIDATION_FRACTION", "0.15")
  )),
  split_by = arg_value(cli_args, "split_by", Sys.getenv("GNN_SPLIT_BY", "microcosm")),
  edge_community_types = parse_csv(arg_value(
    cli_args,
    "edge_community_types",
    Sys.getenv("GNN_EDGE_COMMUNITY_TYPES", "resident,donor")
  )),
  edge_min_weight = as.numeric(arg_value(
    cli_args,
    "edge_min_weight",
    Sys.getenv("GNN_EDGE_MIN_WEIGHT", "0")
  )),
  edge_weight_transform = arg_value(
    cli_args,
    "edge_weight_transform",
    Sys.getenv("GNN_EDGE_WEIGHT_TRANSFORM", "sqrt")
  ),
  abundance_log_scale = as.numeric(arg_value(
    cli_args,
    "abundance_log_scale",
    Sys.getenv("GNN_ABUNDANCE_LOG_SCALE", "10000")
  )),
  presence_loss_weight = as.numeric(arg_value(
    cli_args,
    "presence_loss_weight",
    Sys.getenv("GNN_PRESENCE_LOSS_WEIGHT", "1.0")
  )),
  abundance_loss_weight = as.numeric(arg_value(
    cli_args,
    "abundance_loss_weight",
    Sys.getenv("GNN_ABUNDANCE_LOSS_WEIGHT", "0.5")
  )),
  prediction_threshold = as.numeric(arg_value(
    cli_args,
    "prediction_threshold",
    Sys.getenv("GNN_PREDICTION_THRESHOLD", "0.5")
  )),
  random_seed = as.integer(arg_value(cli_args, "random_seed", Sys.getenv("GNN_RANDOM_SEED", "1"))),
  use_cuda = parse_bool(arg_value(cli_args, "use_cuda", Sys.getenv("GNN_USE_CUDA", "false"))),
  write_predictions = parse_bool(arg_value(
    cli_args,
    "write_predictions",
    Sys.getenv("GNN_WRITE_PREDICTIONS", "true")
  ))
)

sample_taxon_path <- file.path(config$graph_input_dir, "combined_sample_taxon_edges.csv")
triplets_path <- file.path(config$graph_input_dir, "coalescence_triplets.csv")
taxa_path <- file.path(config$graph_input_dir, "nodes_taxa.csv")
spieceasi_edges_path <- file.path(config$graph_input_dir, "taxon_taxon_spieceasi_edges.csv")

for (path in c(sample_taxon_path, triplets_path, taxa_path, spieceasi_edges_path)) {
  if (!file.exists(path)) {
    stop(sprintf("Required input file not found: %s", path), call. = FALSE)
  }
}

if (!dir.exists(config$output_dir)) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
}

set.seed(config$random_seed)
torch_manual_seed(config$random_seed)

if (config$use_cuda && cuda_is_available()) {
  device <- torch_device("cuda")
} else {
  device <- torch_device("cpu")
}

log_message("Using device: ", as.character(device))
log_message("Reading triplets: ", triplets_path)
triplets <- fread(triplets_path, na.strings = c("", "NA"))

log_message("Reading taxon metadata: ", taxa_path)
taxa <- fread(taxa_path, na.strings = c("", "NA"))

log_message("Reading sample-taxon edges: ", sample_taxon_path)
sample_taxon_edges <- fread(
  sample_taxon_path,
  select = c("sample_id", "taxon_id", "abundance", "kingdom", "donor_id", "community_type"),
  na.strings = c("", "NA")
)
sample_taxon_edges[, abundance := fifelse(is.na(abundance), 0, as.numeric(abundance))]

log_message("Reading SpiecEasi taxon-taxon edges: ", spieceasi_edges_path)
spieceasi_edges <- fread(
  spieceasi_edges_path,
  select = c("kingdom", "donor_id", "community_type", "taxon_a", "taxon_b", "weight", "sign", "method"),
  na.strings = c("", "NA")
)
spieceasi_edges[, weight := fifelse(is.na(weight), 0, as.numeric(weight))]
spieceasi_edges[, sign := fifelse(is.na(sign) | sign == 0, 1, as.numeric(sign))]

sample_taxon_edges <- sample_taxon_edges[kingdom %in% config$kingdoms]
spieceasi_edges <- spieceasi_edges[kingdom %in% config$kingdoms]
triplets <- triplets[kingdom %in% config$kingdoms]
taxa <- taxa[kingdom %in% config$kingdoms]

if (nrow(triplets) == 0) {
  stop("No coalescence triplets remain after kingdom filtering.", call. = FALSE)
}

feature_names <- c(
  "resident_abundance",
  "donor_abundance",
  "resident_present",
  "donor_present",
  "resident_log1p_scaled",
  "donor_log1p_scaled",
  "donor_minus_resident"
)

select_taxa_for_kingdom <- function(edges, taxa, kingdom_name, max_taxa) {
  edges_k <- edges[kingdom == kingdom_name]
  taxa_k <- unique(taxa[kingdom == kingdom_name, .(kingdom, taxon_id, name, genus, species, phylum, class, order, family)])

  pre_stats <- edges_k[
    community_type %in% c("donor", "resident"),
    .(
      pre_abundance_sum = sum(abundance, na.rm = TRUE),
      pre_nonzero_count = sum(abundance > 0, na.rm = TRUE)
    ),
    by = taxon_id
  ]

  if (nrow(pre_stats) == 0) {
    stop(sprintf("No donor/resident taxa found for kingdom %s.", kingdom_name), call. = FALSE)
  }

  setorder(pre_stats, -pre_abundance_sum, -pre_nonzero_count, taxon_id)

  if (is.finite(max_taxa)) {
    selected_taxa <- head(pre_stats$taxon_id, as.integer(max_taxa))
  } else {
    selected_taxa <- pre_stats$taxon_id
  }

  selected <- data.table(kingdom = kingdom_name, taxon_id = selected_taxa)
  selected <- merge(selected, pre_stats, by = "taxon_id", all.x = TRUE, sort = FALSE)
  selected <- merge(selected, taxa_k, by = c("kingdom", "taxon_id"), all.x = TRUE, sort = FALSE)
  selected[, taxon_index := seq_len(.N)]
  selected[]
}

transform_edge_weight <- function(weight, transform_name) {
  transform_name <- tolower(transform_name)
  if (transform_name == "none") {
    return(weight)
  }
  if (transform_name == "sqrt") {
    return(sqrt(weight))
  }
  if (transform_name == "log1p") {
    return(log1p(weight))
  }
  stop(
    sprintf("Unknown edge_weight_transform '%s'. Use none, sqrt, or log1p.", transform_name),
    call. = FALSE
  )
}

normalize_adjacency_rows <- function(adj) {
  row_sums <- rowSums(adj)
  nonzero_rows <- which(row_sums > 0)
  if (length(nonzero_rows) > 0) {
    adj[nonzero_rows, ] <- adj[nonzero_rows, , drop = FALSE] / row_sums[nonzero_rows]
  }
  adj
}

build_spieceasi_adjacency <- function(spieceasi_edges, selected_taxa, kingdom_name, config) {
  n_taxa <- nrow(selected_taxa)
  taxon_lookup <- setNames(selected_taxa$taxon_index, selected_taxa$taxon_id)

  edge_dt <- copy(spieceasi_edges[kingdom == kingdom_name])
  if (!("all" %in% tolower(config$edge_community_types))) {
    edge_dt <- edge_dt[community_type %in% config$edge_community_types]
  }

  n_edges_input <- nrow(edge_dt)
  if (n_edges_input == 0) {
    summary <- data.table(
      kingdom = kingdom_name,
      selected_edge_community_types = paste(config$edge_community_types, collapse = ","),
      edge_weight_transform = config$edge_weight_transform,
      edge_min_weight = config$edge_min_weight,
      n_edges_input = 0,
      n_edges_after_taxon_filter = 0,
      n_positive_edges = 0,
      n_negative_edges = 0,
      n_taxa_with_positive_neighbors = 0,
      n_taxa_with_negative_neighbors = 0
    )
    return(list(
      adj_pos = matrix(0, nrow = n_taxa, ncol = n_taxa),
      adj_neg = matrix(0, nrow = n_taxa, ncol = n_taxa),
      edges_used = data.table(),
      summary = summary
    ))
  }

  edge_dt[, from_index := unname(taxon_lookup[taxon_a])]
  edge_dt[, to_index := unname(taxon_lookup[taxon_b])]
  edge_dt <- edge_dt[
    !is.na(from_index) &
      !is.na(to_index) &
      from_index != to_index &
      !is.na(weight) &
      abs(weight) >= config$edge_min_weight
  ]

  if (nrow(edge_dt) == 0) {
    summary <- data.table(
      kingdom = kingdom_name,
      selected_edge_community_types = paste(config$edge_community_types, collapse = ","),
      edge_weight_transform = config$edge_weight_transform,
      edge_min_weight = config$edge_min_weight,
      n_edges_input = n_edges_input,
      n_edges_after_taxon_filter = 0,
      n_positive_edges = 0,
      n_negative_edges = 0,
      n_taxa_with_positive_neighbors = 0,
      n_taxa_with_negative_neighbors = 0
    )
    return(list(
      adj_pos = matrix(0, nrow = n_taxa, ncol = n_taxa),
      adj_neg = matrix(0, nrow = n_taxa, ncol = n_taxa),
      edges_used = data.table(),
      summary = summary
    ))
  }

  edge_dt[, edge_strength := transform_edge_weight(abs(weight), config$edge_weight_transform)]
  edge_dt[, signed_strength := edge_strength * fifelse(sign >= 0, 1, -1)]
  edge_dt[, edge_i := pmin(from_index, to_index)]
  edge_dt[, edge_j := pmax(from_index, to_index)]

  edge_agg <- edge_dt[
    ,
    .(
      signed_strength = mean(signed_strength, na.rm = TRUE),
      mean_abs_strength = mean(edge_strength, na.rm = TRUE),
      n_edge_contexts = .N,
      n_donor_contexts = uniqueN(donor_id),
      community_contexts = paste(sort(unique(community_type)), collapse = ";"),
      methods = paste(sort(unique(method)), collapse = ";")
    ),
    by = .(from_index = edge_i, to_index = edge_j)
  ]
  edge_agg[, from_taxon := selected_taxa$taxon_id[from_index]]
  edge_agg[, to_taxon := selected_taxa$taxon_id[to_index]]
  edge_agg[, edge_sign := fifelse(signed_strength >= 0, 1, -1)]
  edge_agg[, edge_weight := abs(signed_strength)]

  adj_pos <- matrix(0, nrow = n_taxa, ncol = n_taxa)
  adj_neg <- matrix(0, nrow = n_taxa, ncol = n_taxa)

  for (row_idx in seq_len(nrow(edge_agg))) {
    i <- edge_agg$from_index[[row_idx]]
    j <- edge_agg$to_index[[row_idx]]
    signed_weight <- edge_agg$signed_strength[[row_idx]]
    if (is.na(signed_weight) || signed_weight == 0) {
      next
    }
    if (signed_weight > 0) {
      adj_pos[i, j] <- abs(signed_weight)
      adj_pos[j, i] <- abs(signed_weight)
    } else {
      adj_neg[i, j] <- abs(signed_weight)
      adj_neg[j, i] <- abs(signed_weight)
    }
  }

  pos_neighbor_counts <- rowSums(adj_pos > 0)
  neg_neighbor_counts <- rowSums(adj_neg > 0)
  adj_pos <- normalize_adjacency_rows(adj_pos)
  adj_neg <- normalize_adjacency_rows(adj_neg)

  summary <- data.table(
    kingdom = kingdom_name,
    selected_edge_community_types = paste(config$edge_community_types, collapse = ","),
    edge_weight_transform = config$edge_weight_transform,
    edge_min_weight = config$edge_min_weight,
    n_edges_input = n_edges_input,
    n_edges_after_taxon_filter = nrow(edge_agg),
    n_positive_edges = sum(edge_agg$signed_strength > 0, na.rm = TRUE),
    n_negative_edges = sum(edge_agg$signed_strength < 0, na.rm = TRUE),
    n_taxa_with_positive_neighbors = sum(pos_neighbor_counts > 0),
    n_taxa_with_negative_neighbors = sum(neg_neighbor_counts > 0)
  )

  list(
    adj_pos = adj_pos,
    adj_neg = adj_neg,
    edges_used = edge_agg[
      ,
      .(
        kingdom = kingdom_name,
        from_taxon,
        to_taxon,
        edge_weight,
        edge_sign,
        signed_strength,
        mean_abs_strength,
        n_edge_contexts,
        n_donor_contexts,
        community_contexts,
        methods
      )
    ],
    summary = summary
  )
}

split_examples <- function(triplets_k, config) {
  n_examples <- nrow(triplets_k)
  if (n_examples < 3) {
    stop("Need at least three prediction units to create train/validation/test splits.", call. = FALSE)
  }

  split_by <- match.arg(config$split_by, c("microcosm", "donor_id"))

  if (split_by == "donor_id") {
    units <- unique(triplets_k$donor_id)
    if (length(units) < 3) {
      warning("Too few donor_id groups for donor-level split; falling back to microcosm-level split.")
      split_by <- "microcosm"
    }
  }

  if (split_by == "microcosm") {
    unit_dt <- data.table(unit = triplets_k$microcosm_id)
  } else {
    unit_dt <- data.table(unit = unique(triplets_k$donor_id))
  }

  unit_dt <- unit_dt[sample(.N)]
  n_units <- nrow(unit_dt)
  n_train <- max(1, floor(config$train_fraction * n_units))
  n_validation <- max(1, floor(config$validation_fraction * n_units))

  if (n_train + n_validation > n_units - 1) {
    n_train <- max(1, n_units - n_validation - 1)
  }
  if (n_train + n_validation > n_units - 1) {
    n_validation <- max(1, n_units - n_train - 1)
  }

  train_units <- unit_dt$unit[seq_len(n_train)]
  validation_units <- unit_dt$unit[seq.int(n_train + 1, n_train + n_validation)]
  test_units <- setdiff(unit_dt$unit, c(train_units, validation_units))

  if (split_by == "microcosm") {
    train_idx <- which(triplets_k$microcosm_id %in% train_units)
    validation_idx <- which(triplets_k$microcosm_id %in% validation_units)
    test_idx <- which(triplets_k$microcosm_id %in% test_units)
  } else {
    train_idx <- which(triplets_k$donor_id %in% train_units)
    validation_idx <- which(triplets_k$donor_id %in% validation_units)
    test_idx <- which(triplets_k$donor_id %in% test_units)
  }

  list(
    train = train_idx,
    validation = validation_idx,
    test = test_idx
  )
}

fill_array_feature <- function(x, row_idx, col_idx, feature_idx, values) {
  valid <- !is.na(row_idx) & !is.na(col_idx)
  if (any(valid)) {
    x[cbind(row_idx[valid], col_idx[valid], feature_idx)] <- values[valid]
  }
  x
}

fill_matrix_target <- function(y, row_idx, col_idx, values) {
  valid <- !is.na(row_idx) & !is.na(col_idx)
  if (any(valid)) {
    y[cbind(row_idx[valid], col_idx[valid])] <- values[valid]
  }
  y
}

build_kingdom_dataset <- function(edges, triplets, taxa, spieceasi_edges, kingdom_name, config) {
  log_message("Building dataset for ", kingdom_name)

  triplets_k <- copy(triplets[kingdom == kingdom_name])
  if (nrow(triplets_k) == 0) {
    return(NULL)
  }
  triplets_k[, example_index := seq_len(.N)]

  selected_taxa <- select_taxa_for_kingdom(edges, taxa, kingdom_name, config$max_taxa_per_kingdom)
  taxon_lookup <- setNames(selected_taxa$taxon_index, selected_taxa$taxon_id)
  graph_scaffold <- build_spieceasi_adjacency(spieceasi_edges, selected_taxa, kingdom_name, config)
  log_message(
    kingdom_name,
    " SpiecEasi scaffold: ",
    graph_scaffold$summary$n_edges_after_taxon_filter,
    " aggregated edges across community_type=",
    graph_scaffold$summary$selected_edge_community_types
  )

  n_examples <- nrow(triplets_k)
  n_taxa <- nrow(selected_taxa)
  n_features <- length(feature_names)

  x_raw <- array(0, dim = c(n_examples, n_taxa, n_features))
  y_presence <- matrix(0, nrow = n_examples, ncol = n_taxa)
  y_log_abundance <- matrix(0, nrow = n_examples, ncol = n_taxa)
  y_abundance <- matrix(0, nrow = n_examples, ncol = n_taxa)

  edges_k <- edges[kingdom == kingdom_name & taxon_id %in% selected_taxa$taxon_id]

  resident_edges <- edges_k[
    community_type == "resident",
    .(resident_abundance = sum(abundance, na.rm = TRUE)),
    by = .(donor_id, sample_id, taxon_id)
  ]
  resident_join <- merge(
    triplets_k[, .(example_index, donor_id, resident_sample_id)],
    resident_edges,
    by.x = c("donor_id", "resident_sample_id"),
    by.y = c("donor_id", "sample_id"),
    allow.cartesian = TRUE
  )
  resident_join[, taxon_index := unname(taxon_lookup[taxon_id])]
  x_raw <- fill_array_feature(
    x_raw,
    resident_join$example_index,
    resident_join$taxon_index,
    1,
    resident_join$resident_abundance
  )

  donor_edges <- edges_k[
    community_type == "donor",
    .(donor_abundance = mean(abundance, na.rm = TRUE)),
    by = .(donor_id, taxon_id)
  ]
  donor_join <- merge(
    triplets_k[, .(example_index, donor_id)],
    donor_edges,
    by = "donor_id",
    allow.cartesian = TRUE
  )
  donor_join[, taxon_index := unname(taxon_lookup[taxon_id])]
  x_raw <- fill_array_feature(
    x_raw,
    donor_join$example_index,
    donor_join$taxon_index,
    2,
    donor_join$donor_abundance
  )

  final_edges <- edges_k[
    community_type == "final",
    .(final_abundance = sum(abundance, na.rm = TRUE)),
    by = .(donor_id, sample_id, taxon_id)
  ]
  final_join <- merge(
    triplets_k[, .(example_index, donor_id, final_sample_id)],
    final_edges,
    by.x = c("donor_id", "final_sample_id"),
    by.y = c("donor_id", "sample_id"),
    allow.cartesian = TRUE
  )
  final_join[, taxon_index := unname(taxon_lookup[taxon_id])]

  y_abundance <- fill_matrix_target(
    y_abundance,
    final_join$example_index,
    final_join$taxon_index,
    final_join$final_abundance
  )
  y_presence <- fill_matrix_target(
    y_presence,
    final_join$example_index,
    final_join$taxon_index,
    as.numeric(final_join$final_abundance > 0)
  )
  y_log_abundance <- fill_matrix_target(
    y_log_abundance,
    final_join$example_index,
    final_join$taxon_index,
    log1p(final_join$final_abundance * config$abundance_log_scale)
  )

  x_raw[, , 3] <- as.numeric(x_raw[, , 1] > 0)
  x_raw[, , 4] <- as.numeric(x_raw[, , 2] > 0)
  x_raw[, , 5] <- log1p(x_raw[, , 1] * config$abundance_log_scale)
  x_raw[, , 6] <- log1p(x_raw[, , 2] * config$abundance_log_scale)
  x_raw[, , 7] <- x_raw[, , 2] - x_raw[, , 1]

  splits <- split_examples(triplets_k, config)

  list(
    kingdom = kingdom_name,
    triplets = triplets_k,
    taxa = selected_taxa,
    x_raw = x_raw,
    y_presence = y_presence,
    y_log_abundance = y_log_abundance,
    y_abundance = y_abundance,
    adj_pos = graph_scaffold$adj_pos,
    adj_neg = graph_scaffold$adj_neg,
    graph_edges_used = graph_scaffold$edges_used,
    graph_edge_summary = graph_scaffold$summary,
    splits = splits,
    feature_names = feature_names
  )
}

standardize_features <- function(x, train_idx) {
  n_features <- dim(x)[3]
  train_flat <- array(x[train_idx, , , drop = FALSE], dim = c(length(train_idx) * dim(x)[2], n_features))
  means <- colMeans(train_flat, na.rm = TRUE)
  sds <- apply(train_flat, 2, stats::sd, na.rm = TRUE)
  sds[is.na(sds) | sds == 0] <- 1

  x_scaled <- x
  for (feature_idx in seq_len(n_features)) {
    x_scaled[, , feature_idx] <- (x_scaled[, , feature_idx] - means[[feature_idx]]) / sds[[feature_idx]]
  }

  list(
    x = x_scaled,
    means = means,
    sds = sds
  )
}

signed_graph_block <- nn_module(
  "signed_graph_block",
  initialize = function(hidden_dim, dropout) {
    self$self_linear <- nn_linear(hidden_dim, hidden_dim)
    self$positive_neighbor_linear <- nn_linear(hidden_dim, hidden_dim)
    self$negative_neighbor_linear <- nn_linear(hidden_dim, hidden_dim)
    self$norm <- nn_layer_norm(hidden_dim)
    self$dropout <- nn_dropout(dropout)
  },
  forward = function(h, adj_pos, adj_neg) {
    residual <- h
    batch_size <- h$size(1)
    adj_pos_batch <- adj_pos$unsqueeze(1)$expand(c(batch_size, adj_pos$size(1), adj_pos$size(2)))
    adj_neg_batch <- adj_neg$unsqueeze(1)$expand(c(batch_size, adj_neg$size(1), adj_neg$size(2)))

    positive_messages <- torch_bmm(adj_pos_batch, h)
    negative_messages <- torch_bmm(adj_neg_batch, h)

    h <- self$self_linear(h) +
      self$positive_neighbor_linear(positive_messages) +
      self$negative_neighbor_linear(negative_messages)
    h <- nnf_relu(h)
    h <- self$dropout(h)
    self$norm(h + residual)
  }
)

coalescence_gnn <- nn_module(
  "coalescence_gnn",
  initialize = function(n_features, n_taxa, embedding_dim, hidden_dim, dropout) {
    self$n_taxa <- n_taxa
    self$embedding_dim <- embedding_dim
    self$taxon_embedding <- nn_embedding(n_taxa, embedding_dim)
    self$input_projection <- nn_linear(n_features + embedding_dim, hidden_dim)
    self$block1 <- signed_graph_block(hidden_dim, dropout)
    self$block2 <- signed_graph_block(hidden_dim, dropout)
    self$presence_head <- nn_sequential(
      nn_linear(hidden_dim, hidden_dim),
      nn_relu(),
      nn_dropout(dropout),
      nn_linear(hidden_dim, 1)
    )
    self$abundance_head <- nn_sequential(
      nn_linear(hidden_dim, hidden_dim),
      nn_relu(),
      nn_dropout(dropout),
      nn_linear(hidden_dim, 1)
    )
  },
  forward = function(x, adj_pos, adj_neg) {
    batch_size <- x$size(1)
    taxon_ids <- torch_arange(
      start = 1,
      end = self$n_taxa,
      dtype = torch_long(),
      device = x$device
    )
    taxon_embeddings <- self$taxon_embedding(taxon_ids)
    taxon_embeddings <- taxon_embeddings$unsqueeze(1)$expand(c(batch_size, self$n_taxa, self$embedding_dim))

    h <- torch_cat(list(x, taxon_embeddings), dim = 3)
    h <- self$input_projection(h)
    h <- nnf_relu(h)
    h <- self$block1(h, adj_pos, adj_neg)
    h <- self$block2(h, adj_pos, adj_neg)

    presence_logits <- self$presence_head(h)$squeeze(3)
    log_abundance <- nnf_softplus(self$abundance_head(h)$squeeze(3))

    list(
      presence_logits = presence_logits,
      log_abundance = log_abundance
    )
  }
)

make_batches <- function(indices, batch_size, shuffle = TRUE) {
  if (shuffle) {
    indices <- sample(indices)
  }
  split(indices, ceiling(seq_along(indices) / batch_size))
}

compute_loss <- function(output, y_presence, y_log_abundance, pos_weight_tensor, config) {
  presence_loss <- nnf_binary_cross_entropy_with_logits(
    output$presence_logits,
    y_presence,
    pos_weight = pos_weight_tensor
  )

  present_mask <- y_presence
  present_count <- as.numeric(present_mask$sum()$item())
  if (present_count > 0) {
    abundance_diff <- (output$log_abundance - y_log_abundance) * present_mask
    abundance_loss <- abundance_diff$pow(2)$sum() / present_count
  } else {
    abundance_loss <- output$log_abundance$sum() * 0
  }

  total_loss <- config$presence_loss_weight * presence_loss +
    config$abundance_loss_weight * abundance_loss

  list(
    total = total_loss,
    presence = presence_loss,
    abundance = abundance_loss
  )
}

clone_state_dict <- function(state_dict) {
  lapply(state_dict, function(tensor) tensor$clone())
}

predict_for_indices <- function(model, x_tensor, adj_pos_tensor, adj_neg_tensor, indices, batch_size, device) {
  model$eval()
  presence_chunks <- list()
  abundance_chunks <- list()

  with_no_grad({
    batches <- make_batches(indices, batch_size, shuffle = FALSE)
    for (batch_id in seq_along(batches)) {
      batch_idx <- batches[[batch_id]]
      xb <- x_tensor[batch_idx, , ]$to(device = device)
      output <- model(xb, adj_pos_tensor, adj_neg_tensor)
      presence_chunks[[batch_id]] <- as_array(torch_sigmoid(output$presence_logits)$to(device = "cpu"))
      abundance_chunks[[batch_id]] <- as_array(output$log_abundance$to(device = "cpu"))
    }
  })

  list(
    presence_probability = do.call(rbind, presence_chunks),
    predicted_log_abundance = do.call(rbind, abundance_chunks)
  )
}

compute_metrics <- function(prediction, dataset, indices, split_name, config) {
  observed_presence <- dataset$y_presence[indices, , drop = FALSE]
  observed_log_abundance <- dataset$y_log_abundance[indices, , drop = FALSE]
  observed_abundance <- dataset$y_abundance[indices, , drop = FALSE]

  predicted_presence <- prediction$presence_probability >= config$prediction_threshold
  observed_presence_bool <- observed_presence >= 0.5

  tp <- sum(predicted_presence & observed_presence_bool)
  tn <- sum(!predicted_presence & !observed_presence_bool)
  fp <- sum(predicted_presence & !observed_presence_bool)
  fn <- sum(!predicted_presence & observed_presence_bool)

  accuracy <- (tp + tn) / max(1, tp + tn + fp + fn)
  precision <- tp / max(1, tp + fp)
  recall <- tp / max(1, tp + fn)
  specificity <- tn / max(1, tn + fp)
  f1 <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)

  present_mask <- observed_presence_bool
  if (any(present_mask)) {
    rmse_log_present <- sqrt(mean((prediction$predicted_log_abundance[present_mask] - observed_log_abundance[present_mask])^2))
    mae_log_present <- mean(abs(prediction$predicted_log_abundance[present_mask] - observed_log_abundance[present_mask]))
  } else {
    rmse_log_present <- NA_real_
    mae_log_present <- NA_real_
  }

  predicted_abundance <- pmax(expm1(prediction$predicted_log_abundance) / config$abundance_log_scale, 0)
  predicted_abundance[!predicted_presence] <- 0

  bray_curtis <- vapply(
    seq_len(nrow(observed_abundance)),
    function(i) {
      denominator <- sum(predicted_abundance[i, ] + observed_abundance[i, ], na.rm = TRUE)
      if (denominator == 0) {
        return(NA_real_)
      }
      sum(abs(predicted_abundance[i, ] - observed_abundance[i, ]), na.rm = TRUE) / denominator
    },
    numeric(1)
  )

  data.table(
    kingdom = dataset$kingdom,
    split = split_name,
    n_microcosms = length(indices),
    n_taxa = nrow(dataset$taxa),
    presence_accuracy = accuracy,
    presence_precision = precision,
    presence_recall = recall,
    presence_specificity = specificity,
    presence_f1 = f1,
    rmse_log_abundance_present_taxa = rmse_log_present,
    mae_log_abundance_present_taxa = mae_log_present,
    mean_bray_curtis = mean(bray_curtis, na.rm = TRUE),
    median_bray_curtis = stats::median(bray_curtis, na.rm = TRUE)
  )
}

write_predictions <- function(prediction, dataset, split_labels, output_path, config) {
  n_examples <- nrow(dataset$triplets)
  n_taxa <- nrow(dataset$taxa)
  predicted_presence <- prediction$presence_probability >= config$prediction_threshold
  predicted_abundance <- pmax(expm1(prediction$predicted_log_abundance) / config$abundance_log_scale, 0)
  predicted_abundance[!predicted_presence] <- 0

  output <- data.table(
    kingdom = dataset$kingdom,
    microcosm_id = rep(dataset$triplets$microcosm_id, each = n_taxa),
    donor_id = rep(dataset$triplets$donor_id, each = n_taxa),
    resident_sample_id = rep(dataset$triplets$resident_sample_id, each = n_taxa),
    final_sample_id = rep(dataset$triplets$final_sample_id, each = n_taxa),
    split = rep(split_labels, each = n_taxa),
    taxon_id = rep(dataset$taxa$taxon_id, times = n_examples),
    resident_abundance = as.vector(t(dataset$x_raw[, , 1])),
    donor_abundance = as.vector(t(dataset$x_raw[, , 2])),
    observed_final_presence = as.vector(t(dataset$y_presence)),
    observed_final_abundance = as.vector(t(dataset$y_abundance)),
    predicted_final_presence_probability = as.vector(t(prediction$presence_probability)),
    predicted_final_presence = as.vector(t(as.numeric(predicted_presence))),
    predicted_final_abundance = as.vector(t(predicted_abundance))
  )

  fwrite(output, output_path)
}

train_one_kingdom <- function(dataset, config, device) {
  log_message("Training model for ", dataset$kingdom)

  standardized <- standardize_features(dataset$x_raw, dataset$splits$train)
  x_tensor <- torch_tensor(standardized$x, dtype = torch_float())
  y_presence_tensor <- torch_tensor(dataset$y_presence, dtype = torch_float())
  y_log_abundance_tensor <- torch_tensor(dataset$y_log_abundance, dtype = torch_float())
  adj_pos_tensor <- torch_tensor(dataset$adj_pos, dtype = torch_float())$to(device = device)
  adj_neg_tensor <- torch_tensor(dataset$adj_neg, dtype = torch_float())$to(device = device)

  train_presence <- dataset$y_presence[dataset$splits$train, , drop = FALSE]
  positives <- sum(train_presence == 1)
  negatives <- sum(train_presence == 0)
  pos_weight <- if (positives > 0) min(negatives / positives, 50) else 1
  pos_weight_tensor <- torch_tensor(pos_weight, dtype = torch_float())$to(device = device)

  model <- coalescence_gnn(
    n_features = length(dataset$feature_names),
    n_taxa = nrow(dataset$taxa),
    embedding_dim = config$taxon_embedding_dim,
    hidden_dim = config$hidden_dim,
    dropout = config$dropout
  )$to(device = device)

  optimizer <- optim_adamw(
    model$parameters,
    lr = config$learning_rate,
    weight_decay = config$weight_decay
  )

  history <- list()
  best_validation_loss <- Inf
  best_state <- NULL

  for (epoch in seq_len(config$epochs)) {
    model$train()
    epoch_losses <- c()
    epoch_presence_losses <- c()
    epoch_abundance_losses <- c()

    batches <- make_batches(dataset$splits$train, config$batch_size, shuffle = TRUE)
    for (batch_idx in batches) {
      xb <- x_tensor[batch_idx, , ]$to(device = device)
      yp <- y_presence_tensor[batch_idx, ]$to(device = device)
      ya <- y_log_abundance_tensor[batch_idx, ]$to(device = device)

      optimizer$zero_grad()
      output <- model(xb, adj_pos_tensor, adj_neg_tensor)
      losses <- compute_loss(output, yp, ya, pos_weight_tensor, config)
      losses$total$backward()
      optimizer$step()

      epoch_losses <- c(epoch_losses, as.numeric(losses$total$item()))
      epoch_presence_losses <- c(epoch_presence_losses, as.numeric(losses$presence$item()))
      epoch_abundance_losses <- c(epoch_abundance_losses, as.numeric(losses$abundance$item()))
    }

    validation_prediction <- predict_for_indices(
      model,
      x_tensor,
      adj_pos_tensor,
      adj_neg_tensor,
      dataset$splits$validation,
      config$batch_size,
      device
    )
    validation_metrics <- compute_metrics(
      validation_prediction,
      dataset,
      dataset$splits$validation,
      "validation",
      config
    )

    validation_loss_proxy <- validation_metrics$rmse_log_abundance_present_taxa +
      (1 - validation_metrics$presence_f1)
    if (is.na(validation_loss_proxy)) {
      validation_loss_proxy <- 1 - validation_metrics$presence_f1
    }

    if (validation_loss_proxy < best_validation_loss) {
      best_validation_loss <- validation_loss_proxy
      best_state <- clone_state_dict(model$state_dict())
    }

    history[[epoch]] <- data.table(
      kingdom = dataset$kingdom,
      epoch = epoch,
      train_loss = mean(epoch_losses),
      train_presence_loss = mean(epoch_presence_losses),
      train_abundance_loss = mean(epoch_abundance_losses),
      validation_presence_f1 = validation_metrics$presence_f1,
      validation_mean_bray_curtis = validation_metrics$mean_bray_curtis,
      validation_rmse_log_abundance_present_taxa = validation_metrics$rmse_log_abundance_present_taxa
    )

    if (epoch == 1 || epoch %% 5 == 0 || epoch == config$epochs) {
      log_message(
        dataset$kingdom,
        " epoch ",
        epoch,
        "/",
        config$epochs,
        " train_loss=",
        sprintf("%.4f", mean(epoch_losses)),
        " val_f1=",
        sprintf("%.4f", validation_metrics$presence_f1),
        " val_bray=",
        sprintf("%.4f", validation_metrics$mean_bray_curtis)
      )
    }
  }

  if (!is.null(best_state)) {
    model$load_state_dict(best_state)
  }

  all_predictions <- predict_for_indices(
    model,
    x_tensor,
    adj_pos_tensor,
    adj_neg_tensor,
    seq_len(nrow(dataset$triplets)),
    config$batch_size,
    device
  )

  split_labels <- rep("unassigned", nrow(dataset$triplets))
  split_labels[dataset$splits$train] <- "train"
  split_labels[dataset$splits$validation] <- "validation"
  split_labels[dataset$splits$test] <- "test"

  metrics <- rbindlist(list(
    compute_metrics(
      list(
        presence_probability = all_predictions$presence_probability[dataset$splits$train, , drop = FALSE],
        predicted_log_abundance = all_predictions$predicted_log_abundance[dataset$splits$train, , drop = FALSE]
      ),
      dataset,
      dataset$splits$train,
      "train",
      config
    ),
    compute_metrics(
      list(
        presence_probability = all_predictions$presence_probability[dataset$splits$validation, , drop = FALSE],
        predicted_log_abundance = all_predictions$predicted_log_abundance[dataset$splits$validation, , drop = FALSE]
      ),
      dataset,
      dataset$splits$validation,
      "validation",
      config
    ),
    compute_metrics(
      list(
        presence_probability = all_predictions$presence_probability[dataset$splits$test, , drop = FALSE],
        predicted_log_abundance = all_predictions$predicted_log_abundance[dataset$splits$test, , drop = FALSE]
      ),
      dataset,
      dataset$splits$test,
      "test",
      config
    )
  ), fill = TRUE)

  kingdom_dir <- file.path(config$output_dir, dataset$kingdom)
  if (!dir.exists(kingdom_dir)) {
    dir.create(kingdom_dir, recursive = TRUE, showWarnings = FALSE)
  }

  fwrite(rbindlist(history), file.path(kingdom_dir, "training_history.csv"))
  fwrite(metrics, file.path(kingdom_dir, "evaluation_metrics.csv"))
  fwrite(dataset$taxa, file.path(kingdom_dir, "taxa_used.csv"))
  fwrite(dataset$graph_edge_summary, file.path(kingdom_dir, "spieceasi_graph_summary.csv"))
  fwrite(dataset$graph_edges_used, file.path(kingdom_dir, "spieceasi_edges_used.csv"))
  fwrite(
    data.table(
      microcosm_id = dataset$triplets$microcosm_id,
      kingdom = dataset$triplets$kingdom,
      donor_id = dataset$triplets$donor_id,
      resident_sample_id = dataset$triplets$resident_sample_id,
      final_sample_id = dataset$triplets$final_sample_id,
      split = split_labels
    ),
    file.path(kingdom_dir, "data_splits.csv")
  )

  torch_save(model$state_dict(), file.path(kingdom_dir, "model_state.pt"))

  if (config$write_predictions) {
    write_predictions(
      all_predictions,
      dataset,
      split_labels,
      file.path(kingdom_dir, "predictions_long.csv"),
      config
    )
  }

  fwrite(
    data.table(
      feature_name = dataset$feature_names,
      mean = standardized$means,
      sd = standardized$sds
    ),
    file.path(kingdom_dir, "feature_standardization.csv")
  )

  metrics
}

config_table <- data.table(
  parameter = names(config),
  value = vapply(config, function(value) paste(value, collapse = ","), character(1))
)
fwrite(config_table, file.path(config$output_dir, "run_config.csv"))

all_metrics <- list()
for (kingdom_name in config$kingdoms) {
  dataset <- build_kingdom_dataset(sample_taxon_edges, triplets, taxa, spieceasi_edges, kingdom_name, config)
  if (is.null(dataset)) {
    warning(sprintf("Skipping %s because no triplets were available.", kingdom_name))
    next
  }
  all_metrics[[kingdom_name]] <- train_one_kingdom(dataset, config, device)
}

if (length(all_metrics) == 0) {
  stop("No kingdom models were trained.", call. = FALSE)
}

combined_metrics <- rbindlist(all_metrics, fill = TRUE)
fwrite(combined_metrics, file.path(config$output_dir, "evaluation_metrics_all_kingdoms.csv"))

log_message("Done. Outputs written to: ", config$output_dir)
log_message("Key output: evaluation_metrics_all_kingdoms.csv")
