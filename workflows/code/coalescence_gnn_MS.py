# =============================================================================
# coalescence_gnn.py  — IMPROVED VERSION
# Microbial Coalescence GNN — TensorFlow / Keras
#
# Key improvements over v1:
#   - get_config() implemented (fixes serialization warning)
#   - Top-K taxon filtering (fixes p >> n problem)
#   - Learning rate scheduler
#   - Gradient clipping
#   - Naive baseline comparison
#   - Per-kingdom evaluation
#   - Fixed duplicate index warning in predictions CSV
# =============================================================================

import os
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.preprocessing import LabelEncoder, StandardScaler

# ── 0. PATHS (injected from R, fallback for standalone use) ───────────────────
try:
    NODES_SAMPLES_PATH
except NameError:
    GRAPH_INPUTS_DIR   = "workflows/output/graph_inputs"
    NODES_SAMPLES_PATH = os.path.join(GRAPH_INPUTS_DIR, "nodes_samples.csv")
    NODES_TAXA_PATH    = os.path.join(GRAPH_INPUTS_DIR, "nodes_taxa.csv")
    EDGES_SAMPLE_TAXON = os.path.join(GRAPH_INPUTS_DIR, "combined_sample_taxon_edges.csv")
    EDGES_TAXON_TAXON  = os.path.join(GRAPH_INPUTS_DIR, "taxon_taxon_spieceasi_edges.csv")
    TRIPLETS_PATH      = os.path.join(GRAPH_INPUTS_DIR, "coalescence_triplets.csv")
    OUTPUT_DIR         = "gnn_outputs"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── 1. HYPERPARAMETERS ────────────────────────────────────────────────────────
HIDDEN_DIM    = 128      # increased from 64
N_LAYERS      = 3        # increased from 2
DROPOUT       = 0.4
LR            = 5e-4     # slightly lower for stability
LR_PATIENCE   = 5        # reduce LR after 5 epochs no improvement
LR_FACTOR     = 0.5      # multiply LR by this on plateau
MIN_LR        = 1e-5
EPOCHS        = 150
BATCH_SIZE    = 32
PATIENCE      = 20       # increased from 15
RANDOM_SEED   = 42
GRAD_CLIP     = 1.0      # gradient clipping norm

# ── KEY IMPROVEMENT: Limit prediction to top-K most prevalent taxa ────────────
# This reduces the p/n ratio from ~97 to something tractable
# Set to None to use all taxa (original behaviour)
TOP_K_TAXA = 500         # predict only the 500 most prevalent taxa

np.random.seed(RANDOM_SEED)
tf.random.set_seed(RANDOM_SEED)

# ── 2. LOAD DATA ──────────────────────────────────────────────────────────────
print("=" * 60)
print("STEP 1: Loading data")
print("=" * 60)

nodes_samples = pd.read_csv(NODES_SAMPLES_PATH)
nodes_taxa    = pd.read_csv(NODES_TAXA_PATH)
edges_st      = pd.read_csv(EDGES_SAMPLE_TAXON)
triplets      = pd.read_csv(TRIPLETS_PATH)

try:
    edges_tt     = pd.read_csv(EDGES_TAXON_TAXON)
    has_tt_edges = len(edges_tt) > 0
    print(f"  Taxon-taxon edges:   {len(edges_tt):,}")
except Exception:
    edges_tt     = pd.DataFrame()
    has_tt_edges = False
    print("  Taxon-taxon edges:   not available")

print(f"  Sample nodes:        {len(nodes_samples):,}")
print(f"  Taxon nodes:         {len(nodes_taxa):,}")
print(f"  Sample-taxon edges:  {len(edges_st):,}")
print(f"  Triplets:            {len(triplets):,}")

# ── 3. TOP-K TAXON FILTERING ──────────────────────────────────────────────────
# Only keep the most prevalent taxa for prediction
# This is the single most important improvement for p >> n datasets

if TOP_K_TAXA is not None:
    print(f"\nSTEP 1b: Filtering to top {TOP_K_TAXA} most prevalent taxa")

    # Count how many samples each taxon appears in (across final communities only)
    final_edges = edges_st[edges_st["community_type"] == "final"].copy()
    taxon_prevalence = (
        final_edges.groupby(["kingdom", "taxon_id"])["abundance"]
        .apply(lambda x: (x > 0).sum())
        .reset_index()
        .rename(columns={"abundance": "prevalence"})
    )

    # Take top K per kingdom to keep bacteria/fungi balanced
    top_taxa = (
        taxon_prevalence
        .sort_values("prevalence", ascending=False)
        .groupby("kingdom")
        .head(TOP_K_TAXA // 2)   # half from each kingdom
    )

    # Filter nodes_taxa and edges_st
    top_taxa_set = set(
        top_taxa["kingdom"].astype(str) + "::" + top_taxa["taxon_id"].astype(str)
    )

    nodes_taxa_orig = nodes_taxa.copy()
    nodes_taxa["_key"] = nodes_taxa["kingdom"].astype(str) + "::" + nodes_taxa["taxon_id"].astype(str)
    nodes_taxa = nodes_taxa[nodes_taxa["_key"].isin(top_taxa_set)].copy()
    nodes_taxa = nodes_taxa.drop(columns=["_key"]).reset_index(drop=True)

    edges_st["_key"] = edges_st["kingdom"].astype(str) + "::" + edges_st["taxon_id"].astype(str)
    edges_st = edges_st[edges_st["_key"].isin(top_taxa_set)].copy()
    edges_st = edges_st.drop(columns=["_key"]).reset_index(drop=True)

    if has_tt_edges:
        edges_tt["_key_a"] = edges_tt["kingdom"].astype(str) + "::" + edges_tt["taxon_a"].astype(str)
        edges_tt["_key_b"] = edges_tt["kingdom"].astype(str) + "::" + edges_tt["taxon_b"].astype(str)
        edges_tt = edges_tt[
            edges_tt["_key_a"].isin(top_taxa_set) &
            edges_tt["_key_b"].isin(top_taxa_set)
        ].copy()
        edges_tt = edges_tt.drop(columns=["_key_a", "_key_b"]).reset_index(drop=True)

    print(f"  Taxa after filtering: {len(nodes_taxa):,} (from {len(nodes_taxa_orig):,})")
    print(f"  Edges after filtering: {len(edges_st):,}")

# ── 4. BUILD COMPOSITE NODE IDs ───────────────────────────────────────────────
print("\nSTEP 2: Building node indices")

nodes_samples = nodes_samples.copy()
nodes_samples["node_id"] = (
    nodes_samples["kingdom"].astype(str) + "::" +
    nodes_samples["donor_id"].astype(str) + "::" +
    nodes_samples["community_type"].astype(str) + "::" +
    nodes_samples["sample_id"].astype(str)
)

nodes_taxa = nodes_taxa.copy()
nodes_taxa["node_id"] = (
    nodes_taxa["kingdom"].astype(str) + "::" +
    nodes_taxa["taxon_id"].astype(str)
)

all_sample_ids   = nodes_samples["node_id"].tolist()
all_taxon_ids    = nodes_taxa["node_id"].tolist()
sample_id_to_idx = {nid: i for i, nid in enumerate(all_sample_ids)}
taxon_id_to_idx  = {nid: (i + len(all_sample_ids)) for i, nid in enumerate(all_taxon_ids)}
taxon_local_idx  = {nid: i for i, nid in enumerate(all_taxon_ids)}

n_sample_nodes = len(all_sample_ids)
n_taxon_nodes  = len(all_taxon_ids)
n_total_nodes  = n_sample_nodes + n_taxon_nodes

print(f"  Sample nodes:  {n_sample_nodes:,}")
print(f"  Taxon nodes:   {n_taxon_nodes:,}")
print(f"  Total nodes:   {n_total_nodes:,}")

# ── 5. NODE FEATURES ──────────────────────────────────────────────────────────
print("\nSTEP 3: Encoding node features")

le_kingdom   = LabelEncoder()
le_donor     = LabelEncoder()
le_comm_type = LabelEncoder()

nodes_samples["kingdom_enc"]   = le_kingdom.fit_transform(nodes_samples["kingdom"])
nodes_samples["donor_enc"]     = le_donor.fit_transform(nodes_samples["donor_id"])
nodes_samples["comm_type_enc"] = le_comm_type.fit_transform(nodes_samples["community_type"])
sample_features = nodes_samples[["kingdom_enc", "donor_enc", "comm_type_enc"]].values.astype(np.float32)

taxon_feat_cols = []
for col in ["kingdom", "phylum", "class", "order", "family"]:
    if col in nodes_taxa.columns:
        le = LabelEncoder()
        enc_col = f"{col}_enc"
        nodes_taxa[enc_col] = le.fit_transform(
            nodes_taxa[col].fillna("Unknown").astype(str)
        )
        taxon_feat_cols.append(enc_col)

taxon_features = nodes_taxa[taxon_feat_cols].values.astype(np.float32)

max_dim = max(sample_features.shape[1], taxon_features.shape[1])

def pad_to(arr, target):
    if arr.shape[1] < target:
        pad = np.zeros((arr.shape[0], target - arr.shape[1]), dtype=np.float32)
        return np.hstack([arr, pad])
    return arr

sample_features = pad_to(sample_features, max_dim)
taxon_features  = pad_to(taxon_features,  max_dim)
node_features   = np.vstack([sample_features, taxon_features]).astype(np.float32)
scaler          = StandardScaler()
node_features   = scaler.fit_transform(node_features).astype(np.float32)

print(f"  Node feature matrix: {node_features.shape}")

# ── 6. BUILD EDGE LISTS ───────────────────────────────────────────────────────
print("\nSTEP 4: Building edge lists")

edges_st = edges_st.copy()
edges_st["sample_node_id"] = (
    edges_st["kingdom"].astype(str) + "::" +
    edges_st["donor_id"].astype(str) + "::" +
    edges_st["community_type"].astype(str) + "::" +
    edges_st["sample_id"].astype(str)
)
edges_st["taxon_node_id"] = (
    edges_st["kingdom"].astype(str) + "::" +
    edges_st["taxon_id"].astype(str)
)

edges_st["src_idx"] = edges_st["sample_node_id"].map(sample_id_to_idx)
edges_st["dst_idx"] = edges_st["taxon_node_id"].map(taxon_id_to_idx)
edges_st = edges_st.dropna(subset=["src_idx", "dst_idx"]).copy()
edges_st["src_idx"] = edges_st["src_idx"].astype(int)
edges_st["dst_idx"] = edges_st["dst_idx"].astype(int)

edges_st["weight"] = edges_st.groupby("sample_node_id")["abundance"].transform(
    lambda x: x / (x.sum() + 1e-8)
).astype(np.float32)

src_st = np.concatenate([edges_st["src_idx"].values, edges_st["dst_idx"].values])
dst_st = np.concatenate([edges_st["dst_idx"].values, edges_st["src_idx"].values])
wgt_st = np.concatenate([edges_st["weight"].values,  edges_st["weight"].values])

print(f"  Sample-taxon edges (bidirectional): {len(src_st):,}")

if has_tt_edges:
    edges_tt = edges_tt.copy()
    edges_tt["src_node_id"] = edges_tt["kingdom"].astype(str) + "::" + edges_tt["taxon_a"].astype(str)
    edges_tt["dst_node_id"] = edges_tt["kingdom"].astype(str) + "::" + edges_tt["taxon_b"].astype(str)
    edges_tt["src_idx"] = edges_tt["src_node_id"].map(taxon_id_to_idx)
    edges_tt["dst_idx"] = edges_tt["dst_node_id"].map(taxon_id_to_idx)
    edges_tt = edges_tt.dropna(subset=["src_idx", "dst_idx"]).copy()
    edges_tt["src_idx"] = edges_tt["src_idx"].astype(int)
    edges_tt["dst_idx"] = edges_tt["dst_idx"].astype(int)

    src_tt = np.concatenate([edges_tt["src_idx"].values, edges_tt["dst_idx"].values])
    dst_tt = np.concatenate([edges_tt["dst_idx"].values, edges_tt["src_idx"].values])
    wgt_tt = np.ones(len(src_tt), dtype=np.float32)

    all_src = np.concatenate([src_st, src_tt])
    all_dst = np.concatenate([dst_st, dst_tt])
    all_wgt = np.concatenate([wgt_st, wgt_tt])
    print(f"  Taxon-taxon edges (bidirectional):  {len(src_tt):,}")
else:
    all_src, all_dst, all_wgt = src_st, dst_st, wgt_st

print(f"  Total edges:                        {len(all_src):,}")

# ── 7. BUILD TRAINING EXAMPLES ────────────────────────────────────────────────
print("\nSTEP 5: Building training examples")

edges_st_indexed = edges_st.set_index("sample_node_id")

def get_abundance_vector(sample_node_id):
    vec = np.zeros(n_taxon_nodes, dtype=np.float32)
    if sample_node_id not in edges_st_indexed.index:
        return vec
    rows = edges_st_indexed.loc[[sample_node_id]]
    for _, row in rows.iterrows():
        t_id = row["taxon_node_id"]
        if t_id in taxon_local_idx:
            vec[taxon_local_idx[t_id]] = row["weight"]
    return vec

X_resident_idx = []
X_donor_idx    = []
y_final_vec    = []
y_resident_vec = []   # for baseline comparison
triplet_meta   = []
skipped        = 0

for _, row in triplets.iterrows():
    kingdom  = str(row["kingdom"])
    donor_id = str(row["donor_id"])
    res_sid  = str(row["resident_sample_id"])
    fin_sid  = str(row["final_sample_id"])

    res_node_id   = f"{kingdom}::{donor_id}::resident::{res_sid}"
    final_node_id = f"{kingdom}::{donor_id}::final::{fin_sid}"

    if res_node_id not in sample_id_to_idx or final_node_id not in sample_id_to_idx:
        skipped += 1
        continue

    donor_candidates = [k for k in sample_id_to_idx
                        if k.startswith(f"{kingdom}::{donor_id}::donor::")]
    if not donor_candidates:
        skipped += 1
        continue

    res_global_idx   = sample_id_to_idx[res_node_id]
    donor_global_idx = sample_id_to_idx[donor_candidates[0]]

    y_vec  = get_abundance_vector(final_node_id)
    y_res  = get_abundance_vector(res_node_id)

    X_resident_idx.append(res_global_idx)
    X_donor_idx.append(donor_global_idx)
    y_final_vec.append(y_vec)
    y_resident_vec.append(y_res)
    triplet_meta.append({
        "kingdom":  kingdom,
        "donor_id": donor_id,
        "resident_sample_id": res_sid,
        "final_sample_id":    fin_sid
    })

X_resident_idx = np.array(X_resident_idx, dtype=np.int32)
X_donor_idx    = np.array(X_donor_idx,    dtype=np.int32)
y_final_vec    = np.array(y_final_vec,    dtype=np.float32)
y_resident_vec = np.array(y_resident_vec, dtype=np.float32)
triplet_meta_df = pd.DataFrame(triplet_meta)

print(f"  Valid triplets: {len(X_resident_idx):,}  |  Skipped: {skipped:,}")
print(f"  Target shape:   {y_final_vec.shape}")

# ── 8. TRAIN / VAL / TEST SPLIT ───────────────────────────────────────────────
print("\nSTEP 6: Splitting data")

n   = len(X_resident_idx)
idx = np.arange(n)
np.random.shuffle(idx)

n_train = int(0.70 * n)
n_val   = int(0.15 * n)

train_idx = idx[:n_train]
val_idx   = idx[n_train : n_train + n_val]
test_idx  = idx[n_train + n_val :]

print(f"  Train: {len(train_idx)}  |  Val: {len(val_idx)}  |  Test: {len(test_idx)}")

# ── 9. NAIVE BASELINE ─────────────────────────────────────────────────────────
# "Predict that the final community = resident community"
# This is the ecologically meaningful null model for coalescence

def bray_curtis_dissimilarity(y_true, y_pred):
    num = tf.reduce_sum(tf.abs(y_true - y_pred), axis=-1)
    den = tf.reduce_sum(y_true + y_pred, axis=-1) + 1e-8
    return tf.reduce_mean(num / den)

baseline_bc = float(bray_curtis_dissimilarity(
    tf.constant(y_final_vec[test_idx]),
    tf.constant(y_resident_vec[test_idx])
))
print(f"\n  Naive baseline (resident = final) Bray-Curtis: {baseline_bc:.4f}")
print(f"  (GNN must beat this to show it is learning coalescence dynamics)")

# ── 10. MODEL DEFINITION ──────────────────────────────────────────────────────
print("\nSTEP 7: Building model")

class GraphConvLayer(layers.Layer):
    """
    Weighted mean-aggregation graph convolution with residual connection.
    h_v' = FFN( h_v || mean_{u in N(v)}[ h_u * w_{uv} ] )
    """
    def __init__(self, units, dropout_rate=0.3, **kwargs):
        super().__init__(**kwargs)
        self.units        = units
        self.dropout_rate = dropout_rate
        self.ffn = keras.Sequential([
            layers.Dense(units, activation="relu"),
            layers.BatchNormalization(),
            layers.Dropout(dropout_rate),
            layers.Dense(units, activation="relu"),
        ])

    def call(self, node_feats, edge_src, edge_dst, edge_wgt, training=False):
        n_nodes        = tf.shape(node_feats)[0]
        neighbour_feat = tf.gather(node_feats, edge_dst)
        weighted       = neighbour_feat * tf.expand_dims(edge_wgt, -1)
        aggregated     = tf.math.unsorted_segment_mean(
            weighted, edge_src, num_segments=n_nodes
        )
        combined = tf.concat([node_feats, aggregated], axis=-1)
        return self.ffn(combined, training=training)

    # ── FIX: implement get_config() to silence serialization warning ──────────
    def get_config(self):
        config = super().get_config()
        config.update({
            "units":        self.units,
            "dropout_rate": self.dropout_rate
        })
        return config


class CoalescenceGNN(keras.Model):
    """
    Heterogeneous GNN for microbial coalescence prediction.

    Input:  (resident_node_idx, donor_node_idx)  — integer indices
    Output: softmax abundance distribution over top-K taxa
    """
    def __init__(self, n_nodes, n_input_features, n_taxa,
                 hidden_dim=128, n_layers=3, dropout=0.4, **kwargs):
        super().__init__(**kwargs)

        # Store architecture params for get_config()
        self.n_nodes          = n_nodes
        self.n_input_features = n_input_features
        self.n_taxa           = n_taxa
        self.hidden_dim       = hidden_dim
        self.n_layers_        = n_layers
        self.dropout_rate     = dropout

        # Input projection
        self.input_proj = layers.Dense(hidden_dim, activation="relu",
                                       name="input_proj")

        # Graph conv layers
        self.conv_layers = [
            GraphConvLayer(hidden_dim, dropout_rate=dropout, name=f"conv_{i}")
            for i in range(n_layers)
        ]
        self.skip_projs = [
            layers.Dense(hidden_dim, use_bias=False, name=f"skip_{i}")
            for i in range(n_layers)
        ]

        # Decoder
        self.decoder = keras.Sequential([
            layers.Dense(hidden_dim * 2, activation="relu"),
            layers.BatchNormalization(),
            layers.Dropout(dropout),
            layers.Dense(hidden_dim,     activation="relu"),
            layers.Dropout(dropout),
            layers.Dense(n_taxa,         activation="softmax"),
        ], name="decoder")

    def call(self, inputs, node_features, edge_src, edge_dst, edge_wgt,
             training=False):
        """
        inputs       : (resident_idx, donor_idx) — [B] integer tensors
        node_features: [N, F] float tensor (passed explicitly — avoids
                        non-serializable __init__ args warning)
        """
        resident_idx, donor_idx = inputs

        h = self.input_proj(node_features, training=training)

        for conv, skip in zip(self.conv_layers, self.skip_projs):
            h_new = conv(h, edge_src, edge_dst, edge_wgt, training=training)
            h     = h_new + skip(h)

        res_emb   = tf.gather(h, resident_idx)
        donor_emb = tf.gather(h, donor_idx)
        combined  = tf.concat([res_emb, donor_emb], axis=-1)
        return self.decoder(combined, training=training)

    # ── FIX: implement get_config() ───────────────────────────────────────────
    def get_config(self):
        config = super().get_config()
        config.update({
            "n_nodes":          self.n_nodes,
            "n_input_features": self.n_input_features,
            "n_taxa":           self.n_taxa,
            "hidden_dim":       self.hidden_dim,
            "n_layers":         self.n_layers_,
            "dropout":          self.dropout_rate,
        })
        return config


# Convert graph data to TF constants (outside model — avoids serialization issue)
tf_node_features = tf.constant(node_features, dtype=tf.float32)
tf_edge_src      = tf.constant(all_src,        dtype=tf.int32)
tf_edge_dst      = tf.constant(all_dst,        dtype=tf.int32)
tf_edge_wgt      = tf.constant(all_wgt,        dtype=tf.float32)

model = CoalescenceGNN(
    n_nodes          = n_total_nodes,
    n_input_features = node_features.shape[1],
    n_taxa           = n_taxon_nodes,
    hidden_dim       = HIDDEN_DIM,
    n_layers         = N_LAYERS,
    dropout          = DROPOUT,
    name             = "CoalescenceGNN_v2"
)

# ── 11. LOSS AND METRICS ──────────────────────────────────────────────────────

def kl_divergence_loss(y_true, y_pred):
    y_true = tf.clip_by_value(y_true, 1e-8, 1.0)
    y_pred = tf.clip_by_value(y_pred, 1e-8, 1.0)
    return tf.reduce_mean(
        tf.reduce_sum(y_true * tf.math.log(y_true / y_pred), axis=-1)
    )

# ── 12. TRAINING LOOP ─────────────────────────────────────────────────────────
print("\nSTEP 8: Training")
print("=" * 60)
print(f"  p/n ratio: {n_taxon_nodes}/{len(train_idx)} = {n_taxon_nodes/len(train_idx):.1f}")
print(f"  Naive baseline BC: {baseline_bc:.4f}  (target to beat)")
print("=" * 60)

current_lr       = LR
optimizer        = keras.optimizers.Adam(learning_rate=current_lr,
                                         clipnorm=GRAD_CLIP)
best_val_loss    = np.inf
patience_counter = 0
lr_patience_ctr  = 0
weights_path     = os.path.join(OUTPUT_DIR, "best_coalescence_gnn.weights.h5")
history          = {"train_loss": [], "val_loss": [], "val_bc": [], "lr": []}

for epoch in range(EPOCHS):

    # ── Training ──
    shuffled     = np.random.permutation(train_idx)
    epoch_losses = []

    for start in range(0, len(shuffled), BATCH_SIZE):
        batch       = shuffled[start : start + BATCH_SIZE]
        res_batch   = X_resident_idx[batch]
        donor_batch = X_donor_idx[batch]
        y_batch     = y_final_vec[batch]

        with tf.GradientTape() as tape:
            y_pred = model(
                (res_batch, donor_batch),
                node_features = tf_node_features,
                edge_src      = tf_edge_src,
                edge_dst      = tf_edge_dst,
                edge_wgt      = tf_edge_wgt,
                training      = True
            )
            loss = kl_divergence_loss(y_batch, y_pred)

        grads = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(grads, model.trainable_variables))
        epoch_losses.append(float(loss))

    # ── Validation ──
    val_preds = model(
        (X_resident_idx[val_idx], X_donor_idx[val_idx]),
        node_features = tf_node_features,
        edge_src      = tf_edge_src,
        edge_dst      = tf_edge_dst,
        edge_wgt      = tf_edge_wgt,
        training      = False
    )
    val_loss   = float(kl_divergence_loss(y_final_vec[val_idx], val_preds))
    val_bc     = float(bray_curtis_dissimilarity(y_final_vec[val_idx], val_preds))
    train_loss = float(np.mean(epoch_losses))

    history["train_loss"].append(train_loss)
    history["val_loss"].append(val_loss)
    history["val_bc"].append(val_bc)
    history["lr"].append(current_lr)

    is_best = val_loss < best_val_loss
    tag     = " ← best" if is_best else ""
    print(f"Epoch {epoch+1:3d}/{EPOCHS}  "
          f"train_KL={train_loss:.4f}  "
          f"val_KL={val_loss:.4f}  "
          f"val_BC={val_bc:.4f}  "
          f"lr={current_lr:.2e}{tag}")

    # ── Early stopping ──
    if is_best:
        best_val_loss    = val_loss
        patience_counter = 0
        lr_patience_ctr  = 0
        model.save_weights(weights_path)
    else:
        patience_counter += 1
        lr_patience_ctr  += 1

        # Learning rate reduction on plateau
        if lr_patience_ctr >= LR_PATIENCE and current_lr > MIN_LR:
            current_lr = max(current_lr * LR_FACTOR, MIN_LR)
            optimizer.learning_rate.assign(current_lr)
            lr_patience_ctr = 0
            print(f"  → Learning rate reduced to {current_lr:.2e}")

        if patience_counter >= PATIENCE:
            print(f"\nEarly stopping at epoch {epoch + 1}.")
            break

# ── 13. TEST EVALUATION ───────────────────────────────────────────────────────
print("\nSTEP 9: Test evaluation")
print("=" * 60)

model.load_weights(weights_path)

test_preds = model(
    (X_resident_idx[test_idx], X_donor_idx[test_idx]),
    node_features = tf_node_features,
    edge_src      = tf_edge_src,
    edge_dst      = tf_edge_dst,
    edge_wgt      = tf_edge_wgt,
    training      = False
)
test_kl = float(kl_divergence_loss(y_final_vec[test_idx], test_preds))
test_bc = float(bray_curtis_dissimilarity(y_final_vec[test_idx], test_preds))

print(f"  Naive baseline BC:              {baseline_bc:.4f}")
print(f"  GNN Test Bray-Curtis:           {test_bc:.4f}")
print(f"  GNN Test KL Divergence:         {test_kl:.4f}")
improvement = baseline_bc - test_bc
print(f"  Improvement over baseline:      {improvement:+.4f} "
      f"({'better' if improvement > 0 else 'worse'} than naive)")

# ── Per-kingdom evaluation ────────────────────────────────────────────────────
print("\n  Per-kingdom breakdown:")
test_meta = triplet_meta_df.iloc[test_idx].reset_index(drop=True)
test_preds_np = test_preds.numpy()

for kingdom in ["Bacteria", "Fungi"]:
    k_mask = test_meta["kingdom"] == kingdom
    if k_mask.sum() == 0:
        continue
    k_idx = np.where(k_mask.values)[0]
    k_bc  = float(bray_curtis_dissimilarity(
        tf.constant(y_final_vec[test_idx][k_idx]),
        tf.constant(test_preds_np[k_idx])
    ))
    k_base = float(bray_curtis_dissimilarity(
        tf.constant(y_final_vec[test_idx][k_idx]),
        tf.constant(y_resident_vec[test_idx][k_idx])
    ))
    print(f"    {kingdom:10s}  GNN BC={k_bc:.4f}  baseline BC={k_base:.4f}  "
          f"n={k_mask.sum()}")

# ── 14. SAVE OUTPUTS ──────────────────────────────────────────────────────────
print("\nSTEP 10: Saving outputs")

# Predictions — FIX: reset_index() prevents duplicate index warning
pred_df = pd.DataFrame(
    test_preds_np,
    columns = all_taxon_ids
)
# Add metadata columns at the front
pred_df.insert(0, "kingdom",            test_meta["kingdom"].values)
pred_df.insert(1, "donor_id",           test_meta["donor_id"].values)
pred_df.insert(2, "resident_sample_id", test_meta["resident_sample_id"].values)
pred_df.insert(3, "final_sample_id",    test_meta["final_sample_id"].values)
pred_df = pred_df.reset_index(drop=True)   # ← fixes duplicate index warning

pred_path = os.path.join(OUTPUT_DIR, "gnn_predictions_test.csv")
pred_df.to_csv(pred_path, index=False)
print(f"  Predictions:      {pred_path}")

# Training history
hist_df = pd.DataFrame(history)
hist_df.insert(0, "epoch", range(1, len(hist_df) + 1))
hist_path = os.path.join(OUTPUT_DIR, "gnn_training_history.csv")
hist_df.to_csv(hist_path, index=False)
print(f"  Training history: {hist_path}")
print(f"  Model weights:    {weights_path}")

# Summary metrics
summary = {
    "n_taxa_predicted":   n_taxon_nodes,
    "n_train":            len(train_idx),
    "n_val":              len(val_idx),
    "n_test":             len(test_idx),
    "baseline_bc":        round(baseline_bc, 4),
    "test_kl":            round(test_kl, 4),
    "test_bc":            round(test_bc, 4),
    "improvement_vs_baseline": round(improvement, 4),
    "epochs_trained":     len(history["train_loss"]),
    "best_val_kl":        round(best_val_loss, 4),
    "top_k_taxa_filter":  TOP_K_TAXA if TOP_K_TAXA else "none",
    "hidden_dim":         HIDDEN_DIM,
    "n_layers":           N_LAYERS,
}
summary_df = pd.DataFrame([summary])
summary_path = os.path.join(OUTPUT_DIR, "gnn_run_summary.csv")
summary_df.to_csv(summary_path, index=False)
print(f"  Run summary:      {summary_path}")

print("\nPipeline complete.")
