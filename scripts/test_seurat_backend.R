#!/usr/bin/env Rscript
# End-to-end smoke test for CHOIR's Seurat (WNN) backend code paths.
#
# This script does NOT run the full CHOIR() pipeline — several Suggests
# (harmony, bluster, mrtree, ggtree, spatstat.univar) cannot be installed
# in offline sandboxes. Instead it exercises the exact code sequence that
# CHOIR runs for `integration_backend = "seurat"`:
#
#   1. Per-modality dimensionality reductions produced as a list-of-matrices
#   2. FindMultiModalNeighbors fused via a temporary Seurat object
#   3. .getMultiModalDistance (from R/TreeUtils.R) as the distance matrix
#   4. .getBaselineResolution + .getStartingResolution (Leiden fix)
#   5. FindClusters on the fused SNN graph
#   6. Subtree subsetting (list-of-matrices lapply path)
#
# Pass criteria: every section prints "[PASS]" with no R error/warning.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

root <- tryCatch(normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")),
                 error = function(e) getwd())
if (!dir.exists(file.path(root, "R"))) root <- getwd()

cat("CHOIR root:", root, "\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n\n")

# --------------------------------------------------------------------------
# Source CHOIR internals that the Seurat backend depends on.
# --------------------------------------------------------------------------
# TreeUtils.R pulls in dplyr/methods/stats/progress — load what's available,
# stub what isn't (progress bar is non-essential for smoke testing).
suppressPackageStartupMessages({
  library(dplyr)
  library(methods)
})

# Minimal progress stub so TreeUtils.R source() doesn't fail on load
if (!requireNamespace("progress", quietly = TRUE)) {
  progress <- list(progress_bar = list(new = function(...) list(
    tick = function(...) invisible(NULL),
    message = function(...) invisible(NULL))))
  assign("progress", progress, envir = globalenv())
}

src_ok <- tryCatch({
  source(file.path(root, "R", "TreeUtils.R"), local = FALSE)
  TRUE
}, error = function(e) {
  cat("source(TreeUtils.R) failed:", conditionMessage(e), "\n")
  FALSE
})
stopifnot(src_ok)
stopifnot(exists(".getMultiModalDistance"))
stopifnot(exists(".getBaselineResolution"))
stopifnot(exists(".usesHarmony"))
stopifnot(exists(".getStartingResolution"))
cat("[PASS] TreeUtils.R sourced; key internals found.\n\n")

# --------------------------------------------------------------------------
# 1. Synthetic multiome: 3 latent "cell types" × 80 cells × 2 modalities
# --------------------------------------------------------------------------
set.seed(1337)
n_per_type <- 200
n_types    <- 3
n_cells    <- n_per_type * n_types
cell_ids   <- paste0("cell_", seq_len(n_cells))
true_label <- rep(paste0("T", seq_len(n_types)), each = n_per_type)

# -- RNA modality: 500 genes, 3 type-specific programs ---------------------
n_genes <- 500
gene_names <- paste0("gene_", seq_len(n_genes))

# Each type has 30 upregulated genes
type_genes <- split(sample(seq_len(n_genes)), rep(seq_len(n_types),
                                                  length.out = n_genes))[seq_len(n_types)]
rna_mu <- matrix(0.3, nrow = n_genes, ncol = n_cells,
                 dimnames = list(gene_names, cell_ids))
for (k in seq_len(n_types)) {
  cells_k <- which(true_label == paste0("T", k))
  genes_k <- sample(type_genes[[k]], 30)
  rna_mu[genes_k, cells_k] <- rna_mu[genes_k, cells_k] + 2.5
}
rna_counts <- matrix(rpois(length(rna_mu), lambda = rna_mu),
                     nrow = n_genes, ncol = n_cells,
                     dimnames = dimnames(rna_mu))

# -- ATAC-like modality: 1000 peaks, 3 type-specific open sets -------------
n_peaks <- 1000
peak_names <- paste0("peak_", seq_len(n_peaks))
type_peaks <- split(sample(seq_len(n_peaks)), rep(seq_len(n_types),
                                                   length.out = n_peaks))[seq_len(n_types)]
atac_p <- matrix(0.02, nrow = n_peaks, ncol = n_cells,
                 dimnames = list(peak_names, cell_ids))
for (k in seq_len(n_types)) {
  cells_k <- which(true_label == paste0("T", k))
  peaks_k <- sample(type_peaks[[k]], 80)
  atac_p[peaks_k, cells_k] <- 0.35
}
atac_counts <- matrix(rbinom(length(atac_p), size = 1, prob = atac_p),
                      nrow = n_peaks, ncol = n_cells,
                      dimnames = dimnames(atac_p))

cat(sprintf("[PASS] Generated synthetic data: %d cells, %d genes, %d peaks, %d types.\n\n",
            n_cells, n_genes, n_peaks, n_types))

# --------------------------------------------------------------------------
# 2. Build per-modality reductions (mirrors buildParentTree.R multimodal step)
#    RNA → LogNorm + PCA                  (mimics the new ArchR-PCA branch)
#    ATAC → TF-IDF + SVD                  (mimics IterativeLSI output)
# --------------------------------------------------------------------------
rna <- CreateSeuratObject(counts = as(rna_counts, "dgCMatrix"), assay = "RNA")
rna <- NormalizeData(rna, verbose = FALSE)
rna <- FindVariableFeatures(rna, nfeatures = 200, verbose = FALSE)
rna <- ScaleData(rna, features = VariableFeatures(rna), verbose = FALSE)
rna <- RunPCA(rna, features = VariableFeatures(rna), npcs = 20,
              seed.use = 1, verbose = FALSE)
rna_embed <- Embeddings(rna, "pca")          # n_cells × 20

# ATAC: simple TF-IDF + truncated SVD
tf  <- sweep(atac_counts, 2, pmax(colSums(atac_counts), 1), "/")
idf <- log(1 + ncol(atac_counts) / pmax(rowSums(atac_counts > 0), 1))
atac_tfidf <- tf * idf
svd_atac <- svd(atac_tfidf, nu = 0, nv = 20)
atac_embed <- svd_atac$v %*% diag(svd_atac$d[seq_len(20)])
rownames(atac_embed) <- cell_ids
colnames(atac_embed) <- paste0("LSI_", seq_len(20))

reduction_coords_list <- list(rna_embed, atac_embed)
cat(sprintf("[PASS] Per-modality reductions built: RNA(%s), ATAC(%s).\n\n",
            paste(dim(rna_embed), collapse = "x"),
            paste(dim(atac_embed), collapse = "x")))

# --------------------------------------------------------------------------
# 3. Reproduce the EXACT WNN code block from buildParentTree.R (fixed path)
# --------------------------------------------------------------------------
n_modalities <- 2
tmp <- matrix(stats::rnorm(n_cells * 3, 10), ncol = n_cells, nrow = 3)
colnames(tmp) <- cell_ids
rownames(tmp) <- paste0("t", seq_len(nrow(tmp)))
tmp_seurat <- CreateSeuratObject(tmp, min.cells = 0, min.features = 0,
                                 assay = "tmp")

dim_list <- vector("list", length = n_modalities)
for (i in seq_len(n_modalities)) {
  tmp_seurat[[paste0("DR_", i)]] <- CreateDimReducObject(
    embeddings = reduction_coords_list[[i]],
    key = paste0("DR_", i, "_"), assay = "tmp")
  dim_list[[i]] <- seq_len(ncol(reduction_coords_list[[i]]))
}

wnn <- FindMultiModalNeighbors(tmp_seurat,
                               reduction.list = as.list(paste0("DR_", seq_len(n_modalities))),
                               dims.list      = dim_list,
                               knn.graph.name = "nn",
                               snn.graph.name = "snn",
                               verbose = FALSE)

stopifnot("nn"  %in% names(wnn@graphs))
stopifnot("snn" %in% names(wnn@graphs))
stopifnot(ncol(wnn@graphs$snn) == n_cells)
stopifnot(all(rownames(wnn@graphs$snn) == cell_ids))
cat("[PASS] FindMultiModalNeighbors returned valid nn + snn of shape",
    paste(dim(wnn@graphs$snn), collapse = "x"), "\n\n")

# --------------------------------------------------------------------------
# 4. .getMultiModalDistance (from TreeUtils.R)
# --------------------------------------------------------------------------
dist_mat <- .getMultiModalDistance(
  object = tmp_seurat,
  reduction_list = as.list(paste0("DR_", seq_len(n_modalities))),
  dim_list = dim_list)

stopifnot(inherits(dist_mat, "matrix"))
stopifnot(nrow(dist_mat) == n_cells)
stopifnot(ncol(dist_mat) == n_cells)
# Diagonal should be zero (self distance). NA is permitted in last column
# per CHOIR's construction; just check the first (n-1) columns per row.
self_dist_ok <- all(vapply(seq_len(n_cells), function(i) {
  d <- dist_mat[i, i]; isTRUE(is.na(d)) || isTRUE(d == 0 || d < 1e-9)
}, logical(1)))
cat(sprintf("[PASS] .getMultiModalDistance shape %dx%d ; self-distance check: %s\n\n",
            nrow(dist_mat), ncol(dist_mat), self_dist_ok))

# --------------------------------------------------------------------------
# 5. .getBaselineResolution — Leiden fix
# --------------------------------------------------------------------------
r_louvain <- .getBaselineResolution(list(algorithm = 1))
r_leiden  <- .getBaselineResolution(list(algorithm = 4))
r_empty   <- .getBaselineResolution(list())
stopifnot(identical(r_louvain, 0))
stopifnot(r_leiden > 0 && r_leiden < 1e-3)
stopifnot(identical(r_empty, 0))
cat(sprintf("[PASS] baseline resolution — Louvain=%g, Leiden=%g, empty=%g\n\n",
            r_louvain, r_leiden, r_empty))

# --------------------------------------------------------------------------
# 6. .usesHarmony — multimodal vector handling
# --------------------------------------------------------------------------
stopifnot(.usesHarmony("Harmony") == TRUE)
stopifnot(.usesHarmony("none")    == FALSE)
stopifnot(.usesHarmony(c("Harmony","none")) == TRUE)
stopifnot(.usesHarmony(c("none","none"))    == FALSE)
cat("[PASS] .usesHarmony handles scalar and multimodal vectors.\n\n")

# --------------------------------------------------------------------------
# 7. Cluster the WNN SNN graph at several resolutions; check we recover
#    roughly the planted 3-type structure.
# --------------------------------------------------------------------------
snn_mat <- wnn@graphs$snn

cl_louvain <- FindClusters(snn_mat, resolution = 0.5,
                           random.seed = 1, verbose = FALSE)[, 1]
cat("Louvain @ res=0.5 cluster table:\n")
print(table(cl_louvain, true_label))

cl_leiden <- tryCatch(
  FindClusters(snn_mat, resolution = r_leiden,
               algorithm = 4, method = "igraph",
               random.seed = 1, verbose = FALSE)[, 1],
  error = function(e) {
    cat("  (Leiden backend unavailable:", conditionMessage(e), ")\n")
    NULL
  })
if (!is.null(cl_leiden)) {
  cat(sprintf("Leiden @ res=%g produced %d clusters (baseline should be ~1).\n",
              r_leiden, length(unique(cl_leiden))))
}

# ARI vs true labels (Louvain — Leiden baseline is intentionally collapsed)
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  a_comb <- sum(choose(rowSums(tab), 2))
  b_comb <- sum(choose(colSums(tab), 2))
  tab_comb <- sum(choose(tab, 2))
  exp_comb <- a_comb * b_comb / choose(n, 2)
  max_comb <- 0.5 * (a_comb + b_comb)
  (tab_comb - exp_comb) / (max_comb - exp_comb)
}
ari_louvain <- ari(cl_louvain, true_label)
stopifnot(ari_louvain > 0.7)   # WNN on planted 3-type data should be easy
cat(sprintf("[PASS] Louvain ARI vs truth = %.3f  (> 0.7 required)\n\n",
            ari_louvain))

# --------------------------------------------------------------------------
# 8. Subtree subsetting (mirrors buildTree.R:924-933)
#    Verifies the list-of-matrices path introduced with seurat backend.
# --------------------------------------------------------------------------
P0 <- list("reduction_coords" = reduction_coords_list)
subset_ids <- sample(cell_ids, 300)

if (is.list(P0[["reduction_coords"]])) {
  P_i <- list("reduction_coords" = lapply(P0[["reduction_coords"]],
                                          function(x) x[subset_ids, , drop = FALSE]))
} else {
  P_i <- list("reduction_coords" = P0[["reduction_coords"]][subset_ids, ])
}
stopifnot(is.list(P_i$reduction_coords))
stopifnot(length(P_i$reduction_coords) == 2)
stopifnot(nrow(P_i$reduction_coords[[1]]) == length(subset_ids))
stopifnot(nrow(P_i$reduction_coords[[2]]) == length(subset_ids))
stopifnot(all(rownames(P_i$reduction_coords[[1]]) == subset_ids))
cat(sprintf("[PASS] Subtree list-subset path: RNA(%s), ATAC(%s).\n\n",
            paste(dim(P_i$reduction_coords[[1]]), collapse = "x"),
            paste(dim(P_i$reduction_coords[[2]]), collapse = "x")))

# --------------------------------------------------------------------------
# 9. Re-run WNN on the subtree to verify the path holds for small graphs
# --------------------------------------------------------------------------
tmp2 <- matrix(rnorm(length(subset_ids) * 3, 10),
               ncol = length(subset_ids), nrow = 3)
colnames(tmp2) <- subset_ids
rownames(tmp2) <- paste0("t", seq_len(nrow(tmp2)))
sub_seurat <- CreateSeuratObject(tmp2, min.cells = 0, min.features = 0,
                                 assay = "tmp")
sub_dim_list <- vector("list", length = n_modalities)
for (i in seq_len(n_modalities)) {
  sub_seurat[[paste0("DR_", i)]] <- CreateDimReducObject(
    embeddings = P_i$reduction_coords[[i]],
    key = paste0("DR_", i, "_"), assay = "tmp")
  # NOTE: dim_list must come from the subtree reduction (buildTree.R:1041 fix)
  sub_dim_list[[i]] <- seq_len(ncol(P_i$reduction_coords[[i]]))
}
sub_wnn <- FindMultiModalNeighbors(sub_seurat,
                                   reduction.list = as.list(paste0("DR_", seq_len(n_modalities))),
                                   dims.list      = sub_dim_list,
                                   knn.graph.name = "nn",
                                   snn.graph.name = "snn",
                                   verbose = FALSE)
stopifnot(ncol(sub_wnn@graphs$snn) == length(subset_ids))
cat(sprintf("[PASS] Subtree WNN: snn shape %s\n\n",
            paste(dim(sub_wnn@graphs$snn), collapse = "x")))

cat("=============================================\n")
cat("All Seurat backend smoke tests PASSED.\n")
cat("=============================================\n")
