#!/usr/bin/env Rscript
# End-to-end smoke test for integration_backend = "seurat" on an ArchRProject.
#
# We cannot install the full ArchR package stack in this offline sandbox
# (chromVAR / motifmatchr / presto / chromVARmotifs require CRAN/Bioc access),
# and ArchR's getTestProject() downloads from a blocked S3 bucket.
#
# Strategy:
#   1. Define the real ArchRProject S4 class (from ArchR/R/AllClasses.R) using
#      only S4Vectors + GenomicRanges from apt.
#   2. Populate the S4 slots with synthetic multiome counts so
#      `methods::is(obj, "ArchRProject")` is TRUE — the exact dispatch CHOIR
#      keys off.
#   3. Stub the handful of ArchR:: functions that CHOIR actually calls
#      (getMatrixFromProject, addIterativeLSI, addCombinedDims,
#      getCellNames) with in-R implementations that return data in the
#      same shape ArchR does.
#   4. Drive CHOIR's .runDimReduction() with
#        normalization_method = c("none",   "LogNorm"),
#        reduction_method    = c("IterativeLSI", "PCA"),
#        ArchR_matrix        = c("TileMatrix", "GeneExpressionMatrix"),
#        atac                = c(TRUE, FALSE),
#        integration_backend = "seurat"
#      so it produces a list-of-matrices (one reduction per modality).
#   5. Run that list through the EXACT WNN block copied from
#      buildParentTree.R lines 540–572 (post-fix).
#   6. Louvain-cluster the WNN SNN and verify ARI vs planted labels.

suppressPackageStartupMessages({
  library(S4Vectors)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

root <- tryCatch(normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")),
                 error = function(e) getwd())
if (!dir.exists(file.path(root, "R"))) root <- getwd()
cat("CHOIR root:", root, "\n")
cat("R:", R.version.string, "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n\n")

# --------------------------------------------------------------------------
# 1. Define the real ArchRProject S4 class (copied verbatim from
#    ArchR/R/AllClasses.R, with the validity check relaxed — we do not have
#    Arrow files on disk).
# --------------------------------------------------------------------------
suppressMessages(setClass("ArchRProject",
  representation(
    projectMetadata   = "SimpleList",
    projectSummary    = "SimpleList",
    sampleColData     = "DataFrame",
    sampleMetadata    = "SimpleList",
    cellColData       = "DataFrame",
    cellMetadata      = "SimpleList",
    reducedDims       = "SimpleList",
    embeddings        = "SimpleList",
    peakSet           = "ANY",
    peakAnnotation    = "SimpleList",
    geneAnnotation    = "SimpleList",
    genomeAnnotation  = "SimpleList",
    imputeWeights     = "SimpleList"
  )))
suppressMessages(setValidity("ArchRProject", function(object) TRUE))

# --------------------------------------------------------------------------
# 2. Synthesize a tiny but realistic multiome dataset.
#    600 cells across 3 planted types, 2 batches (for Harmony-adjacent tests
#    later), RNA + ATAC modalities driven by the shared latent.
# --------------------------------------------------------------------------
set.seed(0xC401)
n_per_type <- 250
n_types    <- 3
n_cells    <- n_per_type * n_types
cell_ids   <- paste0("cell_", seq_len(n_cells))
true_label <- rep(paste0("T", seq_len(n_types)), each = n_per_type)
batch_lbl  <- rep(c("B1", "B2"), length.out = n_cells)

# RNA ---------------------------------------------------------------------
n_genes    <- 600
gene_names <- paste0("ENSG", sprintf("%05d", seq_len(n_genes)))
type_genes <- split(sample(seq_len(n_genes)),
                    rep(seq_len(n_types), length.out = n_genes))[seq_len(n_types)]
rna_mu <- matrix(0.3, nrow = n_genes, ncol = n_cells,
                 dimnames = list(gene_names, cell_ids))
for (k in seq_len(n_types)) {
  idx <- which(true_label == paste0("T", k))
  g   <- sample(type_genes[[k]], 40)
  rna_mu[g, idx] <- rna_mu[g, idx] + 2.5
}
rna_counts <- matrix(rpois(length(rna_mu), lambda = rna_mu),
                     nrow = n_genes, dimnames = dimnames(rna_mu))

# ATAC --------------------------------------------------------------------
n_tiles    <- 1200
tile_names <- paste0("chr1:", seq_len(n_tiles))
type_tiles <- split(sample(seq_len(n_tiles)),
                    rep(seq_len(n_types), length.out = n_tiles))[seq_len(n_types)]
atac_p <- matrix(0.02, nrow = n_tiles, ncol = n_cells,
                 dimnames = list(tile_names, cell_ids))
for (k in seq_len(n_types)) {
  idx <- which(true_label == paste0("T", k))
  p   <- sample(type_tiles[[k]], 100)
  atac_p[p, idx] <- 0.35
}
atac_counts <- matrix(rbinom(length(atac_p), size = 1, prob = atac_p),
                      nrow = n_tiles, dimnames = dimnames(atac_p))

# Cell metadata (per-cell depth, batch) -----------------------------------
cellColData <- DataFrame(
  Sample   = rep("mock_sample", n_cells),
  Gex_nUMI = colSums(rna_counts),
  nFrags   = colSums(atac_counts),
  Batch    = batch_lbl,
  row.names = cell_ids
)
cat(sprintf("[PASS] Synthesized: %d cells, %d types, %d batches, %d genes, %d tiles.\n\n",
            n_cells, n_types, length(unique(batch_lbl)),
            n_genes, n_tiles))

# --------------------------------------------------------------------------
# 3. Construct the ArchRProject object.
# --------------------------------------------------------------------------
sampleColData <- DataFrame(ArrowFiles = "mock.arrow", row.names = "mock_sample")
archr_obj <- new("ArchRProject",
  projectMetadata  = SimpleList(outputDirectory = tempdir()),
  projectSummary   = SimpleList(),
  sampleColData    = sampleColData,
  sampleMetadata   = SimpleList(),
  cellColData      = cellColData,
  cellMetadata     = SimpleList(),
  reducedDims      = SimpleList(),
  embeddings       = SimpleList(),
  peakSet          = NULL,
  peakAnnotation   = SimpleList(),
  geneAnnotation   = SimpleList(),
  genomeAnnotation = SimpleList(),
  imputeWeights    = SimpleList()
)
stopifnot(methods::is(archr_obj, "ArchRProject"))
cat("[PASS] ArchRProject instance built: class =", class(archr_obj), "\n\n")

# --------------------------------------------------------------------------
# 4. Stub the ArchR:: functions that CHOIR calls. We attach these to a
#    fake ArchR namespace so `ArchR::fn()` dispatches to our stubs.
# --------------------------------------------------------------------------
# Per-modality count matrices live here; getMatrixFromProject looks them up.
.mock_matrices <- list(
  TileMatrix           = atac_counts,
  GeneExpressionMatrix = rna_counts
)

# Build a minimal installable "ArchR" shim package on disk so that
# `ArchR::fn()` works through the normal namespace machinery.
build_fake_archr_pkg <- function(pkg_dir, funs_src_path) {
  dir.create(pkg_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(pkg_dir, "R"), showWarnings = FALSE)
  writeLines(c(
    "Package: ArchR",
    "Type: Package",
    "Title: Mock ArchR for CHOIR tests",
    "Version: 0.0.0",
    "Depends: R (>= 4.0)",
    "Imports: S4Vectors, SummarizedExperiment, matrixStats, methods",
    "Description: Stubbed API for CHOIR smoke testing.",
    "License: MIT"
  ), file.path(pkg_dir, "DESCRIPTION"))
  writeLines(c(
    paste0("export(", paste(names(archr_funs), collapse = ", "), ")"),
    "import(S4Vectors)",
    "import(SummarizedExperiment)",
    "importFrom(matrixStats, rowVars)",
    "importFrom(methods, new)"
  ), file.path(pkg_dir, "NAMESPACE"))
  file.copy(funs_src_path, file.path(pkg_dir, "R", "funs.R"), overwrite = TRUE)
}

# Build the mocked ArchR function set
`%||%` <- function(a, b) if (is.null(a)) b else a

archr_funs <- list(
  getMatrixFromProject = function(ArchRProj, useMatrix, ...) {
    m <- .mock_matrices[[useMatrix]]
    if (is.null(m)) stop("mock getMatrixFromProject: unknown matrix '", useMatrix, "'")
    SummarizedExperiment::SummarizedExperiment(
      assays  = setNames(list(m), useMatrix),
      rowData = DataFrame(name = rownames(m)),
      colData = ArchRProj@cellColData[colnames(m), , drop = FALSE])
  },
  addIterativeLSI = function(ArchRProj, useMatrix, name, varFeatures = 1000,
                             depthCol = "nFrags", force = TRUE, threads = 1,
                             seed = 1, ...) {
    m <- .mock_matrices[[useMatrix]]
    cells_in <- rownames(ArchRProj@cellColData)
    m <- m[, cells_in, drop = FALSE]
    feat_var <- matrixStats::rowVars(as.matrix(m))
    top_idx  <- order(feat_var, decreasing = TRUE)[seq_len(min(varFeatures, length(feat_var)))]
    m_sel <- m[top_idx, , drop = FALSE]
    cs <- pmax(colSums(m_sel), 1)
    rs <- pmax(rowSums(m_sel > 0), 1)
    tf  <- sweep(m_sel, 2, cs, "/")
    idf <- log(1 + ncol(m_sel) / rs)
    mat <- tf * idf
    sv  <- svd(mat, nu = 0, nv = 30)
    mSVD <- sv$v %*% diag(sv$d[seq_len(30)])
    rownames(mSVD) <- colnames(m_sel)
    colnames(mSVD) <- paste0("LSI", seq_len(30))
    depth_vec <- as.numeric(ArchRProj@cellColData[cells_in, depthCol])
    ctd <- abs(cor(mSVD, depth_vec, use = "pairwise.complete.obs"))[, 1]
    ArchRProj@reducedDims[[name]] <- SimpleList(
      matSVD      = mSVD,
      matDR       = mSVD,
      LSIFeatures = DataFrame(name = rownames(m_sel)),
      scaleDims   = TRUE,
      corToDepth  = list(scaled = ctd, unscaled = ctd),
      params      = list(reduction_method = "IterativeLSI"),
      date        = Sys.time())
    ArchRProj
  },
  addCombinedDims = function(ArchRProj, reducedDims, name, ...) {
    mats <- lapply(reducedDims, function(rd) {
      rd_obj <- ArchRProj@reducedDims[[rd]]
      if (is.null(rd_obj)) stop("mock addCombinedDims: reducedDim '", rd, "' not found")
      m <- rd_obj$matSVD
      if (is.null(m)) m <- rd_obj$matDR
      if (is.null(m)) m <- rd_obj[[1]]
      scale(m)
    })
    combined <- do.call(cbind, mats)
    colnames(combined) <- paste0("CD", seq_len(ncol(combined)))
    rownames(combined) <- rownames(mats[[1]])
    ArchRProj@reducedDims[[name]] <- SimpleList(
      matDR = combined, matSVD = combined,
      params = list(reduction_method = "Combined"),
      scaleDims = NA, corToDepth = NA, date = Sys.time())
    ArchRProj
  },
  getCellNames = function(ArchRProj, ...) rownames(ArchRProj@cellColData)
)

# Persist the stub functions to a file (for R CMD INSTALL) and install.
pkg_tmp <- file.path(tempdir(), "ArchR_stub")
funs_file <- file.path(tempdir(), "archr_stub_funs.R")
.mock_matrices_global <- .mock_matrices
# Write stub function definitions to the package R/ dir. We reference
# `.mock_matrices` through a global getter so the installed package can
# reach the data stored in this session.
dump(c("getMatrixFromProject", "addIterativeLSI", "addCombinedDims", "getCellNames"),
     envir = list2env(archr_funs), file = funs_file)
# getMatrixFromProject/addIterativeLSI reference `.mock_matrices` from the
# caller's global env — prefix with `.GlobalEnv$`.
src <- readLines(funs_file)
src <- gsub("\\.mock_matrices\\b", ".GlobalEnv$.mock_matrices", src)
writeLines(src, funs_file)

build_fake_archr_pkg(pkg_tmp, funs_file)
install.packages(pkg_tmp, repos = NULL, type = "source",
                 lib = Sys.getenv("R_LIBS_USER", unset = .libPaths()[1]),
                 quiet = TRUE)
stopifnot(requireNamespace("ArchR", quietly = TRUE))
stopifnot(is.function(ArchR::getMatrixFromProject))
cat("[PASS] ArchR stub package installed &",
    length(archr_funs), "functions exported.\n\n")

# --------------------------------------------------------------------------
# 5. Source CHOIR TreeUtils.R so .runDimReduction and helpers are available.
#    We also need an inline `.requirePackage` shim because CHOIR's internal
#    helper checks for package availability.
# --------------------------------------------------------------------------
.requirePackage <- function(x, source = NULL, installInfo = NULL) {
  invisible(TRUE)
}
assign(".requirePackage", .requirePackage, envir = globalenv())

# `%||%` is used throughout CHOIR
`%||%` <- function(a, b) if (is.null(a)) b else a
assign("%||%", `%||%`, envir = globalenv())

# progress stub (TreeUtils uses it in tree-building; not needed for DR path)
if (!requireNamespace("progress", quietly = TRUE)) {
  progress_bar_new <- function(...) list(
    tick   = function(...) invisible(NULL),
    message = function(...) invisible(NULL))
  progress <- list(progress_bar = list(new = progress_bar_new))
  assign("progress", progress, envir = globalenv())
}

# HelperUtils.R has .matchArg, .getCellIDs, .retrieveData, .storeData that
# .runDimReduction uses. Source it too — it has a .requirePackage that would
# overwrite ours, so we re-stub afterwards.
source(file.path(root, "R", "HelperUtils.R"), local = FALSE)
.requirePackage <- function(x, source = NULL, installInfo = NULL) invisible(TRUE)
assign(".requirePackage", .requirePackage, envir = globalenv())

source(file.path(root, "R", "TreeUtils.R"), local = FALSE)
stopifnot(exists(".runDimReduction"))
stopifnot(exists(".getMultiModalDistance"))
cat("[PASS] CHOIR TreeUtils.R sourced.\n\n")

# --------------------------------------------------------------------------
# 6. Call .runDimReduction() — the function is recursive for multimodal,
#    so n_modalities=2 triggers per-modality DR then combines.
#    With integration_backend="seurat", the ArchR branch should return a
#    LIST of per-modality reduction matrices (not a single combined matrix).
# --------------------------------------------------------------------------
dr_out <- .runDimReduction(
  object                  = archr_obj,
  normalization_method    = c("none",         "LogNorm"),
  reduction_method        = c("IterativeLSI", "PCA"),
  reduction_params        = list(),
  n_var_features          = c(500, 300),
  batch_correction_method = c("none", "none"),
  batch_correction_params = list(),
  batch_labels            = NULL,
  use_assay               = NULL,
  use_slot                = NULL,
  ArchR_matrix            = c("TileMatrix", "GeneExpressionMatrix"),
  ArchR_depthcol          = c("nFrags", "Gex_nUMI"),
  atac                    = c(TRUE, FALSE),
  integration_backend     = "seurat",
  use_cells               = NULL,
  return_full             = TRUE,
  n_cores                 = 1,
  random_seed             = 1,
  verbose                 = FALSE
)

stopifnot(is.list(dr_out))
stopifnot("reduction_coords" %in% names(dr_out))
rc <- dr_out[["reduction_coords"]]
stopifnot(is.list(rc))
stopifnot(length(rc) == 2)
stopifnot(is.matrix(rc[[1]]) || inherits(rc[[1]], "Matrix"))
stopifnot(is.matrix(rc[[2]]) || inherits(rc[[2]], "Matrix"))
stopifnot(nrow(rc[[1]]) == n_cells)
stopifnot(nrow(rc[[2]]) == n_cells)
cat(sprintf("[PASS] .runDimReduction on ArchRProject + integration_backend='seurat'\n"))
cat(sprintf("       Modality 1 (ATAC/IterativeLSI stub): %s\n",
            paste(dim(rc[[1]]), collapse = "x")))
cat(sprintf("       Modality 2 (RNA/PCA new branch):     %s\n\n",
            paste(dim(rc[[2]]), collapse = "x")))

# --------------------------------------------------------------------------
# 7. Run the EXACT WNN code block from buildParentTree.R (post fix)
# --------------------------------------------------------------------------
n_modalities <- 2
reduction_coords_list <- rc
tmp <- matrix(stats::rnorm(n_cells * 3, 10), ncol = n_cells, nrow = 3)
colnames(tmp) <- rownames(rc[[1]])
rownames(tmp) <- paste0("t", seq_len(nrow(tmp)))
tmp_seurat <- CreateSeuratObject(tmp, min.cells = 0, min.features = 0, assay = "tmp")

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
cat(sprintf("[PASS] WNN fusion over ArchRProject reductions: snn %s\n\n",
            paste(dim(wnn@graphs$snn), collapse = "x")))

# --------------------------------------------------------------------------
# 8. Louvain on the fused SNN — should recover the planted 3 types
# --------------------------------------------------------------------------
cl <- FindClusters(wnn@graphs$snn, resolution = 0.5,
                   random.seed = 1, verbose = FALSE)[, 1]
cat("Cluster × truth table:\n")
print(table(cl, true_label))

ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  a_c <- sum(choose(rowSums(tab), 2))
  b_c <- sum(choose(colSums(tab), 2))
  t_c <- sum(choose(tab, 2))
  exp <- a_c * b_c / choose(n, 2)
  mx  <- 0.5 * (a_c + b_c)
  (t_c - exp) / (mx - exp)
}
ari_val <- ari(cl, true_label)
stopifnot(ari_val > 0.8)
cat(sprintf("\n[PASS] Louvain ARI vs planted labels = %.3f  (threshold > 0.8)\n\n",
            ari_val))

# --------------------------------------------------------------------------
# 9. Also verify .getMultiModalDistance works on ArchR-sourced reductions
# --------------------------------------------------------------------------
dist_mat <- .getMultiModalDistance(
  tmp_seurat,
  reduction_list = as.list(paste0("DR_", seq_len(n_modalities))),
  dim_list       = dim_list)
stopifnot(nrow(dist_mat) == n_cells && ncol(dist_mat) == n_cells)
cat(sprintf("[PASS] .getMultiModalDistance on ArchR object: shape %dx%d\n\n",
            nrow(dist_mat), ncol(dist_mat)))

cat("=============================================================\n")
cat("ArchRProject + integration_backend='seurat' smoke test PASSED.\n")
cat("=============================================================\n")
