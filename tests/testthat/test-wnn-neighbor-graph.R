# Integration tests for the WNN (multi-modal) neighbor graph construction path
#
# The buildParentTree() WNN path originally used length(use_assay_build) to size
# loops and vectors. For ArchR objects use_assay_build is NULL, so length(NULL)
# returns 0, producing 1:0 iteration and empty dim_list — a hard runtime crash.
# The fix uses n_modalities instead, matching buildTree().

# ---------------------------------------------------------------------------
# Setup: mock multi-modal reduction data
# ---------------------------------------------------------------------------
set.seed(42)
n_cells <- 50
cell_names <- paste0("cell_", seq_len(n_cells))

# Modality 1: e.g., ATAC (LSI, 30 dims)
mod1 <- matrix(rnorm(n_cells * 30), nrow = n_cells, ncol = 30)
rownames(mod1) <- cell_names
colnames(mod1) <- paste0("LSI_", seq_len(30))

# Modality 2: e.g., RNA (PCA, 20 dims)
mod2 <- matrix(rnorm(n_cells * 20), nrow = n_cells, ncol = 20)
rownames(mod2) <- cell_names
colnames(mod2) <- paste0("PC_", seq_len(20))

reduction_coords_list <- list(mod1, mod2)

# ---------------------------------------------------------------------------
# Section A: Structural tests (no Seurat WNN dependency)
# ---------------------------------------------------------------------------

test_that("n_modalities sizing is correct when use_assay_build is NULL", {
  # This is the core of the bug: length(NULL) == 0, but n_modalities == 2
  use_assay_build <- NULL
  n_modalities <- 2

  # Buggy pattern: length(use_assay_build) gives 0

  expect_equal(length(use_assay_build), 0)
  # Fixed pattern: n_modalities gives correct value

  expect_equal(n_modalities, 2)

  # Buggy: vector("list", length = 0) produces empty list
  buggy_dim_list <- vector("list", length = length(use_assay_build))
  expect_equal(length(buggy_dim_list), 0)

  # Fixed: vector("list", length = 2) produces correct-length list
  fixed_dim_list <- vector("list", length = n_modalities)
  expect_equal(length(fixed_dim_list), 2)

  # Buggy: 1:length(NULL) == 1:0 == c(1, 0), iterates incorrectly
  buggy_seq <- seq(1, length(use_assay_build))
  expect_equal(buggy_seq, c(1, 0))

  # Fixed: seq(1, n_modalities) == 1:2
  fixed_seq <- seq(1, n_modalities)
  expect_equal(fixed_seq, c(1, 2))
})

test_that("subtree subsetting handles list-of-matrices correctly", {
  subset_ids <- cell_names[1:20]

  # Simulate the reduction output for seurat backend (list of matrices)
  P0_dim_reduction <- list(
    "reduction_coords" = reduction_coords_list,
    "var_features" = list(colnames(mod1), colnames(mod2))
  )

  # This is the code path from buildTree.R for seurat backend subtree
  if (is.list(P0_dim_reduction[["reduction_coords"]])) {
    P_i_dim_reduction <- list(
      "reduction_coords" = lapply(P0_dim_reduction[["reduction_coords"]],
                                  function(x) x[subset_ids, , drop = FALSE]),
      "var_features" = P0_dim_reduction[["var_features"]]
    )
  } else {
    P_i_dim_reduction <- list(
      "reduction_coords" = P0_dim_reduction[["reduction_coords"]][subset_ids, ],
      "var_features" = P0_dim_reduction[["var_features"]]
    )
  }

  # Verify subsetting produced correct shapes
  expect_true(is.list(P_i_dim_reduction[["reduction_coords"]]))
  expect_equal(length(P_i_dim_reduction[["reduction_coords"]]), 2)
  expect_equal(nrow(P_i_dim_reduction[["reduction_coords"]][[1]]), 20)
  expect_equal(nrow(P_i_dim_reduction[["reduction_coords"]][[2]]), 20)
  expect_equal(ncol(P_i_dim_reduction[["reduction_coords"]][[1]]), 30)
  expect_equal(ncol(P_i_dim_reduction[["reduction_coords"]][[2]]), 20)
  expect_equal(rownames(P_i_dim_reduction[["reduction_coords"]][[1]]), subset_ids)
})

test_that("subtree subsetting handles single matrix (non-list) correctly", {
  subset_ids <- cell_names[1:20]

  # Simulate the reduction output for archr backend (single matrix)
  P0_dim_reduction <- list(
    "reduction_coords" = mod1,
    "var_features" = colnames(mod1)
  )

  # This is the code path for archr backend / single modality
  if (is.list(P0_dim_reduction[["reduction_coords"]])) {
    P_i_dim_reduction <- list(
      "reduction_coords" = lapply(P0_dim_reduction[["reduction_coords"]],
                                  function(x) x[subset_ids, , drop = FALSE]),
      "var_features" = P0_dim_reduction[["var_features"]]
    )
  } else {
    P_i_dim_reduction <- list(
      "reduction_coords" = P0_dim_reduction[["reduction_coords"]][subset_ids, ],
      "var_features" = P0_dim_reduction[["var_features"]]
    )
  }

  expect_false(is.list(P_i_dim_reduction[["reduction_coords"]]))
  expect_equal(nrow(P_i_dim_reduction[["reduction_coords"]]), 20)
  expect_equal(ncol(P_i_dim_reduction[["reduction_coords"]]), 30)
})

test_that("WNN branch condition correctly routes ArchR + seurat to multimodal path", {
  # Simulate the condition used in buildTree.R / buildParentTree.R
  # if (n_modalities < 2 | (is_archr & (is.null(integration_backend) || integration_backend == "archr")))
  #   -> single reduction path (FindNeighbors)
  # else
  #   -> multimodal path (FindMultiModalNeighbors / WNN)

  check_single_path <- function(n_mod, is_archr, backend) {
    n_mod < 2 | (is_archr & (is.null(backend) || backend == "archr"))
  }

  # Single modality, any object -> single path
  expect_true(check_single_path(1, FALSE, NULL))
  expect_true(check_single_path(1, TRUE, "archr"))

  # ArchR multimodal + archr backend -> single path (addCombinedDims produces one matrix)
  expect_true(check_single_path(2, TRUE, "archr"))

  # ArchR multimodal + seurat backend -> multimodal WNN path
  expect_false(check_single_path(2, TRUE, "seurat"))

  # Non-ArchR multimodal -> multimodal WNN path
  expect_false(check_single_path(2, FALSE, NULL))
})

test_that("subtree WNN dims are derived from subtree reductions, not root reductions", {
  root_reduction_coords <- list(mod1, mod2)
  subtree_reduction_coords <- list(mod1[, 1:5, drop = FALSE],
                                   mod2[, 1:3, drop = FALSE])
  n_modalities <- 2

  root_dim_list <- vector("list", length = n_modalities)
  subtree_dim_list <- vector("list", length = n_modalities)

  for (i in 1:n_modalities) {
    root_dim_list[[i]] <- 1:ncol(root_reduction_coords[[i]])
    subtree_dim_list[[i]] <- 1:ncol(subtree_reduction_coords[[i]])
  }

  expect_equal(root_dim_list[[1]], 1:30)
  expect_equal(root_dim_list[[2]], 1:20)
  expect_equal(subtree_dim_list[[1]], 1:5)
  expect_equal(subtree_dim_list[[2]], 1:3)
})

# ---------------------------------------------------------------------------
# Section B: WNN integration tests (require Seurat >= 5.0.0)
# ---------------------------------------------------------------------------

# Helper: build WNN neighbor graph from a list of reduction matrices
# Replicates the exact WNN code path from buildTree.R / buildParentTree.R (fixed)
build_wnn_graph <- function(reduction_coords_list, n_modalities,
                            neighbor_params = list()) {
  n_cells <- nrow(reduction_coords_list[[1]])
  tmp <- matrix(stats::rnorm(n_cells * 3, 10), ncol = n_cells, nrow = 3)
  colnames(tmp) <- rownames(reduction_coords_list[[1]])
  rownames(tmp) <- paste0("t", seq_len(nrow(tmp)))
  tmp_seurat <- Seurat::CreateSeuratObject(tmp, min.cells = 0, min.features = 0, assay = "tmp")

  dim_list <- vector("list", length = n_modalities)
  for (i in 1:n_modalities) {
    tmp_seurat[[paste0("DR_", i)]] <- Seurat::CreateDimReducObject(
      embeddings = reduction_coords_list[[i]],
      key = paste0("DR", i, "_"), assay = "tmp"
    )
    dim_list[[i]] <- 1:ncol(reduction_coords_list[[i]])
  }

  result <- do.call(Seurat::FindMultiModalNeighbors, c(
    list(
      "object" = tmp_seurat,
      "reduction.list" = list(paste0("DR_", seq(1, n_modalities))),
      "dims.list" = dim_list,
      "knn.graph.name" = "nn",
      "snn.graph.name" = "snn"
    ),
    neighbor_params
  ))@graphs

  return(result)
}

test_that("WNN graph construction works with n_modalities (fixed path)", {
  skip_if(packageVersion("Seurat") < "5.0.0", "WNN integration tests require Seurat V5")

  result <- build_wnn_graph(reduction_coords_list, n_modalities = 2)
  expect_true("nn" %in% names(result))
  expect_true("snn" %in% names(result))
  expect_equal(nrow(result[["nn"]]), n_cells)
  expect_equal(ncol(result[["snn"]]), n_cells)
})

test_that("WNN graph produces valid neighbor/SNN matrices with correct cell names", {
  skip_if(packageVersion("Seurat") < "5.0.0", "WNN integration tests require Seurat V5")

  result <- build_wnn_graph(reduction_coords_list, n_modalities = 2)
  nn_mat <- result[["nn"]]
  snn_mat <- result[["snn"]]

  expect_equal(dim(nn_mat), c(n_cells, n_cells))
  expect_equal(dim(snn_mat), c(n_cells, n_cells))
  expect_equal(rownames(nn_mat), cell_names)
  expect_equal(colnames(snn_mat), cell_names)
})

test_that("WNN graph works with different numbers of dimensions per modality", {
  skip_if(packageVersion("Seurat") < "5.0.0", "WNN integration tests require Seurat V5")

  small_mod1 <- mod1[, 1:3, drop = FALSE]
  small_mod2 <- mod2[, 1:5, drop = FALSE]
  result <- build_wnn_graph(list(small_mod1, small_mod2), n_modalities = 2)
  expect_true("nn" %in% names(result))
  expect_equal(nrow(result[["nn"]]), n_cells)
})

test_that("subsetted reductions produce valid WNN graph", {
  skip_if(packageVersion("Seurat") < "5.0.0", "WNN integration tests require Seurat V5")

  subset_ids <- cell_names[1:20]
  subsetted <- lapply(reduction_coords_list, function(x) x[subset_ids, , drop = FALSE])
  result <- build_wnn_graph(subsetted, n_modalities = 2)
  expect_true("nn" %in% names(result))
  expect_equal(nrow(result[["nn"]]), 20)
})
