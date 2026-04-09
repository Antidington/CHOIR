# Tests for integration_backend validation and related parameter relaxations
#
# These tests exercise .validInput() directly. We create a minimal mock
# ArchRProject class so that methods::is(obj, "ArchRProject") returns TRUE
# without requiring the full ArchR package.

# ---------------------------------------------------------------------------
# Setup: mock ArchRProject class
# ---------------------------------------------------------------------------

# Only define mock class if the real ArchR package is not available
if (!isClass("ArchRProject")) {
  setClass("ArchRProject", contains = "environment")
  mock_archr <- new("ArchRProject")
} else {
  # If ArchR is available, we still need a lightweight object
  # We will skip the tests that require a real ArchR object
  mock_archr <- NULL
}

# Create a simple mock Seurat-like object
mock_seurat <- list()
class(mock_seurat) <- "Seurat"

# ---------------------------------------------------------------------------
# Helper: call .validInput for integration_backend
# ---------------------------------------------------------------------------

validate_ib <- function(input, object, n_mod, arch_matrix, atac_vec,
                        red_method, norm_method, batch_method,
                        batch_labels, countsplit) {
  CHOIR:::.validInput(
    input = input,
    name = "integration_backend",
    other = list(object, n_mod, arch_matrix, atac_vec,
                 red_method, norm_method, batch_method,
                 batch_labels, countsplit)
  )
}

# ---------------------------------------------------------------------------
# Tests: integration_backend
# ---------------------------------------------------------------------------

test_that("non-ArchR object with integration_backend warns", {
  expect_warning(
    validate_ib("archr", mock_seurat, 1, NULL, FALSE, NULL, "none", "none", NULL, FALSE),
    "ignored for non-ArchR"
  )
})

test_that("non-ArchR object with NULL integration_backend is silent", {
  expect_silent(
    validate_ib(NULL, mock_seurat, 1, NULL, FALSE, NULL, "none", "none", NULL, FALSE)
  )
})

if (!is.null(mock_archr)) {

  test_that("ArchR object requires integration_backend", {
    expect_error(
      validate_ib(NULL, mock_archr, 1, "TileMatrix", TRUE, "IterativeLSI", "none", "none", NULL, FALSE),
      "required for ArchRProject"
    )
  })

  test_that("invalid integration_backend value errors", {
    expect_error(
      validate_ib("auto", mock_archr, 1, "TileMatrix", TRUE, "IterativeLSI", "none", "none", NULL, FALSE),
      "not among the permitted values"
    )
  })

  test_that("single-modality ArchR rejects seurat backend", {
    expect_error(
      validate_ib("seurat", mock_archr, 1, "TileMatrix", TRUE, "IterativeLSI", "none", "none", NULL, FALSE),
      "only support integration_backend = 'archr'"
    )
  })

  test_that("single-modality ArchR accepts archr backend", {
    expect_silent(
      validate_ib("archr", mock_archr, 1, "TileMatrix", TRUE, "IterativeLSI", "none", "none", NULL, FALSE)
    )
  })

  test_that("ArchR multimodal requires exactly 2 modalities", {
    expect_error(
      validate_ib("archr", mock_archr, 3,
                   c("TileMatrix", "GeneExpressionMatrix", "OtherMatrix"),
                   c(TRUE, FALSE, FALSE),
                   c("IterativeLSI", "PCA", "PCA"),
                   c("none", "LogNorm", "none"),
                   c("none", "none", "none"),
                   NULL, FALSE),
      "exactly 2 modalities"
    )
  })

  test_that("ArchR multimodal rejects broadcasting for atac", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   TRUE,  # length 1, not 2
                   c("IterativeLSI", "PCA"),
                   c("none", "LogNorm"),
                   c("none", "none"),
                   NULL, FALSE),
      "'atac' must be a vector of length 2"
    )
  })

  test_that("ArchR multimodal rejects broadcasting for reduction_method", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   "IterativeLSI",  # length 1, not 2
                   c("none", "LogNorm"),
                   c("none", "none"),
                   NULL, FALSE),
      "'reduction_method' must be a vector of length 2"
    )
  })

  test_that("ArchR multimodal requires one atac=TRUE and one atac=FALSE", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, TRUE),
                   c("IterativeLSI", "IterativeLSI"),
                   c("none", "none"),
                   c("none", "none"),
                   NULL, FALSE),
      "exactly one atac = TRUE and one atac = FALSE"
    )
  })

  test_that("ATAC modality must use TileMatrix", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("GeneScoreMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   c("IterativeLSI", "PCA"),
                   c("none", "LogNorm"),
                   c("none", "none"),
                   NULL, FALSE),
      "requires ArchR_matrix = 'TileMatrix'"
    )
  })

  test_that("RNA modality must use GeneExpressionMatrix", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "TileMatrix"),
                   c(TRUE, FALSE),
                   c("IterativeLSI", "PCA"),
                   c("none", "LogNorm"),
                   c("none", "none"),
                   NULL, FALSE),
      "requires ArchR_matrix = 'GeneExpressionMatrix'"
    )
  })

  test_that("ATAC modality must use IterativeLSI", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   c("PCA", "PCA"),
                   c("none", "LogNorm"),
                   c("none", "none"),
                   NULL, FALSE),
      "ATAC modality reduction_method must be 'IterativeLSI'"
    )
  })

  test_that("RNA modality must use PCA", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   c("IterativeLSI", "IterativeLSI"),
                   c("none", "none"),
                   c("none", "none"),
                   NULL, FALSE),
      "RNA modality reduction_method must be 'PCA'"
    )
  })

  test_that("LogNorm on ATAC modality errors", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   c("IterativeLSI", "PCA"),
                   c("LogNorm", "LogNorm"),
                   c("none", "none"),
                   NULL, FALSE),
      "only permitted for RNA"
    )
  })

  test_that("countsplit with ArchR multimodal errors", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   c("IterativeLSI", "PCA"),
                   c("none", "LogNorm"),
                   c("none", "none"),
                   NULL, TRUE),
      "not supported for ArchR multimodal"
    )
  })

  test_that("Harmony without batch_labels errors", {
    expect_error(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   c("IterativeLSI", "PCA"),
                   c("none", "LogNorm"),
                   c("Harmony", "Harmony"),
                   NULL, FALSE),
      "'batch_labels' is required"
    )
  })

  test_that("only one modality with Harmony warns", {
    expect_warning(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   c("IterativeLSI", "PCA"),
                   c("none", "LogNorm"),
                   c("Harmony", "none"),
                   "batch_col", FALSE),
      "Only one modality uses Harmony"
    )
  })

  test_that("valid archr multimodal config passes", {
    expect_silent(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   c("IterativeLSI", "PCA"),
                   c("none", "LogNorm"),
                   c("none", "none"),
                   NULL, FALSE)
    )
  })

  test_that("valid seurat multimodal config passes", {
    expect_silent(
      validate_ib("seurat", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   c("IterativeLSI", "PCA"),
                   c("none", "LogNorm"),
                   c("none", "none"),
                   NULL, FALSE)
    )
  })

  test_that("NULL reduction_method in multimodal passes (defaults applied later)", {
    expect_silent(
      validate_ib("archr", mock_archr, 2,
                   c("TileMatrix", "GeneExpressionMatrix"),
                   c(TRUE, FALSE),
                   NULL,
                   c("none", "LogNorm"),
                   c("none", "none"),
                   NULL, FALSE)
    )
  })
}

# ---------------------------------------------------------------------------
# Tests: relaxed normalization_method for ArchR
# ---------------------------------------------------------------------------

if (!is.null(mock_archr)) {
  test_that("ArchR normalization_method allows 'none'", {
    expect_silent(
      CHOIR:::.validInput("none", "normalization_method", list(mock_archr, 1, NULL))
    )
  })

  test_that("ArchR normalization_method allows 'LogNorm'", {
    expect_silent(
      CHOIR:::.validInput("LogNorm", "normalization_method", list(mock_archr, 1, NULL))
    )
  })

  test_that("ArchR normalization_method rejects 'SCTransform'", {
    expect_error(
      CHOIR:::.validInput("SCTransform", "normalization_method", list(mock_archr, 1, NULL)),
      "not among the permitted values"
    )
  })
}

# ---------------------------------------------------------------------------
# Tests: relaxed reduction_method for ArchR
# ---------------------------------------------------------------------------

if (!is.null(mock_archr)) {
  test_that("ArchR reduction_method allows 'IterativeLSI'", {
    expect_silent(
      CHOIR:::.validInput("IterativeLSI", "reduction_method", list(mock_archr, 1))
    )
  })

  test_that("ArchR reduction_method allows 'PCA'", {
    expect_silent(
      CHOIR:::.validInput("PCA", "reduction_method", list(mock_archr, 1))
    )
  })

  test_that("ArchR reduction_method rejects 'LSI'", {
    expect_error(
      CHOIR:::.validInput("LSI", "reduction_method", list(mock_archr, 1)),
      "not among the permitted values"
    )
  })
}
