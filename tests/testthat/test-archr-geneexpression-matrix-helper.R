# Tests for ArchR GeneExpressionMatrix retrieval in CHOIR helpers

source(testthat::test_path("..", "..", "R", "HelperUtils.R"), local = TRUE)

if (!isClass("ArchRProject")) {
  setClass("ArchRProject", contains = "environment")
}

test_that("GeneExpressionMatrix helper path subsets by character feature names", {
  mock_archr <- new("ArchRProject")

  original_get_matrix <- get("getMatrixFromProject", envir = asNamespace("ArchR"))
  mock_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(GeneExpressionMatrix = Matrix::Matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(NULL, c("cell1", "cell2"))
    )),
    rowData = S4Vectors::DataFrame(name = c("geneA", "geneB"))
  )

  unlockBinding("getMatrixFromProject", asNamespace("ArchR"))
  assign("getMatrixFromProject", function(object, useMatrix) mock_se, envir = asNamespace("ArchR"))
  lockBinding("getMatrixFromProject", asNamespace("ArchR"))

  on.exit({
    unlockBinding("getMatrixFromProject", asNamespace("ArchR"))
    assign("getMatrixFromProject", original_get_matrix, envir = asNamespace("ArchR"))
    lockBinding("getMatrixFromProject", asNamespace("ArchR"))
  }, add = TRUE)

  result <- .getMatrix.ArchR(
    object = mock_archr,
    use_matrix = NULL,
    ArchR_matrix = "GeneExpressionMatrix",
    use_features = "geneB",
    exclude_features = NULL,
    use_cells = "cell2",
    verbose = FALSE
  )

  expect_equal(rownames(result), "geneB")
  expect_equal(colnames(result), "cell2")
  expect_equal(as.numeric(result[1, 1]), 4)
})
