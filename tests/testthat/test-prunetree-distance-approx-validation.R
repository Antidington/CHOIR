# Tests for pruneTree distance_approx validation context

source(testthat::test_path("..", "..", "R", "ValidationUtils.R"), local = TRUE)

if (!isClass("ArchRProject")) {
  setClass("ArchRProject", contains = "environment")
}

test_that("distance_approx validator accepts 4-element pruneTree context", {
  mock_archr <- new("ArchRProject")

  expect_error(
    .validInput(TRUE, "distance_approx", list(1000, mock_archr, 2, "seurat")),
    "incompatible with ArchR multimodal data using integration_backend = 'seurat'"
  )
})
