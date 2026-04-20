# Tests for Leiden-compatible baseline resolution selection

source(testthat::test_path("..", "..", "R", "TreeUtils.R"), local = TRUE)

test_that("non-Leiden baseline resolution remains zero", {
  expect_equal(
    .getBaselineResolution(list(algorithm = 1)),
    0
  )
})

test_that("Leiden baseline resolution is positive", {
  expect_gt(
    .getBaselineResolution(list(algorithm = 4)),
    0
  )
})

test_that("missing algorithm defaults to zero baseline resolution", {
  expect_equal(
    .getBaselineResolution(list()),
    0
  )
})
