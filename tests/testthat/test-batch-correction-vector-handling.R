# Tests for vector-valued batch_correction_method handling

source(testthat::test_path("..", "..", "R", "TreeUtils.R"), local = TRUE)

test_that("scalar non-Harmony is false", {
  expect_false(.usesHarmony("none"))
})

test_that("scalar Harmony is true", {
  expect_true(.usesHarmony("Harmony"))
})

test_that("multimodal vector with one Harmony is true", {
  expect_true(.usesHarmony(c("Harmony", "none")))
})

test_that("multimodal vector with no Harmony is false", {
  expect_false(.usesHarmony(c("none", "none")))
})
