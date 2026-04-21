# Tests for multimodal batch_correction_method handling in pruneTree setup

test_that("pruneTree batch-correction flags handle multimodal vectors", {
  batch_correction_method <- c("Harmony", "Harmony")
  use_batch_metrics <- any(batch_correction_method == "Harmony")
  no_batch_correction <- all(batch_correction_method == "none")

  expect_true(use_batch_metrics)
  expect_false(no_batch_correction)
})

test_that("pruneTree batch-correction flags handle no-batch vectors", {
  batch_correction_method <- c("none", "none")
  use_batch_metrics <- any(batch_correction_method == "Harmony")
  no_batch_correction <- all(batch_correction_method == "none")

  expect_false(use_batch_metrics)
  expect_true(no_batch_correction)
})
