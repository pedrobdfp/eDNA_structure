test_that("get_example_data() returns the correct structure", {
  d <- get_example_data()
  expect_type(d, "list")
  expect_named(d, c("counts", "covariates", "true_params"), ignore.order = TRUE)
  expect_true(is.matrix(d$counts))
  expect_equal(nrow(d$counts), 60L)
  expect_equal(ncol(d$counts), 25L)
  expect_s3_class(d$covariates, "data.frame")
  expect_equal(nrow(d$covariates), 60L)
  expect_true(all(c("sample_id", "latitude", "depth", "year", "depth_bin") %in%
                    names(d$covariates)))
})

test_that("validate_counts() catches bad inputs", {
  # Non-matrix input
  expect_error(validate_counts("hello"), "must be a numeric matrix")

  # Negative values
  bad_counts <- matrix(c(1L, -1L, 2L, 3L), nrow = 2)
  expect_error(validate_counts(bad_counts), "negative values")

  # NA values
  na_counts <- matrix(c(1L, NA_integer_, 2L, 3L), nrow = 2)
  expect_error(validate_counts(na_counts), "NA value")

  # Empty sample
  empty_counts <- matrix(c(0L, 0L, 2L, 3L), nrow = 2)
  expect_error(validate_counts(empty_counts), "zero total reads")

  # Too few samples
  tiny <- matrix(c(1L, 2L), nrow = 1)
  expect_error(validate_counts(tiny), "at least 3 samples")
})

test_that("validate_covariates() catches bad inputs", {
  good_counts <- matrix(rpois(30, 100), nrow = 6)
  good_counts <- good_counts + 1L

  # Wrong number of rows
  bad_cov <- matrix(1:10, nrow = 5)
  expect_error(validate_covariates(bad_cov, good_counts, TRUE), "same number of rows")

  # NA in covariates
  na_cov <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 6)
  expect_error(validate_covariates(na_cov, good_counts, TRUE), "NA value")

  # Zero variance covariate
  zero_var_cov <- matrix(rep(1, 12), nrow = 6)
  expect_error(validate_covariates(zero_var_cov, good_counts, TRUE), "zero variance")
})

test_that("validate_K() catches bad K values", {
  # Too small
  expect_error(validate_K(1, 10), "at least 2")

  # K >= N
  expect_error(validate_K(10, 10), "less than the number of samples")

  # Non-integer
  expect_error(validate_K(2.5, 10), "single positive integer")
})

test_that("make_community_colors() returns correct number of colors", {
  for (k in 2:7) {
    cols <- make_community_colors(k)
    expect_length(cols, k)
    expect_equal(names(cols), paste0("Community ", seq_len(k)))
  }
})
