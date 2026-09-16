test_that("get_example_data() returns the correct structure", {
  d <- get_example_data()
  expect_type(d, "list")
  expect_named(d, c("counts", "covariates", "metadata"), ignore.order = TRUE)

  expect_true(is.matrix(d$counts))
  expect_type(d$counts, "integer")
  expect_equal(nrow(d$counts), 20L)
  expect_gt(ncol(d$counts), 1L)
  expect_false(is.null(rownames(d$counts)))
  expect_false(is.null(colnames(d$counts)))

  expect_s3_class(d$covariates, "data.frame")
  expect_equal(nrow(d$covariates), nrow(d$counts))
  expect_named(d$covariates, c("Depth", "Distance_shore"), ignore.order = TRUE)

  expect_s3_class(d$metadata, "data.frame")
  expect_equal(nrow(d$metadata), nrow(d$counts))
  expect_true(all(c("sample_id", "TrueCommunity", "Depth", "Distance_shore") %in%
                    names(d$metadata)))
  expect_equal(d$metadata$sample_id, rownames(d$counts))
})

test_that("get_example_data() covariates go straight into eDNA_dmm()", {
  # Regression test for a real usability bug: `covariates` must be all-numeric
  # so the obvious call -- eDNA_dmm(counts = d$counts, covariates = d$covariates)
  # -- works without the user first dropping an ID or ground-truth column.
  d <- get_example_data()
  expect_true(all(vapply(d$covariates, is.numeric, logical(1))))
  expect_no_error(validate_covariates(d$covariates, d$counts, TRUE))
})

test_that("get_example_replicates() returns the correct structure", {
  r <- get_example_replicates()
  expect_type(r, "list")
  expect_named(r, c("counts", "replication", "covariates", "metadata"),
               ignore.order = TRUE)

  expect_true(is.matrix(r$counts))
  expect_type(r$counts, "integer")
  expect_type(r$replication, "character")
  expect_equal(length(r$replication), nrow(r$counts))

  n_stations <- length(unique(r$replication))
  expect_gt(nrow(r$counts), n_stations)          # more rows than stations
  expect_equal(nrow(r$covariates), n_stations)   # covariates are station-level
  expect_equal(nrow(r$metadata), n_stations)

  expect_true(all(vapply(r$covariates, is.numeric, logical(1))))

  # Replication is balanced at 3 bottles per station
  expect_true(all(table(r$replication) == 3L))
})

test_that("validate_counts() catches bad inputs", {
  # Non-matrix input
  expect_error(validate_counts("hello"), "must be a numeric matrix")

  # Fixtures need >= 3 rows and >= 2 cols to clear the dimension checks,
  # which run before any of the value checks below.

  # Negative values
  bad_counts <- matrix(c(1L, -1L, 2L, 3L, 4L, 5L), nrow = 3)
  expect_error(validate_counts(bad_counts), "negative values")

  # NA values
  na_counts <- matrix(c(1L, NA_integer_, 2L, 3L, 4L, 5L), nrow = 3)
  expect_error(validate_counts(na_counts), "NA value")

  # Empty sample (row 1 sums to zero)
  empty_counts <- matrix(c(0L, 1L, 2L, 0L, 4L, 5L), nrow = 3)
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

test_that("shipped raw CSVs split into the documented tables", {
  # Guards the split recipe in the Getting Started vignette: splitting the
  # shipped CSV must reproduce the objects get_example_data() returns.

  path <- system.file("extdata", "example_edna_raw.csv", package = "eDNAstructure")
  expect_true(nzchar(path))

  raw <- utils::read.csv(path, check.names = FALSE)
  expect_true(all(c("sample_id", "Depth", "Distance_shore") %in% names(raw)))

  taxon_cols <- grep("^Sp_", names(raw))
  expect_gt(length(taxon_cols), 1L)

  counts <- as.matrix(raw[, taxon_cols])
  rownames(counts) <- raw$sample_id
  storage.mode(counts) <- "integer"

  expect_identical(counts, get_example_data()$counts)

  covariates <- raw[, c("Depth", "Distance_shore")]
  expect_true(all(vapply(covariates, is.numeric, logical(1))))
})

test_that("shipped replicate CSV carries a replication column", {
  path <- system.file("extdata", "example_edna_replicates_raw.csv",
                      package = "eDNAstructure")
  expect_true(nzchar(path))

  raw <- utils::read.csv(path, check.names = FALSE)
  expect_true(all(c("replicate_id", "replication") %in% names(raw)))
  expect_equal(nrow(raw), 60L)
  expect_equal(length(unique(raw$replication)), 20L)

  # Covariates are constant within a station, so collapsing gives one row each
  cov <- unique(raw[, c("replication", "Depth", "Distance_shore")])
  expect_equal(nrow(cov), 20L)
})
