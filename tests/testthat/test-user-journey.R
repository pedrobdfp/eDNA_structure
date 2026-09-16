# Tests written against what a NEW USER actually types, not against whatever
# the README happens to say today. Each block is a plausible first attempt.
# If one of these breaks, a first-time user hits an error on first contact.

test_that("Option A and Option B are interchangeable", {
  # The README offers get_example_data() and simulate_eDNA_survey() as two
  # routes to the same place. Downstream code must not care which was used.
  a <- get_example_data()
  b <- simulate_eDNA_survey(n_communities = 4, n_species = 40,
                            samples_per_community = 5, seed = 2026)

  for (d in list(a, b)) {
    expect_true(all(c("counts", "covariates", "metadata") %in% names(d)))

    # covariates: numeric only, straight into the model
    expect_true(all(vapply(d$covariates, is.numeric, logical(1))))
    expect_equal(nrow(d$covariates), nrow(d$counts))

    # metadata: carries the labels the plotting functions need
    expect_true("sample_id" %in% names(d$metadata))
    expect_equal(nrow(d$metadata), nrow(d$counts))

    # the identical plotting call works for both
    expect_no_error(
      plot_true_compositions(d$counts, metadata = d$metadata,
                             facet_var = "TrueCommunity")
    )
  }
})

test_that("metadata rows are matched without a sample_id column", {
  d <- get_example_data()

  # 1. covariates has no sample_id column, but has row names -> resolves
  expect_no_error(plot_true_compositions(d$counts, metadata = d$covariates))
  expect_no_error(eDNA_dmm_structure_metadata_ok <- TRUE)

  # 2. no sample_id column and no row names -> positional, with a message
  plain <- d$covariates
  rownames(plain) <- NULL
  expect_message(
    plot_true_compositions(d$counts, metadata = plain),
    "matching rows to samples by position"
  )

  # 3. a differently named ID column, declared explicitly
  renamed <- d$metadata
  names(renamed)[names(renamed) == "sample_id"] <- "station"
  expect_no_error(
    plot_true_compositions(d$counts, metadata = renamed,
                           sample_id_col = "station",
                           facet_var = "TrueCommunity")
  )
})

test_that("unresolvable metadata gives an actionable error", {
  d <- get_example_data()
  bad <- data.frame(site = paste0("X", 1:7), temperature = rnorm(7))

  expect_error(
    plot_true_compositions(d$counts, metadata = bad),
    "Can't tell which sample each row of `metadata` describes"
  )
  # The message must name the shape contract and the escape hatch
  err <- tryCatch(plot_true_compositions(d$counts, metadata = bad),
                  error = function(e) conditionMessage(e))
  expect_true(grepl("samples (rows) x taxa (columns)", err, fixed = TRUE))
  expect_true(grepl("sample_id_col", err, fixed = TRUE))
  expect_true(grepl("site, temperature", err, fixed = TRUE))   # lists what was actually there
})

test_that("covariates from either option go straight into eDNA_dmm()", {
  a <- get_example_data()
  b <- simulate_eDNA_survey(n_communities = 4, n_species = 40,
                            samples_per_community = 5, seed = 2026)
  for (d in list(a, b)) {
    expect_no_error(validate_covariates(d$covariates, d$counts, TRUE))
  }
})

test_that("simulate_eDNA_survey() keeps its ground-truth fields", {
  # The slimming must not remove what the validation workflow depends on.
  b <- simulate_eDNA_survey(n_communities = 4, n_species = 40,
                            samples_per_community = 5, seed = 2026)
  expect_true(all(c("community_compositions", "metab_df", "sample_metadata",
                    "contributors") %in% names(b)))
  expect_false("sample_id" %in% names(b$covariates))       # ID kept out
  expect_false("TrueCommunity" %in% names(b$covariates))   # truth kept out
  expect_true("TrueCommunity" %in% names(b$metadata))      # but available
})
