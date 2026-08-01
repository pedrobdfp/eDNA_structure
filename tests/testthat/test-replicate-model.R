# =============================================================================
# Replicate-aware model: station_id handling
# =============================================================================
# The fast tests here cover validation and reshaping only — no Stan fitting, so
# they run in the normal test suite. The fitting tests at the bottom are opt-in
# via EDNA_RUN_SLOW_TESTS=true because each one compiles and samples a model.

make_counts <- function(n_rows, n_taxa = 5, seed = 1) {
  set.seed(seed)
  m <- matrix(rpois(n_rows * n_taxa, 200) + 1L, nrow = n_rows)
  storage.mode(m) <- "integer"
  m
}

# ── validate_station_id() ────────────────────────────────────────────────────

test_that("validate_station_id() maps rows to stations in first-appearance order", {
  counts <- make_counts(6)
  st <- validate_station_id(c("B", "B", "A", "A", "C", "C"), counts)

  expect_equal(st$levels, c("B", "A", "C"))
  expect_equal(st$index, c(1L, 1L, 2L, 2L, 3L, 3L))
  expect_equal(st$n_stations, 3L)
  expect_false(st$fallback)
  expect_equal(as.integer(st$reps), c(2L, 2L, 2L))
})

test_that("validate_station_id() handles ragged replication", {
  counts <- make_counts(6)
  st <- validate_station_id(c("A", "A", "A", "B", "C", "C"), counts)

  expect_equal(st$n_stations, 3L)
  expect_equal(as.integer(st$reps), c(3L, 1L, 2L))
  expect_false(st$fallback)
})

test_that("validate_station_id() accepts factors and numeric labels", {
  counts <- make_counts(4)
  st_fac <- validate_station_id(factor(c("a", "a", "b", "b")), counts)
  st_num <- validate_station_id(c(10, 10, 20, 20), counts)

  expect_equal(st_fac$index, c(1L, 1L, 2L, 2L))
  expect_equal(st_num$index, c(1L, 1L, 2L, 2L))
  expect_equal(st_num$levels, c("10", "20"))
})

test_that("validate_station_id() rejects malformed input", {
  counts <- make_counts(6)

  expect_error(validate_station_id(c("A", "A", "B"), counts), "length 3")
  expect_error(validate_station_id(c("A", "A", NA, "B", "B", "C"), counts), "NA value")
  expect_error(validate_station_id(rep("A", 6), counts), "only 1 station")
  expect_error(
    validate_station_id(data.frame(a = 1:6, b = 1:6), counts),
    "must be a vector"
  )
})

test_that("validate_station_id() signals fallback when nothing is replicated", {
  counts <- make_counts(4)
  st <- validate_station_id(c("A", "B", "C", "D"), counts)

  # Not an error and not a warning — unreplicated data stays fully supported,
  # the caller just uses the standard model instead.
  expect_true(st$fallback)
  expect_equal(st$n_stations, 4L)
})

test_that("validate_station_id() notes when very few stations are replicated", {
  counts <- make_counts(7)
  expect_message(
    validate_station_id(c("A", "A", "B", "C", "D", "E", "F"), counts),
    "more than one replicate"
  )
})

# ── collapse_station_covariates() ────────────────────────────────────────────

test_that("collapse_station_covariates() passes station-level covariates through", {
  counts <- make_counts(6)
  st  <- validate_station_id(c("A", "A", "B", "B", "C", "C"), counts)
  cov <- data.frame(depth = c(10, 20, 30))

  expect_equal(collapse_station_covariates(cov, st), cov)
  expect_null(collapse_station_covariates(NULL, st))
})

test_that("collapse_station_covariates() collapses replicate-level covariates", {
  counts <- make_counts(6)
  st  <- validate_station_id(c("A", "A", "B", "B", "C", "C"), counts)
  cov <- data.frame(depth = c(10, 10, 20, 20, 30, 30),
                    lat   = c(1, 1, 2, 2, 3, 3))

  out <- collapse_station_covariates(cov, st)
  expect_equal(nrow(out), 3L)
  expect_equal(out$depth, c(10, 20, 30))
  expect_equal(out$lat, c(1, 2, 3))
})

test_that("collapse_station_covariates() collapses matrices too", {
  counts <- make_counts(6)
  st  <- validate_station_id(c("A", "A", "B", "B", "C", "C"), counts)
  cov <- cbind(depth = c(10, 10, 20, 20, 30, 30))

  out <- collapse_station_covariates(cov, st)
  expect_equal(nrow(out), 3L)
  expect_equal(as.numeric(out[, "depth"]), c(10, 20, 30))
})

test_that("collapse_station_covariates() rejects covariates varying within a station", {
  counts <- make_counts(6)
  st  <- validate_station_id(c("A", "A", "B", "B", "C", "C"), counts)
  cov <- data.frame(depth = c(10, 999, 20, 20, 30, 30))

  expect_error(collapse_station_covariates(cov, st), "vary between replicates")
})

test_that("collapse_station_covariates() rejects an unusable row count", {
  counts <- make_counts(6)
  st  <- validate_station_id(c("A", "A", "B", "B", "C", "C"), counts)

  expect_error(
    collapse_station_covariates(data.frame(depth = 1:4), st),
    "matches neither"
  )
})

# ── Reshaping simulated replicates ───────────────────────────────────────────

test_that("replicate rows from metab_df sum to the collapsed count matrix", {
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")

  sim <- simulate_eDNA_survey(
    n_communities         = 2,
    n_species             = 6,
    samples_per_community = 3,
    bio_reps              = 3,
    seq_reps              = 1,
    mean_read_depth       = 2000,
    seed                  = 99
  )

  # One row per (station, replicate)
  rep_long <- sim$metab_df
  rep_long$station <- paste0("STN_", sprintf("%03d", rep_long$SampleID))
  rep_long$rep_key <- paste0(rep_long$station, "_B", rep_long$BioRep,
                             "_S", rep_long$SeqRep)

  rep_wide <- tidyr::pivot_wider(
    rep_long[, c("rep_key", "station", "Species", "Counts")],
    names_from  = "Species",
    values_from = "Counts",
    values_fill = 0L
  )

  station_id <- rep_wide$station
  rep_counts <- as.matrix(rep_wide[, !(names(rep_wide) %in% c("rep_key", "station"))])

  # The whole point: replicates are separate rows, and their per-station sums
  # reproduce exactly the matrix the summed workflow would have built.
  expect_gt(nrow(rep_counts), nrow(sim$counts))
  summed <- rowsum(rep_counts, group = station_id)
  summed <- summed[rownames(sim$counts), colnames(summed), drop = FALSE]

  expect_equal(
    as.vector(summed[, order(colnames(summed))]),
    as.vector(sim$counts[, order(sub("^Sp_", "", colnames(sim$counts)))]),
    ignore_attr = TRUE
  )
})

# ── Fitting (opt-in: these compile and sample Stan models) ───────────────────

skip_unless_slow <- function() {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("EDNA_RUN_SLOW_TESTS"), "true"),
    "Set EDNA_RUN_SLOW_TESTS=true to run Stan fitting tests."
  )
}

test_that("eDNA_dmm() falls back to the standard model when nothing is replicated", {
  skip_unless_slow()

  sim <- simulate_eDNA_survey(n_communities = 2, n_species = 6,
                              samples_per_community = 4, seed = 7)
  cov <- sim$covariates[, c("Depth", "Distance_shore")]

  expect_message(
    fit <- eDNA_dmm(sim$counts, covariates = cov, K = 2,
                    station_id = rownames(sim$counts),
                    iter = 400, warmup = 200, verbose = TRUE),
    "exactly 1 replicate"
  )

  expect_false(fit$replicate)
  expect_null(fit$theta_mean)
  expect_null(fit$phi_mean)
  expect_equal(fit$N, fit$R)
})

test_that("eDNA_dmm() fits the replicate model and returns station compositions", {
  skip_unless_slow()
  skip_if_not_installed("tidyr")

  sim <- simulate_eDNA_survey(
    n_communities = 2, n_species = 6, samples_per_community = 4,
    bio_reps = 3, mean_read_depth = 2000, seed = 11
  )

  rep_long <- sim$metab_df
  rep_long$station <- paste0("STN_", sprintf("%03d", rep_long$SampleID))
  rep_long$rep_key <- paste0(rep_long$station, "_B", rep_long$BioRep,
                             "_S", rep_long$SeqRep)
  rep_wide <- tidyr::pivot_wider(
    rep_long[, c("rep_key", "station", "Species", "Counts")],
    names_from = "Species", values_from = "Counts", values_fill = 0L
  )
  station_id <- rep_wide$station
  rep_counts <- as.matrix(rep_wide[, !(names(rep_wide) %in% c("rep_key", "station"))])

  station_cov <- sim$covariates[match(unique(station_id), sim$covariates$sample_id),
                                c("Depth", "Distance_shore")]

  fit <- eDNA_dmm(rep_counts, covariates = station_cov, K = 2,
                  station_id = station_id,
                  iter = 600, warmup = 300, verbose = FALSE)

  expect_true(fit$replicate)
  expect_equal(fit$N, length(unique(station_id)))
  expect_equal(fit$R, nrow(rep_counts))

  # theta is a composition per station
  expect_equal(dim(fit$theta_mean), c(fit$N, fit$S))
  expect_equal(unname(rowSums(fit$theta_mean)), rep(1, fit$N), tolerance = 1e-6)

  # phi is estimated and distinct from alpha
  expect_true(is.finite(fit$phi_mean))
  expect_gt(fit$phi_mean, 0)

  # sample_info stays one row per station, so plotting code is unaffected
  expect_equal(nrow(fit$sample_info), fit$N)
  expect_true(all(paste0("prob_comm", seq_len(2)) %in% names(fit$sample_info)))
})
