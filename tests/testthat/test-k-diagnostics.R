# Fast, Stan-free regression tests for eDNA_dmm_k_diagnostics(). These use
# minimal fake fit objects rather than real eDNA_dmm() fits, so they run in
# milliseconds; they exist because the real bug they guard against, `fits`
# named "K2"/"K3"/... (as eDNA_loo() actually returns them) silently failing
# to resolve, was invisible without ever rendering the plot against real
# output.

fake_fit <- function(k, certainty = 0.95) {
  pm <- matrix(1 / 5, nrow = k, ncol = 5)
  pm[, 1] <- pm[, 1] + seq_len(k) * 0.01   # break ties, keep rows positive
  pm <- pm / rowSums(pm)
  list(
    pi_mean = pm,
    sample_info = data.frame(assignment_certainty = rep(certainty, 10))
  )
}

test_that("fits named K2/K3/... (eDNA_loo()'s own convention) resolve correctly", {
  fits <- list(K2 = fake_fit(2), K3 = fake_fit(3), K4 = fake_fit(4))

  d <- eDNA_dmm_k_diagnostics(
    fits        = fits,
    K_values    = 2:4,
    elpd_by_run = data.frame(K = rep(2:4, each = 2), elpd = rnorm(6, -100, 1))
  )

  # Every K should have contributed a row with real (non-NA) panel (e)/(f)
  # data. Before the fix, get_fit() looked up fits[["2"]] etc., always got
  # NULL, and these columns were silently all-NA / the merge dropped the rows.
  expect_true(all(c("min_distance", "mean_certainty") %in% names(d$table)))
  expect_equal(nrow(d$table), 3)
  expect_false(anyNA(d$table$min_distance))
  expect_false(anyNA(d$table$mean_certainty))
  expect_true(all(d$table$mean_certainty > 0.9))
})

test_that("fits named plainly as \"2\"/\"3\"/... still work", {
  fits <- list(`2` = fake_fit(2), `3` = fake_fit(3))

  d <- eDNA_dmm_k_diagnostics(
    fits        = fits,
    K_values    = 2:3,
    elpd_by_run = data.frame(K = rep(2:3, each = 2), elpd = rnorm(4, -100, 1))
  )

  expect_false(anyNA(d$table$min_distance))
})

test_that("K_values inferred from a K-prefixed fits list matches the K it names", {
  fits <- list(K2 = fake_fit(2), K5 = fake_fit(5))

  d <- eDNA_dmm_k_diagnostics(
    fits        = fits,
    elpd_by_run = data.frame(K = rep(c(2, 5), each = 2), elpd = rnorm(4, -100, 1))
  )

  expect_setequal(d$table$K, c(2, 5))
})

test_that("a full eDNA_loo()-shaped result unpacks automatically from `fits` alone", {
  loo_result <- list(
    fits         = list(K2 = fake_fit(2), K3 = fake_fit(3)),
    loo_table    = data.frame(K = 2:3, n_obs = c(20, 20), lp_rhat = c(1.0, 1.0)),
    loo_by_chain = data.frame(K = rep(2:3, each = 2), chain = rep(1:2, 2),
                              elpd = rnorm(4, -100, 1))
  )

  manual <- eDNA_dmm_k_diagnostics(
    fits        = loo_result$fits,
    elpd_by_run = loo_result$loo_by_chain,
    convergence = loo_result$loo_table
  )
  auto <- eDNA_dmm_k_diagnostics(loo_result)

  expect_identical(auto$table, manual$table)

  # An explicit elpd_by_run/convergence still wins over the ones inside
  # `fits`, rather than being silently overridden.
  override <- eDNA_dmm_k_diagnostics(
    loo_result,
    elpd_by_run = data.frame(K = 2:3, elpd = c(-1, -2))
  )
  expect_equal(nrow(override$table), 2)
})

test_that("panel (c) auto-builds from elpd_by_run$pareto_bad + convergence$n_obs", {
  fits <- list(K2 = fake_fit(2), K3 = fake_fit(3))
  elpd_by_run <- data.frame(
    K = rep(2:3, each = 2), chain = rep(1:2, 2),
    elpd = rnorm(4, -100, 1), pareto_bad = c(1, 2, 0, 1))
  convergence <- data.frame(K = 2:3, n_obs = c(20, 20), lp_rhat = c(1.0, 1.0))

  d <- eDNA_dmm_k_diagnostics(
    fits        = fits,
    elpd_by_run = elpd_by_run,
    convergence = convergence
  )

  expect_true("pct_pareto_bad" %in% names(d$table))
  expect_false(anyNA(d$table$pct_pareto_bad))
})
