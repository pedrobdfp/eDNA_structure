# =============================================================================
# dmm_align_labels(): put every posterior draw on a common community labelling
# =============================================================================


#' Align community labels across posterior draws
#'
#' @description
#'
#' Permutes each posterior draw's community labels onto a common reference, and
#' rebuilds every community-indexed quantity from the aligned draws.
#'
#' @section Why any fit may need this:
#' A mixture likelihood is invariant to permuting its components, so the labels
#' carry no meaning and a sampler is free to renumber them at any point. Across
#' chains this is guaranteed. *Within* a chain it is rarer but not rare enough
#' to ignore: when two communities are similar, the barrier between their two
#' labellings is shallow and a chain can cross it mid-run.
#'
#' When that happens, `pi[3, ]` means one community for part of the run and a
#' different one for the rest, so the posterior mean composition is an average
#' of two communities and describes neither. The symptom is a high Rhat or a
#' low ESS on `pi` while permutation-invariant quantities such as `lp__` and
#' the log-likelihood look healthy -- a chain that sampled perfectly well and
#' merely renamed things partway through.
#'
#' Aligning the draws repairs it: the estimates become meaningful, and the
#' diagnostics start measuring sampling rather than naming.
#'
#' @section What alignment cannot repair:
#' Relabelling fixes disagreement about *names*. It cannot fix disagreement
#' about the *partition*. With several chains, two further numbers separate the
#' two cases, and both are reported:
#'
#' - `lp__` is invariant to permuting the components, so its Rhat ignores
#'   labelling entirely. A high value means the chains are sampling different
#'   regions of the posterior, full stop.
#' - `chain_agreement` is the share of samples that two chains assign to the
#'   same community *after* alignment, when the labels already mean the same
#'   thing. Values near one mean one mode found repeatedly; a block structure
#'   in the matrix means several modes, each found by some of the chains.
#'
#' When these two say the chains are in different modes, a high `pi` Rhat is
#' real and no amount of relabelling will lower it. The usual cause is a K the
#' data do not identify, and the fix is a different K, not a longer run. That
#' cuts both ways: K too large invites the classic overfit failure, but K too
#' small forces genuinely distinct communities to merge, which is itself
#' ambiguous about which samples go together and fails the same way. Compare
#' a range of K (e.g. [eDNA_loo()]) rather than assuming the fix is smaller.
#'
#' @param fit An `edna_dmm_fit` from [eDNA_dmm()], with
#'   any number of chains.
#' @param method `"STEPHENS"` (default) minimises KL divergence against the mean
#'   membership matrix, using the full membership probabilities.
#'   `"ECR-pivot"` anchors ECR to the allocation of the highest-density draw.
#'   `"ECR-iterative"` is faster but can settle into a local optimum in which
#'   two reference labellings coexist, leaving draws aligned within groups and
#'   mismatched between them; that looks exactly like residual non-convergence.
#' @param verbose Print progress. Default `TRUE`.
#'
#' @return The fit, with `pi_mean`, the membership columns of `sample_info`,
#'   `z_hat`, `assignment_certainty` and `beta_summary` recomputed from the
#'   aligned draws, plus an `alignment` element carrying the `method`, the
#'   `permutations`, `pct_permuted`, and `diagnostics` (Rhat and ESS for `pi`
#'   and `beta` **after** alignment, computed with the chain structure intact).
#'   With several chains it also carries `lp_rhat`, the permutation-invariant
#'   Rhat of the log posterior density, `chain_agreement`, the chain-by-chain
#'   share of samples assigned to the same community after alignment, and
#'   `min_agreement`, its smallest off-diagonal entry.
#'
#' @seealso [eDNA_dmm()], which calls this automatically after fitting.
#'
#' @examples
#' \dontrun{
#' fit     <- eDNA_dmm(counts, covariates, K = 8)
#' aligned <- dmm_align_labels(fit)
#'
#' # Did this chain relabel itself mid-run?
#' aligned$alignment$pct_permuted
#'
#' # How much did it matter?
#' max(abs(aligned$pi_mean - fit$pi_mean))
#' }
#'
#' @export
dmm_align_labels <- function(
    fit,
    method  = c("STEPHENS", "ECR-pivot", "ECR-iterative"),
    verbose = TRUE
) {
  check_stan_fit_object(fit, "dmm_align_labels")
  method <- match.arg(method)
  if (!requireNamespace("label.switching", quietly = TRUE))
    rlang::abort(c("Package `label.switching` is required.",
                   i = 'install.packages("label.switching")'))

  K <- fit$K

  # permuted = FALSE keeps [iteration, chain, parameter]. With permuted = TRUE
  # rstan merges the chains and shuffles the draw order, destroying exactly the
  # structure cross-chain Rhat is computed from.
  cp_a <- rstan::extract(fit$stan_fit, pars = "community_probs", permuted = FALSE)
  n_iter <- dim(cp_a)[1]; n_chain <- dim(cp_a)[2]; n_draws <- n_iter * n_chain
  N <- nrow(fit$sample_info)

  # Stan flattens matrix[N, K] column-major (n fastest), which R's own
  # column-major fill reproduces. Collapsing (iter, chain) puts chain slowest.
  cp <- array(cp_a, dim = c(n_iter, n_chain, N, K)); dim(cp) <- c(n_draws, N, K)
  z  <- apply(cp, c(1, 2), which.max)

  if (verbose) message(sprintf("Aligning %d draws (%d chain%s) by %s...",
                               n_draws, n_chain, if (n_chain > 1) "s" else "", method))

  perms <- switch(
    method,
    "STEPHENS"  = label.switching::stephens(cp)$permutations,
    "ECR-pivot" = {
      lp <- rstan::extract(fit$stan_fit, pars = "lp__", permuted = FALSE)
      label.switching::ecr(zpivot = z[which.max(as.vector(lp)), ],
                           z = z, K = K)$permutations
    },
    "ECR-iterative" = label.switching::ecr.iterative.1(z = z, K = K)$permutations)
  storage.mode(perms) <- "integer"
  pct_permuted <- 100 * mean(apply(perms, 1, function(p) !identical(p, seq_len(K))))
  if (verbose) message(sprintf("  %.1f%% of draws relabelled", pct_permuted))

  pi_a <- rstan::extract(fit$stan_fit, pars = "pi", permuted = FALSE)
  S    <- ncol(fit$pi_mean)
  pi_d <- array(pi_a, dim = c(n_iter, n_chain, K, S)); dim(pi_d) <- c(n_draws, K, S)

  beta_d <- if (length(fit$covariate_names) > 0) dmm_beta_draws(fit) else NULL
  if (!is.null(beta_d)) stopifnot(dim(beta_d)[1] == n_draws)

  for (i in seq_len(n_draws)) {
    p <- perms[i, ]
    pi_d[i, , ] <- pi_d[i, p, ]
    cp[i, , ]   <- cp[i, , p]
    if (!is.null(beta_d)) beta_d[i, , ] <- beta_d[i, p, ]
  }

  # The coefficients are identified by centring: within each covariate the K
  # values sum to zero. Permuting rows does not change a column sum, so the
  # constraint survives alignment, but it is reapplied here so that draws taken
  # from any source end up on the same scale and small numerical drift cannot
  # accumulate.
  if (!is.null(beta_d)) {
    for (j in seq_len(dim(beta_d)[3]))
      beta_d[, , j] <- beta_d[, , j] - rowMeans(beta_d[, , j])
  }

  fit$pi_mean <- apply(pi_d, c(2, 3), mean)
  dimnames(fit$pi_mean) <- list(paste0("Community ", seq_len(K)), fit$taxa_names)

  cp_mean <- apply(cp, c(2, 3), mean)
  prob_cols <- paste0("prob_comm", seq_len(K))
  fit$sample_info[prob_cols] <- as.data.frame(cp_mean)
  fit$sample_info$z_hat <- factor(paste0("Community ", apply(cp_mean, 1, which.max)),
                                  levels = paste0("Community ", seq_len(K)))
  fit$sample_info$assignment_certainty <- apply(cp_mean, 1, max)

  # Each parameter is folded back into [iteration, chain] so Rhat sees the
  # chain structure. With one chain this is split-Rhat, which still detects a
  # mid-run relabelling because the two halves disagree.
  # Per-parameter, then summarised. The summary reports the WORST parameter in
  # each block, not a single coefficient: "beta" covers every community x
  # covariate combination, so a bad row there names one coefficient out of
  # K x (P+1), and `diagnostics_full` says which.
  diag_full <- function(arr, label, labeller) {
    m <- matrix(arr, nrow = n_draws)
    st <- vapply(seq_len(ncol(m)), function(j) {
      dm <- matrix(m[, j], nrow = n_iter, ncol = n_chain)
      c(posterior::rhat(dm), posterior::ess_bulk(dm))
    }, numeric(2))
    data.frame(block = label, parameter = labeller(seq_len(ncol(m))),
               rhat = st[1, ], ess = st[2, ], stringsAsFactors = FALSE)
  }
  summarise_block <- function(d) {
    worst <- d$parameter[which.max(d$rhat)]
    data.frame(quantity = d$block[1], rhat = max(d$rhat, na.rm = TRUE),
               ess = min(d$ess, na.rm = TRUE),
               n_parameters = nrow(d), worst_parameter = worst,
               stringsAsFactors = FALSE)
  }

  # pi is K x S, flattened with k fastest.
  pi_lab <- function(j) sprintf("pi[community %d, %s]",
                                ((j - 1) %% K) + 1L,
                                fit$taxa_names[((j - 1) %/% K) + 1L])
  full <- diag_full(pi_d, "pi", pi_lab)

  if (!is.null(beta_d)) {
    fit$beta_draws_aligned <- beta_d
    fit$beta_summary <- .summarise_aligned_beta(beta_d, K, n_iter, n_chain)
    cov_labels <- dimnames(beta_d)[[3]]
    beta_lab <- function(j) sprintf("beta[community %d, %s]",
                                    ((j - 1) %% K) + 1L,
                                    cov_labels[((j - 1) %/% K) + 1L])
    full <- rbind(full, diag_full(beta_d, "beta", beta_lab))
  }
  dg <- do.call(rbind, lapply(split(full, full$block), summarise_block))
  dg <- dg[match(c("pi", "beta"), dg$quantity), , drop = FALSE]
  dg <- dg[!is.na(dg$quantity), , drop = FALSE]
  rownames(dg) <- NULL

  # Two checks that alignment cannot flatter, because neither depends on the
  # labels. lp__ is permutation invariant by construction. The agreement matrix
  # is computed AFTER alignment, when a community index means the same thing in
  # every chain, so what is left is disagreement about the partition itself.
  lp_a <- rstan::extract(fit$stan_fit, pars = "lp__", permuted = FALSE)
  lp_rhat <- posterior::rhat(matrix(lp_a, nrow = n_iter, ncol = n_chain))

  chain_agreement <- NULL
  if (n_chain > 1) {
    chain_of <- rep(seq_len(n_chain), each = n_iter)
    hard <- matrix(NA_integer_, N, n_chain)
    for (ci in seq_len(n_chain))
      hard[, ci] <- max.col(apply(cp[chain_of == ci, , , drop = FALSE],
                                  c(2, 3), mean), ties.method = "first")
    chain_agreement <- matrix(NA_real_, n_chain, n_chain,
                              dimnames = list(paste0("chain", seq_len(n_chain)),
                                              paste0("chain", seq_len(n_chain))))
    for (a in seq_len(n_chain)) for (b in seq_len(n_chain))
      chain_agreement[a, b] <- mean(hard[, a] == hard[, b])
    if (verbose)
      message(sprintf("  lp__ Rhat %.3f | chains agree on %.1f%%-%.1f%% of samples",
                      lp_rhat, 100 * min(chain_agreement[lower.tri(chain_agreement)]),
                      100 * max(chain_agreement[lower.tri(chain_agreement)])))
  }

  fit$alignment <- list(method = method, permutations = perms,
                        pct_permuted = pct_permuted, diagnostics = dg,
                        diagnostics_full = full, chains = n_chain,
                        lp_rhat = lp_rhat, chain_agreement = chain_agreement,
                        min_agreement = if (is.null(chain_agreement)) NA_real_
                                        else min(chain_agreement[lower.tri(chain_agreement)]))
  fit
}


# Summarise the aligned, centred coefficient draws.
#
# Every community is summarised: with a sum-to-zero constraint none of them is
# a reference, and each value is that community's deviation from the average
# community. `reliability` translates the effective sample size into the plain
# reading used by summary(), because a coefficient can be well estimated and
# still be badly sampled.
#' @keywords internal
.summarise_aligned_beta <- function(d, K, n_iter, n_chain) {
  cov_labels <- dimnames(d)[[3]]                  # "intercept" then covariates
  rows <- list()
  for (k in seq_len(K)) {
    for (j in seq_along(cov_labels)) {
      draws <- d[, k, j]                          # all draws, one coefficient
      ess   <- posterior::ess_bulk(
        matrix(draws, nrow = n_iter, ncol = n_chain))
      rows[[length(rows) + 1L]] <- data.frame(
        community     = paste0("Community ", k),
        covariate     = cov_labels[j],
        mean          = mean(draws),
        median        = stats::median(draws),
        ci_5          = stats::quantile(draws, 0.05, names = FALSE),
        ci_95         = stats::quantile(draws, 0.95, names = FALSE),
        ci_10         = stats::quantile(draws, 0.10, names = FALSE),
        ci_90         = stats::quantile(draws, 0.90, names = FALSE),
        prob_positive = mean(draws > 0),          # posterior P(coefficient > 0)
        prob_negative = mean(draws < 0),          # posterior P(coefficient < 0)
        ess           = ess,
        reliability   = if (ess > 400) "trustworthy"
                        else if (ess >= 100) "cautious" else "unreliable",
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
