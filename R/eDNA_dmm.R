# =============================================================================
# eDNA_dmm(): Fit a Dirichlet-Multinomial Mixture model
# =============================================================================

#' Fit a Dirichlet-Multinomial Mixture model to eDNA count data
#'
#' @description
#' `eDNA_dmm()` fits a Bayesian Dirichlet-Multinomial Mixture (DMM) model to a
#' sample by taxon read count matrix. It identifies `K` latent ecological
#' communities, estimates their taxonomic compositions, and models how
#' environmental covariates such as depth or latitude drive community
#' membership through a softmax regression.
#'
#' Inference is performed with Stan via [rstan::sampling()]. Several chains are
#' run by default and their community labels are aligned automatically before
#' anything is summarised, so the returned estimates and convergence
#' diagnostics are directly interpretable. The result is an `edna_dmm_fit`
#' object that can be passed to [eDNA_dmm_structure()], [eDNA_dmm_nmds()],
#' [eDNA_dmm_compositions()], [eDNA_dmm_beta()] and
#' [eDNA_dmm_beta_intervals()].
#'
#' @section Input format:
#' `counts` should be a sample by taxon matrix of non-negative integer read
#' counts. Rows are samples and columns are taxa or ASVs; the model does not
#' require taxonomic annotation. Row names are used as sample identifiers in
#' plots and column names as taxon labels. Taxa with no reads in any sample are
#' dropped.
#'
#' `covariates` should be a sample by covariate matrix or data frame with the
#' same number of rows as `counts`, in the same order. Covariates are
#' standardised internally by default (`scale_covariates = TRUE`), which is
#' strongly recommended for interpretable coefficients and good mixing. Pass
#' `NULL` to fit an intercept-only model.
#'
#' @section The model:
#' For each sample `i` the model marginalises over a latent community
#' assignment `z_i`:
#'
#' - Community compositions: `pi_k ~ Dirichlet(conc * 1_S)` for each community
#'   `k`.
#' - Community membership: `P(z_i = k) = softmax(x_i' beta_k)`, where `x_i` is
#'   the covariate row for sample `i` with an intercept prepended.
#' - Observed counts: `x_i | z_i = k ~ DirichletMultinomial(N_i, alpha * pi_k)`.
#'
#' `alpha` is a global overdispersion scalar estimated from the data. Values
#' much greater than one approach a multinomial; values near one indicate the
#' strong overdispersion typical of eDNA read counts.
#'
#' @section Identifying the softmax coefficients:
#' Adding the same constant to every community's linear predictor leaves the
#' softmax unchanged, so for each covariate the `K` coefficients hold one
#' redundant degree of freedom. The model removes it by constraining the
#' coefficients to sum to zero across communities within each covariate. Each
#' coefficient is therefore the deviation of that community from the average
#' community, and `beta` has `K` rows rather than `K - 1`.
#'
#' This keeps the prior exchangeable across communities, which matches the fact
#' that mixture components have no inherent identity, and leaves the posterior
#' less correlated and easier to sample. To read coefficients against one
#' particular community instead, pass `reference` to [dmm_beta_draws()] or
#' [eDNA_dmm_beta_intervals()]; that is a change of coordinates applied to the
#' finished draws and requires no refitting.
#'
#' @section Why the labels are aligned:
#' A mixture likelihood is invariant to permuting its components, so chains
#' that have converged on the same posterior can still disagree about which
#' community is labelled 1, 2, 3 and so on, and a single chain can renumber its
#' components partway through a run. Estimates and cross-chain diagnostics
#' computed in that state describe the labelling rather than the fit.
#'
#' `eDNA_dmm()` therefore aligns every posterior draw onto a common labelling
#' with [dmm_align_labels()] before computing anything, using the algorithm of
#' Stephens (2000). Community compositions, membership probabilities and
#' coefficients are all rebuilt from the aligned draws.
#'
#' Alignment repairs disagreement about names. It cannot repair disagreement
#' about the partition. Two quantities that relabelling cannot affect are
#' reported alongside it in `fit$alignment`: `lp_rhat`, the convergence
#' statistic of the log posterior density, which is invariant to permutation by
#' construction, and `min_agreement`, the share of samples that the two least
#' similar chains place in the same community. Values close to one for both
#' mean the chains found the same solution; a clear departure in either means
#' they did not, which usually indicates a `K` the data do not support.
#'
#' @param counts Integer matrix of read counts, samples in rows and taxa in
#'   columns.
#' @param covariates Matrix or data frame of covariates, one row per sample, or
#'   `NULL` for an intercept-only model.
#' @param K Number of latent communities. Default `2`.
#' @param scale_covariates Standardise covariates internally. Default `TRUE`.
#' @param chains Number of chains. Default `4`.
#' @param cores Cores for running chains in parallel. `NULL` (default) uses one
#'   per chain, capped at the number of physical cores.
#' @param method Label alignment algorithm, passed to [dmm_align_labels()].
#'   `"STEPHENS"` (default) minimises Kullback-Leibler divergence against the
#'   mean membership matrix. `"ECR-pivot"` anchors to the allocation of the
#'   highest-density draw. `"ECR-iterative"` is faster but can settle into a
#'   local optimum that mimics non-convergence, so prefer the first two.
#' @param iter Total iterations per chain. Default `4000`.
#' @param warmup Warmup iterations per chain. Default `2000`.
#' @param adapt_delta Target acceptance rate. Raise towards `0.99` if divergent
#'   transitions appear. Default `0.95`.
#' @param max_treedepth Maximum tree depth. Default `12`.
#' @param seed Random seed. Default `13`.
#' @param conc Dirichlet concentration for the community compositions. Values
#'   below one favour sparse compositions, which suits eDNA. Default `0.5`.
#' @param alpha_shape,alpha_rate Gamma prior on the overdispersion `alpha`; the
#'   prior mean is `alpha_shape / alpha_rate`. Defaults `5` and `2`.
#' @param rhat_threshold,ess_threshold Thresholds above and below which a
#'   warning is issued for quantities that remain problematic after alignment.
#'   Defaults `1.05` and `100`.
#' @param verbose Print progress and diagnostics. Default `TRUE`.
#'
#' @return An object of class `edna_dmm_fit`, a list containing:
#' \describe{
#'   \item{`stan_fit`}{The underlying [rstan::stanfit-class] object.}
#'   \item{`sample_info`}{One row per sample: `sample_id`, the membership
#'     probability of each community, the most probable community `z_hat`, and
#'     `assignment_certainty`, the largest membership probability.}
#'   \item{`pi_mean`}{`K` by `S` matrix of posterior mean compositions.}
#'   \item{`beta_summary`}{One row per community and covariate, summarising the
#'     centred coefficients.}
#'   \item{`beta_draws_aligned`}{`draws` by `K` by `P + 1` array of aligned,
#'     centred coefficient draws.}
#'   \item{`alignment`}{Alignment method, the permutations applied, the share
#'     of draws relabelled, post-alignment diagnostics, `lp_rhat`,
#'     `chain_agreement` and `min_agreement`.}
#'   \item{`alpha_mean`}{Posterior mean overdispersion.}
#'   \item{`K`,`N`,`S`}{Communities, samples and taxa.}
#'   \item{`taxa_names`,`covariate_names`,`scale_info`,`counts`,`stan_data`,`call`}{
#'     Inputs and metadata retained for downstream functions.}
#' }
#'
#' @references
#' Stephens, M. (2000). Dealing with label switching in mixture models.
#' *Journal of the Royal Statistical Society B*, 62(4), 795-809.
#'
#' @seealso [eDNA_loo()] to compare models across `K`,
#'   [eDNA_dmm_k_diagnostics()] for the wider set of `K` diagnostics,
#'   [dmm_align_labels()] for the alignment step on its own.
#'
#' @examples
#' \dontrun{
#' d   <- get_example_data()
#' fit <- eDNA_dmm(d$counts, d$covariates, K = 4)
#'
#' fit                              # dimensions and community sizes
#' summary(fit)                     # convergence and covariate effects
#' fit$alignment$min_agreement      # did the chains find the same solution?
#'
#' eDNA_dmm_structure(fit, metadata = d$metadata)
#' eDNA_dmm_beta_intervals(fit)$plot
#' }
#'
#' @export
eDNA_dmm <- function(
    counts,
    covariates       = NULL,
    K                = 2,
    scale_covariates = TRUE,
    chains           = 4,
    cores            = NULL,
    method           = c("STEPHENS", "ECR-pivot", "ECR-iterative"),
    iter             = 4000,
    warmup           = 2000,
    adapt_delta      = 0.95,
    max_treedepth    = 12,
    seed             = 13,
    conc             = 0.5,
    alpha_shape      = 5,
    alpha_rate       = 2,
    rhat_threshold   = 1.05,
    ess_threshold    = 100,
    verbose          = TRUE
) {
  cl     <- match.call()
  method <- match.arg(method)

  if (!requireNamespace("label.switching", quietly = TRUE))
    rlang::abort(c(
      "Package `label.switching` is required for community label alignment.",
      i = 'Install it with install.packages("label.switching").'
    ))

  # ── Input validation ────────────────────────────────────────────────────────
  counts <- validate_counts(counts)
  N      <- nrow(counts)
  K      <- validate_K(K, N)
  counts <- counts[, colSums(counts) > 0, drop = FALSE]   # drop all-zero taxa
  covariates <- validate_covariates(covariates, counts, scale_covariates,
                                    n_expected = N)

  if (!is.numeric(chains) || chains < 1 || chains != round(chains))
    rlang::abort("`chains` must be a positive integer.")
  if (!is.numeric(iter) || iter < 100)
    rlang::abort("`iter` must be a positive integer >= 100.")
  if (!is.numeric(warmup) || warmup < 1 || warmup >= iter)
    rlang::abort("`warmup` must be a positive integer less than `iter`.")
  if (!is.numeric(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1)
    rlang::abort("`adapt_delta` must be strictly between 0 and 1.")
  if (!is.numeric(conc) || conc <= 0)
    rlang::abort("`conc` must be a positive number.")
  if (!is.numeric(alpha_shape) || alpha_shape <= 0)
    rlang::abort("`alpha_shape` must be a positive number.")
  if (!is.numeric(alpha_rate) || alpha_rate <= 0)
    rlang::abort("`alpha_rate` must be a positive number.")

  # ── Covariate scaling, recorded before the attributes are stripped ──────────
  scale_info <- NULL
  P <- ncol(covariates)
  if (scale_covariates && P > 0) {
    scale_info <- list(center = attr(covariates, "scale_center"),
                       scale  = attr(covariates, "scale_scale"))
    attr(covariates, "scale_center") <- NULL   # Stan rejects extra attributes
    attr(covariates, "scale_scale")  <- NULL
  }

  S               <- ncol(counts)
  taxa_names      <- colnames(counts)
  covariate_names <- if (P > 0) colnames(covariates) else character(0)
  sample_ids      <- rownames(counts)

  # ── Stan data ───────────────────────────────────────────────────────────────
  stan_data <- list(
    N           = N,                    # number of samples
    S           = S,                    # number of taxa
    K           = K,                    # number of communities
    P           = P,                    # number of covariates
    X           = counts,               # count matrix, samples by taxa
    covariates  = if (P > 0) covariates else matrix(numeric(0), nrow = N, ncol = 0),
    conc        = conc,                 # Dirichlet concentration on each pi_k
    alpha_shape = alpha_shape,          # Gamma prior shape for alpha
    alpha_rate  = alpha_rate            # Gamma prior rate for alpha
  )

  # ── Sampling ────────────────────────────────────────────────────────────────
  # Chains are independent, so one core each by default. Running them
  # sequentially multiplies wall time by `chains` for no benefit.
  if (is.null(cores))
    cores <- max(1L, min(as.integer(chains), parallel::detectCores(logical = FALSE)))

  if (verbose) {
    message(sprintf("Fitting DMM: N=%d samples, S=%d taxa, K=%d communities, P=%d covariates",
                    N, S, K, P))
    message(sprintf("MCMC: %d chain%s on %d core%s, %d iterations (%d warmup)",
                    chains, if (chains > 1) "s" else "",
                    cores,  if (cores  > 1) "s" else "", iter, warmup))
  }

  stan_fit <- withCallingHandlers(
    rstan::sampling(
      .get_stanmodel("dmm"), data = stan_data,
      chains = chains, cores = cores,
      iter = iter, warmup = warmup, seed = seed,
      refresh = if (verbose) max(1, floor(iter / 20)) else 0,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)),
    warning = function(w) {
      # Convergence warnings raised here describe the unaligned draws, in which
      # community-indexed parameters still carry each chain's own arbitrary
      # labelling. They are expected, and alignment below is what resolves
      # them. Anything genuine is reported again afterwards.
      if (grepl("R-?hat|Bulk|Tail|ESS|divergent|convergence", conditionMessage(w),
                ignore.case = TRUE))
        invokeRestart("muffleWarning")
    })

  if (verbose) rstan::check_hmc_diagnostics(stan_fit)

  # ── Assemble the fit ────────────────────────────────────────────────────────
  # Community-indexed values here are provisional: they are computed before
  # alignment, and dmm_align_labels() replaces every one of them below.
  cp_mean <- apply(rstan::extract(stan_fit, pars = "community_probs")[[1]],
                   c(2, 3), mean)
  pi_mean <- apply(rstan::extract(stan_fit, pars = "pi")[[1]], c(2, 3), mean)
  colnames(pi_mean) <- taxa_names
  rownames(pi_mean) <- paste0("Community ", seq_len(K))

  prob_df <- as.data.frame(cp_mean)
  colnames(prob_df) <- paste0("prob_comm", seq_len(K))
  sample_info <- cbind(
    data.frame(sample_id = if (!is.null(sample_ids)) sample_ids
                           else paste0("Sample_", seq_len(N)),
               stringsAsFactors = FALSE),
    prob_df)
  sample_info$z_hat <- factor(
    apply(cp_mean, 1, which.max), levels = seq_len(K),
    labels = paste0("Community ", seq_len(K)))
  sample_info$assignment_certainty <- apply(cp_mean, 1, max)

  fit <- structure(
    list(
      stan_fit        = stan_fit,
      sample_info     = sample_info,
      pi_mean         = pi_mean,
      beta_summary    = NULL,
      alpha_mean      = mean(rstan::extract(stan_fit, pars = "alpha")[[1]]),
      K               = K,
      N               = N,
      S               = S,
      taxa_names      = taxa_names,
      covariate_names = covariate_names,
      scale_info      = scale_info,
      counts          = counts,
      stan_data       = stan_data,
      call            = cl
    ),
    class = c("edna_dmm_fit", "list"))

  # ── Align labels, then rebuild every community-indexed quantity ─────────────
  fit <- dmm_align_labels(fit, method = method, verbose = verbose)

  # ── Diagnostics that survive alignment ──────────────────────────────────────
  dg  <- fit$alignment$diagnostics
  bad <- dg$rhat > rhat_threshold | dg$ess < ess_threshold
  if (any(bad)) {
    warning(sprintf(
      paste0("Convergence problems remain after label alignment (%s). ",
             "This is not label switching: either the chains have found ",
             "different partitions of the samples, or a chain has failed. ",
             "Inspect fit$alignment$diagnostics, and consider a smaller K."),
      paste(sprintf("%s Rhat %.3f, ESS %.0f", dg$quantity[bad],
                    dg$rhat[bad], dg$ess[bad]), collapse = "; ")),
      call. = FALSE)
  } else if (verbose) {
    message(sprintf("  after alignment: max Rhat %.3f, min ESS %.0f",
                    max(dg$rhat), min(dg$ess)))
  }

  fit
}


# =============================================================================
# S3 methods for edna_dmm_fit
# =============================================================================

#' @export
print.edna_dmm_fit <- function(x, ...) {
  cat("eDNA Dirichlet-Multinomial Mixture Model\n")
  cat("=========================================\n")
  cat(sprintf("  K (communities)  : %d\n", x$K))
  cat(sprintf("  N (samples)      : %d\n", x$N))
  cat(sprintf("  S (taxa)         : %d\n", x$S))
  cat(sprintf("  Covariates       : %s\n",
              if (length(x$covariate_names) > 0)
                paste(x$covariate_names, collapse = ", ")
              else "(none, intercept-only)"))
  cat(sprintf("  Mean alpha       : %.2f\n", x$alpha_mean))

  if (!is.null(x$alignment)) {
    cat(sprintf("  Chains           : %d, labels aligned by %s\n",
                x$alignment$chains, x$alignment$method))
    if (!is.na(x$alignment$min_agreement)) {
      cat(sprintf("  Chain agreement  : %.3f (lp Rhat %.3f)\n",
                  x$alignment$min_agreement, x$alignment$lp_rhat))
    }
  }

  cat("\nCommunity sizes (most probable assignment):\n")
  tbl <- table(x$sample_info$z_hat)
  for (nm in names(tbl)) {
    cat(sprintf("  %-15s: %d samples (%.1f%%)\n",
                nm, tbl[[nm]], 100 * tbl[[nm]] / x$N))
  }

  cat("\nUse summary() for convergence diagnostics, or pass this object to:\n")
  cat("  eDNA_dmm_structure()      structure bar plots\n")
  cat("  eDNA_dmm_compositions()   community composition plots\n")
  cat("  eDNA_dmm_nmds()           NMDS ordination\n")
  cat("  eDNA_dmm_beta_intervals() covariate coefficient plots\n")
  invisible(x)
}

#' @export
summary.edna_dmm_fit <- function(object, ...) {
  cat("eDNA DMM: fit summary\n")
  cat("=====================\n\n")

  cat("Call:\n  ")
  print(object$call)
  cat("\n")

  cat(sprintf("Dimensions: N=%d samples, S=%d taxa, K=%d communities\n\n",
              object$N, object$S, object$K))

  # Overdispersion
  alpha_draws <- rstan::extract(object$stan_fit, pars = "alpha")$alpha
  cat("Overdispersion (alpha):\n")
  cat(sprintf("  Mean = %.2f, Median = %.2f, 90%% CI = [%.2f, %.2f]\n",
              mean(alpha_draws), stats::median(alpha_draws),
              stats::quantile(alpha_draws, 0.05),
              stats::quantile(alpha_draws, 0.95)))
  cat("  (alpha much greater than 1: low overdispersion, near-multinomial;\n")
  cat("   alpha near 1: high overdispersion, typical of eDNA)\n\n")

  # Label alignment, and whether the chains agree once labels are common
  if (!is.null(object$alignment)) {
    al <- object$alignment
    cat("Label alignment:\n")
    cat(sprintf("  Method = %s, %.1f%% of draws relabelled\n",
                al$method, al$pct_permuted))
    cat(sprintf("  lp Rhat = %.3f (invariant to labelling)\n", al$lp_rhat))
    if (!is.na(al$min_agreement)) {
      cat(sprintf("  Chains agree on %.1f%% of sample assignments\n",
                  100 * al$min_agreement))
      cat("  (near 100%: one solution found repeatedly; well below that,\n")
      cat("   the chains found different partitions, which relabelling\n")
      cat("   cannot repair and a longer run will not fix)\n")
    }
    cat("\n  Diagnostics after alignment:\n")
    for (i in seq_len(nrow(al$diagnostics))) {
      r <- al$diagnostics[i, ]
      cat(sprintf("    %-6s Rhat = %.3f, ESS = %.0f  (worst: %s)\n",
                  r$quantity, r$rhat, r$ess, r$worst_parameter))
    }
    cat("\n")
  }

  # Assignment certainty
  cert <- object$sample_info$assignment_certainty
  cat("Assignment certainty across samples:\n")
  cat(sprintf("  Mean = %.2f, Min = %.2f, Max = %.2f\n",
              mean(cert), min(cert), max(cert)))
  cat(sprintf("  Decisive (>=80%% certainty): %d/%d samples (%.0f%%)\n\n",
              sum(cert >= 0.8), object$N, 100 * mean(cert >= 0.8)))

  # Covariate effects, on the centred scale
  if (!is.null(object$beta_summary) && nrow(object$beta_summary) > 0 &&
      length(object$covariate_names) > 0) {
    beta_show <- object$beta_summary[object$beta_summary$covariate != "intercept", ]
    if (nrow(beta_show) > 0) {
      cat("Covariate effects (deviation from the average community, 90% CI):\n")
      fmt <- "  %-14s | %-16s: %6.2f [%6.2f, %6.2f]  P(dir)=%3.0f%%  [%s]\n"
      for (i in seq_len(nrow(beta_show))) {
        r <- beta_show[i, ]
        p_dir <- max(r$prob_positive, r$prob_negative)
        cat(sprintf(fmt, r$community, r$covariate, r$mean, r$ci_5, r$ci_95,
                    100 * p_dir, r$reliability))
      }
      cat("\n  Reliability: trustworthy (ESS>400), cautious (100-400), unreliable (<100)\n\n")
    }
  }

  cat("HMC diagnostics:\n")
  rstan::check_hmc_diagnostics(object$stan_fit)

  invisible(object)
}
