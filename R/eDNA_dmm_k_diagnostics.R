# =============================================================================
# eDNA_dmm_k_diagnostics(): how many communities the data support
# =============================================================================


#' Diagnostics for choosing the number of communities
#'
#' @description
#' Summarises a set of fits across K and draws them as one figure. No single
#' panel picks K; together they usually make the choice obvious, and they fail
#' in different ways, which is the point of showing them side by side.
#'
#' `level = "simple"` draws the two predictive panels alone. `level =
#' "advanced"` (the default) adds the reliability, convergence and
#' interpretability panels beneath them.
#'
#' @section The panels:
#' **(a) Predictive fit.** LOO ELPD against K, averaged over runs, with the
#' spread across runs. Rising ELPD means the extra community earns its keep
#' predictively. This panel alone is often ambiguous: ELPD frequently keeps
#' creeping up long after the extra communities have stopped being
#' interpretable.
#'
#' **(b) Marginal gain.** The paired difference in pointwise ELPD at each step
#' in K, with twice its standard error. Because the models are fitted to the
#' same samples, sample-to-sample variation cancels in the difference, which
#' makes this a far sharper comparison than the ELPD totals in (a). Filled
#' points clear twice their standard error; open points do not.
#'
#' **(c) Reliability of (a) and (b).** Directly beneath the ELPD it qualifies:
#' the share of samples whose PSIS-LOO importance weights have an estimated
#' Pareto shape k above `pareto_threshold`. Above that threshold the
#' leave-one-out prediction for a sample rests on a handful of draws and its
#' ELPD contribution is not trustworthy, usually because the sample is
#' influential enough that dropping it would move the posterior. The share
#' grows with K, since each added community is supported by fewer samples, so
#' the approximation degrades fastest exactly where the ELPD curve is flattest.
#'
#' **(d) Convergence.** Whether the fits at each K found the same posterior at
#' all. With several chains per K this is the share of samples that the two
#' least-similar chains assign to the same community once labels have been
#' aligned ([dmm_align_labels()]); open points mark K where the
#' permutation-invariant `lp__` Rhat also exceeds `rhat_threshold`. Agreement
#' near one means one solution found repeatedly. A drop means the chains found
#' genuinely different partitions, which relabelling cannot reconcile and a
#' longer run will not fix -- the usual cause is a K the data do not identify.
#' Where chain agreement is unavailable the panel falls back to `lp__` Rhat.
#'
#' **(e) Assignment certainty.** The distribution across samples of the largest
#' membership probability, shown in full rather than as a count above a
#' threshold, which would hide where the mass sits and depend on an arbitrary
#' cut. When K exceeds what the data support, membership spreads thinly over
#' several communities instead of concentrating, and a "community" a sample
#' belongs to with probability 0.3 is not a useful description of that sample.
#' Plotted on the floor-corrected scale: raw certainty cannot fall below 1/K,
#' and that floor moves with K, so raw distributions are not comparable across
#' panels and part of any decline would be arithmetic. Rescaling puts 0 at
#' "spread evenly over every community" and 1 at "assigned with certainty" for
#' every K. The red diamond is the mean, reported as `mean_excess` in the table.
#'
#' **(f) Community distinctness.** The smallest distance between any two
#' community compositions. Compositions are compared with Aitchison distance by
#' default -- Euclidean distance between centred log-ratio transforms, the
#' standard metric for compositional data, which weighs a ratio between two rare
#' taxa as heavily as one between two abundant taxa. When an added community is
#' a near-copy of one already present, this collapses -- the model has begun
#' splitting hairs rather than finding structure. A sharp drop marks the K at
#' which that starts.
#'
#' The left column asks whether the model predicts well and whether that
#' judgement can be trusted; the right column asks whether the fit is a single
#' answer and whether that answer is interpretable. Panels (d) to (f) can all
#' deteriorate while (a) is still rising, and when they do, the higher K is
#' buying predictive accuracy at the cost of everything else.
#'
#' @param fits The return value of [eDNA_loo()] itself, the easiest way to
#'   call this function: `elpd_by_run`, `convergence`, and `fits` proper are
#'   all pulled out of it automatically (anything you pass explicitly for
#'   those still wins). Otherwise, a named list of `edna_dmm_fit` objects, or
#'   a function taking a single integer K and returning a fit (useful when
#'   fits are cached on disk and you would rather not hold them all in memory
#'   at once). The fit list itself is only needed for panels (e) and (f); pass
#'   `NULL` with `level = "simple"`.
#' @param K_values Integer vector of K to include. Required when `fits` is a
#'   function or `NULL`; inferred from `fits` otherwise.
#' @param level `"advanced"` (default) draws six panels; `"simple"` draws the
#'   two predictive panels, (a) and (b), alone.
#' @param elpd_by_run Optional data frame with columns `K`, `elpd` and a run
#'   identifier (`seed` or `chain`), one row per run, used for panels (a) and
#'   (b). When `NULL` (default), LOO is computed from each fit, which is slower.
#' @param adjacent Optional data frame with columns `K_from`, `K_to`, `gain`,
#'   `se_gain` and `worth_it` for panel (b). When `NULL`, panel (b) is drawn
#'   from `elpd_by_run` differences without the paired standard error, which is
#'   a weaker comparison; supply the paired table where you have it.
#' @param pareto_by_run Optional data frame with one row per run and columns
#'   `K` and either `pct_bad` (percentage of samples above the threshold) or
#'   `n_pareto_bad` and `n_obs`, for panel (c). When `NULL` and `elpd_by_run`
#'   is also `NULL`, it is computed alongside LOO; when `NULL` and
#'   `elpd_by_run` is supplied, panel (c) is left blank.
#' @param convergence Optional data frame with column `K` and at least one of
#'   `min_agreement` and `lp_rhat`, for panel (d); one row per K, as returned
#'   in the `loo_table` of [eDNA_loo()]. When `NULL`, panel (d) is left
#'   blank.
#' @param pareto_threshold Pareto shape k above which a sample's LOO
#'   contribution is treated as unreliable. Default `0.7`.
#' @param rhat_threshold `lp__` Rhat above which a K is marked unconverged in
#'   panel (d). Default `1.1`, the classical Gelman-Rubin rule of thumb, rather
#'   than the stricter `1.01` now recommended for final inference: the panel's
#'   substantive signal is the agreement on the vertical axis, which separates
#'   cleanly, and the fill is a secondary flag. Be aware that between roughly
#'   `1.01` and `1.1` this choice can decide a borderline K on its own; where it
#'   does, read the chain agreement rather than the fill.
#' @param agreement_threshold Chain agreement drawn as a reference line in
#'   panel (d). Default `0.9`.
#' @param certainty_threshold Membership probability above which a sample
#'   counts as confidently assigned. Reported as `pct_confident` in the table;
#'   panel (e) shows the whole distribution and does not use it. Default `0.8`.
#' @param distance Distance between community compositions in panel (f):
#'   `"aitchison"` (default) or `"tv"` for total variation. Aitchison respects
#'   the simplex geometry; total variation is bounded on `[0, 1]` and easier to
#'   read but is dominated by the abundant taxa.
#' @param base_size Base font size. Default `12`.
#'
#' @return A list with:
#' \describe{
#'   \item{`plot`}{The assembled figure.}
#'   \item{`panels`}{The panels separately, named `a` to `f`. Panels that could
#'     not be drawn are `NULL`.}
#'   \item{`table`}{One row per K: `elpd_mean`, `elpd_sd` and, where the inputs
#'     allow, `pct_pareto_bad`, `min_agreement`, `lp_rhat`, `min_distance`,
#'     `mean_certainty`, `median_certainty`, `q10_certainty`, `mean_excess`,
#'     `pct_confident`.}
#' }
#'
#' @seealso [eDNA_loo()], which fits across K and produces the `elpd_by_run`,
#'   `pareto_by_run` and `convergence` inputs used here.
#'
#' @examples
#' \dontrun{
#' # The easy way: hand it eDNA_loo()'s return value directly
#' res <- eDNA_loo(counts, covariates, K_range = 2:6)
#' eDNA_dmm_k_diagnostics(res)$plot
#'
#' # The two predictive panels only, from a table of runs alone
#' eDNA_dmm_k_diagnostics(NULL, K_values = 2:10, level = "simple",
#'                        elpd_by_run = runs)$plot
#'
#' # Everything, from fits cached on disk instead of held in memory
#' d <- eDNA_dmm_k_diagnostics(
#'   fits        = function(k) readRDS(sprintf("fits/fit_K%d.rds", k)),
#'   K_values    = 2:10,
#'   elpd_by_run = runs,
#'   convergence = res$loo_table)
#' d$plot
#' }
#'
#' @export
eDNA_dmm_k_diagnostics <- function(
    fits,
    K_values            = NULL,
    level               = c("advanced", "simple"),
    elpd_by_run         = NULL,
    adjacent            = NULL,
    pareto_by_run       = NULL,
    convergence         = NULL,
    pareto_threshold    = 0.7,
    rhat_threshold      = 1.1,
    agreement_threshold = 0.9,
    certainty_threshold = 0.8,
    distance            = c("aitchison", "tv"),
    base_size           = 12
) {
  level    <- match.arg(level)
  distance <- match.arg(distance)
  advanced <- level == "advanced"

  # Accept eDNA_loo()'s return value directly as `fits`, so
  # eDNA_dmm_k_diagnostics(loo_result) is enough on its own. Detected by the
  # specific combination eDNA_loo() actually returns (a $fits list alongside
  # $loo_table and $loo_by_chain, and no $stan_fit at the top level, which
  # rules out someone having passed a single edna_dmm_fit by mistake).
  # Anything explicitly passed for elpd_by_run/convergence still wins, so this
  # never overrides a caller who wants to supply their own.
  if (is.list(fits) && !is.function(fits) && is.null(fits$stan_fit) &&
      !is.null(fits$fits) && !is.null(fits$loo_table)) {
    loo_result <- fits
    if (is.null(elpd_by_run)) elpd_by_run <- loo_result$loo_by_chain
    if (is.null(convergence)) convergence <- loo_result$loo_table
    fits <- loo_result$fits
  }

  # eDNA_loo()'s `fits`, the documented source for this argument, names its
  # elements "K2", "K3", ... (see eDNA_loo()'s return value). Try that
  # convention first and fall back to a plain integer name, so a list named
  # either way (or K_values passed as strings) resolves correctly instead of
  # silently returning NULL and leaving panels (e)/(f) empty.
  get_fit <- if (is.function(fits)) fits else function(k) {
    f <- fits[[paste0("K", k)]]
    if (is.null(f)) f <- fits[[as.character(k)]]
    f
  }
  if (is.null(K_values)) {
    if (is.null(fits) || is.function(fits))
      rlang::abort("`K_values` is required unless `fits` is a named list.")
    K_values <- sort(as.integer(sub("^K", "", names(fits))))
  }
  K_values <- sort(as.integer(K_values))

  # -- Per-K structural summaries, for panels (e) and (f) only -----------------
  tbl     <- data.frame(K = K_values)
  cert_df <- NULL
  if (advanced && !is.null(fits)) {
    rows <- list(); cert_all <- list()
    for (k in K_values) {
      f <- get_fit(k)
      if (is.null(f)) next

      d <- .community_distance(f$pi_mean, distance)
      diag(d) <- NA_real_

      cert <- f$sample_info$assignment_certainty
      cert_all[[as.character(k)]] <- data.frame(
        K = k, certainty = cert,
        # Rescaled between the 1/K floor and one, so distributions from
        # different K are on the same axis.
        excess = (cert - 1 / k) / (1 - 1 / k))
      rows[[length(rows) + 1L]] <- data.frame(
        K                = k,
        min_distance     = min(d, na.rm = TRUE),
        mean_certainty   = mean(cert),
        median_certainty = stats::median(cert),
        q10_certainty    = stats::quantile(cert, 0.10, names = FALSE),
        # Certainty cannot fall below 1/K, and that floor moves as K grows, so
        # a raw certainty is not comparable across K. The excess above the
        # floor, rescaled to [0, 1], is: 0 for a sample spread evenly over all
        # communities, 1 for a sample assigned with probability one.
        mean_excess      = mean((cert - 1 / k) / (1 - 1 / k)),
        pct_confident    = 100 * mean(cert >= certainty_threshold))
    }
    if (length(rows)) {
      tbl     <- merge(tbl, do.call(rbind, rows), by = "K", all.x = TRUE)
      cert_df <- do.call(rbind, cert_all)
      cert_df$Kf <- factor(cert_df$K, levels = K_values)
    }
  }

  # -- ELPD, and the reliability of the approximation behind it ----------------
  if (is.null(elpd_by_run)) {
    if (is.null(fits))
      rlang::abort("Supply `elpd_by_run`, or `fits` to compute it from.")
    er <- list(); pkr <- list()
    for (k in K_values) {
      f <- get_fit(k); if (is.null(f)) next
      ll <- loo::extract_log_lik(f$stan_fit, parameter_name = "log_lik")
      lo <- loo::loo(ll)
      er[[length(er) + 1L]] <- data.frame(
        K = k, elpd = sum(lo$pointwise[, "elpd_loo"]))
      pkr[[length(pkr) + 1L]] <- data.frame(
        K = k, pct_bad = 100 * mean(loo::pareto_k_values(lo) > pareto_threshold))
    }
    elpd_by_run <- do.call(rbind, er)
    if (is.null(pareto_by_run)) pareto_by_run <- do.call(rbind, pkr)
  }
  elpd_by_run <- elpd_by_run[elpd_by_run$K %in% K_values, , drop = FALSE]
  summ <- stats::aggregate(elpd ~ K, elpd_by_run,
                           function(v) c(m = mean(v), s = stats::sd(v)))
  summ <- data.frame(K = summ$K, elpd_mean = summ$elpd[, "m"],
                     elpd_sd = summ$elpd[, "s"])
  summ$elpd_sd[is.na(summ$elpd_sd)] <- 0
  tbl <- merge(tbl, summ, by = "K", all = TRUE)

  # If not supplied, try building it from what elpd_by_run and convergence
  # already carry: eDNA_loo()'s loo_by_chain has a per-chain `pareto_bad`
  # count, and its loo_table has `n_obs`, so the common eDNA_loo() ->
  # eDNA_dmm_k_diagnostics() pipeline gets panel (c) for free without the
  # caller having to assemble a separate table by hand.
  if (advanced && is.null(pareto_by_run) &&
      !is.null(elpd_by_run$pareto_bad) && !is.null(convergence$n_obs)) {
    pareto_by_run <- merge(
      elpd_by_run[, intersect(c("K", "pareto_bad"), names(elpd_by_run)),
                 drop = FALSE],
      unique(convergence[, c("K", "n_obs")]),
      by = "K")
    names(pareto_by_run)[names(pareto_by_run) == "pareto_bad"] <- "n_pareto_bad"
  }

  # Either a ready-made percentage or the counts it is built from is accepted,
  # since LOO tables store one or the other depending on how they were written.
  pk <- pr <- NULL
  if (advanced && !is.null(pareto_by_run)) {
    pr <- pareto_by_run[pareto_by_run$K %in% K_values, , drop = FALSE]
    if (is.null(pr$pct_bad)) {
      n_bad <- if (!is.null(pr$n_pareto_bad)) pr$n_pareto_bad else pr$n_high_k
      if (is.null(n_bad) || is.null(pr$n_obs))
        rlang::abort(c(
          "`pareto_by_run` needs a `pct_bad` column, or `n_pareto_bad` and `n_obs`.",
          i = "One row per run, with the K it was fitted at."))
      pr$pct_bad <- 100 * n_bad / pr$n_obs
    }
    pk <- stats::aggregate(pct_bad ~ K, pr, mean)
    names(pk)[2] <- "pct_pareto_bad"
    tbl <- merge(tbl, pk, by = "K", all.x = TRUE)
  }

  # Chain agreement is the readable quantity and lp__ Rhat the rigorous one;
  # both are kept, and the panel prefers agreement when it is available.
  cv <- NULL
  if (advanced && !is.null(convergence)) {
    cv   <- convergence[convergence$K %in% K_values, , drop = FALSE]
    keep <- intersect(c("K", "min_agreement", "lp_rhat"), names(cv))
    if (length(keep) < 2)
      rlang::abort("`convergence` needs a `K` column and `min_agreement` or `lp_rhat`.")
    cv  <- cv[, keep, drop = FALSE]
    tbl <- merge(tbl, cv, by = "K", all.x = TRUE)
  }

  blue <- "#1D3557"; red <- "#C1666B"
  base <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  # (a) predictive fit
  pa <- ggplot2::ggplot(summ, ggplot2::aes(x = .data$K, y = .data$elpd_mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$elpd_mean - .data$elpd_sd,
                                      ymax = .data$elpd_mean + .data$elpd_sd),
                         fill = "grey80", alpha = 0.7) +
    ggplot2::geom_point(data = elpd_by_run, ggplot2::aes(y = .data$elpd),
                        colour = "grey45", size = 1.5, alpha = 0.7) +
    ggplot2::geom_line(linewidth = 0.9, colour = blue) +
    ggplot2::geom_point(size = 3, colour = blue) +
    ggplot2::scale_x_continuous(breaks = K_values) +
    ggplot2::labs(x = "Number of communities (K)", y = "LOO ELPD") + base

  # (b) marginal gain. Filled versus open carries the 2 SE verdict; no legend,
  # so the distinction belongs in the figure caption.
  #
  # When possible this is computed properly: for each adjacent (K, K+1) pair,
  # loo::loo_compare() on that pair's own best chain (the same chain
  # elpd_by_run/panel (a) already treat as the valid predictive score, see
  # eDNA_loo()). That gives a real paired SE, unlike differencing the
  # already-averaged per-K means in `summ`, which throws the pairing away.
  # This needs one stan_fit per K (i.e. `fits` as a real list, not an
  # on-disk-loading function) and a `chain` column in `elpd_by_run` to know
  # which chain is "best" at each K; without both, it falls back to an
  # unpaired difference with no SE, exactly as before.
  if (is.null(adjacent)) {
    can_pair <- advanced && !is.null(fits) && !is.function(fits) &&
      !is.null(elpd_by_run$chain) && requireNamespace("loo", quietly = TRUE)
    paired <- NULL
    if (can_pair) {
      paired <- tryCatch({
        best_chain_at <- function(k) {
          sub <- elpd_by_run[elpd_by_run$K == k, , drop = FALSE]
          sub$chain[which.max(sub$elpd)]
        }
        rows <- list()
        for (i in seq_len(length(K_values) - 1)) {
          k1 <- K_values[i]; k2 <- K_values[i + 1]
          f1 <- get_fit(k1); f2 <- get_fit(k2)
          if (is.null(f1) || is.null(f2)) next
          ll1 <- loo::extract_log_lik(f1$stan_fit, "log_lik",
                                      merge_chains = FALSE)[, best_chain_at(k1), ]
          ll2 <- loo::extract_log_lik(f2$stan_fit, "log_lik",
                                      merge_chains = FALSE)[, best_chain_at(k2), ]
          cmp <- loo::loo_compare(list(.k1 = loo::loo(ll1), .k2 = loo::loo(ll2)))
          cmp <- as.data.frame(cmp)
          r1 <- cmp[cmp$model == ".k1", , drop = FALSE]
          r2 <- cmp[cmp$model == ".k2", , drop = FALSE]
          # loo_compare() zeroes out whichever model is better; read the gain
          # (elpd at k2 minus elpd at k1) off whichever row is NOT the zero.
          if (isTRUE(all.equal(r1$elpd_diff, 0))) {
            gain <- r2$elpd_diff; se <- r2$se_diff
          } else {
            gain <- -r1$elpd_diff; se <- r1$se_diff
          }
          rows[[length(rows) + 1L]] <- data.frame(
            K_from = k1, K_to = k2, gain = gain, se_gain = se,
            worth_it = (gain - 2 * se) > 0)
        }
        if (length(rows)) do.call(rbind, rows) else NULL
      }, error = function(e) NULL)
    }
    adjacent <- if (!is.null(paired)) paired else data.frame(
      K_from = utils::head(summ$K, -1), K_to = summ$K[-1],
      gain = diff(summ$elpd_mean), se_gain = NA_real_, worth_it = NA)
  }
  # Panel (b) is only readable as the increments of panel (a) if it was
  # computed from the same posteriors. It is easy to break that without
  # noticing: supply paired gains from one fit per K while panel (a) averages
  # several runs per K, and the two panels disagree, in places by enough to
  # flip a step's sign. Worse, differencing pooled multi-chain posteriors is
  # not merely inconsistent but invalid, because PSIS-LOO on a multimodal
  # posterior is dominated by its worst mode. Check and say so.
  implied <- summ$elpd_mean[match(adjacent$K_to, summ$K)] -
             summ$elpd_mean[match(adjacent$K_from, summ$K)]
  off_by  <- abs(implied - adjacent$gain)
  tol     <- pmax(2 * ifelse(is.na(adjacent$se_gain), 0, adjacent$se_gain), 1)
  if (any(!is.na(off_by) & off_by > tol)) {
    i <- which(!is.na(off_by) & off_by > tol)
    warning(sprintf(
      paste0("Panel (b) is not the increments of panel (a). Steps %s differ by ",
             "%s. `adjacent` and `elpd_by_run` must describe the same posterior ",
             "at each K, and neither may come from pooling chains that sit in ",
             "different modes."),
      paste(sprintf("%d->%d", adjacent$K_from[i], adjacent$K_to[i]), collapse = ", "),
      paste(sprintf("%.1f", off_by[i]), collapse = ", ")), call. = FALSE)
  }

  adjacent$step <- factor(sprintf("%d->%d", adjacent$K_from, adjacent$K_to),
                          levels = sprintf("%d->%d", adjacent$K_from, adjacent$K_to))
  pb <- ggplot2::ggplot(adjacent, ggplot2::aes(x = .data$step, y = .data$gain)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    { if (any(!is.na(adjacent$se_gain)))
        ggplot2::geom_linerange(
          ggplot2::aes(ymin = .data$gain - 2 * .data$se_gain,
                       ymax = .data$gain + 2 * .data$se_gain),
          linewidth = 0.8, colour = "grey30") else NULL } +
    ggplot2::geom_point(ggplot2::aes(fill = .data$worth_it), shape = 21,
                        size = 3.2, stroke = 0.9, colour = "grey20") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = blue, `FALSE` = "white"),
                               na.value = "grey70", guide = "none") +
    ggplot2::labs(x = "Step in K", y = "Paired gain in ELPD") + base +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # (c) reliability of the LOO approximation, under the ELPD it qualifies
  pc <- NULL
  if (!is.null(pk)) {
    pc <- ggplot2::ggplot(pk, ggplot2::aes(x = .data$K,
                                           y = .data$pct_pareto_bad)) +
      ggplot2::geom_point(data = pr, ggplot2::aes(x = .data$K, y = .data$pct_bad),
                          colour = "grey45", size = 1.5, alpha = 0.7,
                          inherit.aes = FALSE) +
      ggplot2::geom_line(linewidth = 0.9, colour = red) +
      ggplot2::geom_point(size = 3, colour = red) +
      ggplot2::scale_x_continuous(breaks = K_values) +
      ggplot2::expand_limits(y = 0) +
      ggplot2::labs(x = "Number of communities (K)",
                    y = sprintf("%% of samples with\nPareto k > %.1f",
                                pareto_threshold)) +
      base
  }

  # (d) convergence: one posterior found repeatedly, or several found once each?
  pd <- NULL
  if (!is.null(cv)) {
    has_agr <- !is.null(cv$min_agreement) && any(!is.na(cv$min_agreement))
    if (has_agr) {
      cv$converged <- if (is.null(cv$lp_rhat)) TRUE else
        !is.na(cv$lp_rhat) & cv$lp_rhat < rhat_threshold
      pd <- ggplot2::ggplot(cv[!is.na(cv$min_agreement), , drop = FALSE],
                            ggplot2::aes(x = .data$K, y = .data$min_agreement)) +
        ggplot2::geom_hline(yintercept = agreement_threshold, linetype = "dashed",
                            colour = "grey40") +
        ggplot2::geom_line(linewidth = 0.9, colour = blue) +
        ggplot2::geom_point(ggplot2::aes(fill = .data$converged), shape = 21,
                            size = 3.2, stroke = 0.9, colour = blue) +
        ggplot2::scale_fill_manual(values = c(`TRUE` = blue, `FALSE` = "white"),
                                   na.value = "grey70", guide = "none") +
        ggplot2::scale_y_continuous(limits = c(0, 1),
                                    labels = scales::percent_format(accuracy = 1)) +
        ggplot2::scale_x_continuous(breaks = K_values) +
        ggplot2::expand_limits(x = range(K_values)) +
        ggplot2::labs(x = "Number of communities (K)",
                      y = "Samples the least-similar\nchains assign alike") +
        base
    } else {
      pd <- ggplot2::ggplot(cv[!is.na(cv$lp_rhat), , drop = FALSE],
                            ggplot2::aes(x = .data$K, y = .data$lp_rhat)) +
        ggplot2::geom_hline(yintercept = rhat_threshold, linetype = "dashed",
                            colour = "grey40") +
        ggplot2::geom_line(linewidth = 0.9, colour = blue) +
        ggplot2::geom_point(size = 3, colour = blue) +
        ggplot2::scale_x_continuous(breaks = K_values) +
        ggplot2::expand_limits(x = range(K_values)) +
        ggplot2::labs(x = "Number of communities (K)",
                      y = "R-hat of the log posterior density") +
        base
    }
  }

  # (e) assignability, as a distribution rather than a count above a threshold:
  # a threshold hides where the mass sits and depends on an arbitrary cut.
  # Plotted on the floor-corrected scale so panels are comparable across K.
  pe <- NULL
  if (!is.null(cert_df)) {
    cert_mean <- stats::aggregate(excess ~ Kf, cert_df, mean)
    pe <- ggplot2::ggplot(cert_df, ggplot2::aes(x = .data$Kf, y = .data$excess)) +
      ggplot2::geom_violin(fill = "grey85", colour = NA, scale = "width",
                           width = 0.9, trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.16, outlier.size = 0.5, outlier.alpha = 0.4,
                            fill = "white", colour = blue, linewidth = 0.4) +
      ggplot2::geom_point(data = cert_mean,
                          ggplot2::aes(x = .data$Kf, y = .data$excess),
                          shape = 23, size = 2.2, fill = red, colour = red,
                          inherit.aes = FALSE) +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(x = "Number of communities (K)",
                    y = "Assignment certainty") +
      base
  }

  # (f) distinctness
  pf <- NULL
  if (!is.null(tbl$min_distance) && any(!is.na(tbl$min_distance))) {
    pf <- ggplot2::ggplot(tbl[!is.na(tbl$min_distance), , drop = FALSE],
                          ggplot2::aes(x = .data$K, y = .data$min_distance)) +
      ggplot2::geom_line(linewidth = 0.9, colour = blue) +
      ggplot2::geom_point(size = 3, colour = blue) +
      ggplot2::scale_x_continuous(breaks = K_values) +
      ggplot2::expand_limits(y = 0) +
      ggplot2::labs(
        x = "Number of communities (K)",
        y = sprintf("%s distance between\nthe closest two communities",
                    if (distance == "aitchison") "Aitchison" else "Total variation")) +
      base
  }

  # A panel that could not be drawn keeps its slot rather than reflowing the
  # grid, so a letter always refers to the same quantity whatever was supplied.
  p <- if (advanced)
    cowplot::plot_grid(pa, pb, pc, pd, pe, pf, ncol = 2,
                       labels = c("a", "b", "c", "d", "e", "f"), label_size = 15)
  else
    cowplot::plot_grid(pa, pb, ncol = 2, labels = c("a", "b"), label_size = 15)

  list(plot = p,
       panels = list(a = pa, b = pb, c = pc, d = pd, e = pe, f = pf),
       table = tbl)
}


# Pairwise distance between community compositions.
#
#   "aitchison"  Euclidean distance between centred log-ratio transforms. The
#                standard metric for compositional data: scale-invariant and
#                subcompositionally coherent, and it weighs a ratio between two
#                rare taxa as heavily as one between two abundant taxa. That is
#                usually what "are these two communities the same?" means, since
#                a community defined by its rare taxa is still a distinct
#                community.
#
#   "tv"         Total variation, half the L1 distance. Bounded on [0, 1] and
#                easy to read, but dominated by the abundant taxa: two
#                communities agreeing on their dominants look close however much
#                their rare taxa differ.
#
# Aitchison is undefined at zero. Posterior mean compositions are strictly
# positive here, so no flooring is applied, but note that very small entries
# carry large log-ratios and can dominate the distance.
#' @keywords internal
.community_distance <- function(pm, distance = c("aitchison", "tv")) {
  distance <- match.arg(distance)
  if (distance == "tv") return(as.matrix(stats::dist(pm, method = "manhattan")) / 2)
  if (any(pm <= 0))
    rlang::abort(c(
      "Aitchison distance is undefined for zero or negative proportions.",
      i = 'Use distance = "tv", or floor the compositions before calling.'))
  lg  <- log(pm)
  clr <- lg - rowMeans(lg)
  as.matrix(stats::dist(clr, method = "euclidean"))
}
