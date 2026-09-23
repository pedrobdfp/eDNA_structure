# =============================================================================
# eDNA_loo(): Leave-one-out cross-validation across K
# =============================================================================

#' Compare DMM models across values of K
#'
#' @description
#' Fits a DMM for each value of `K` with [eDNA_dmm()] and compares them by
#' leave-one-out cross-validation (PSIS-LOO) via the \pkg{loo} package.
#' Because every fit runs several chains whose labels are aligned, the sweep
#' also reports whether the chains agreed on the same solution at each `K`,
#' which is a separate question from how well the model predicts.
#'
#' @section Reading the two answers separately:
#' LOO ELPD says how well a model predicts a held-out sample. It almost always
#' improves as `K` grows, because a more flexible mixture fits anything better,
#' so it rarely identifies a stopping point on its own.
#'
#' The alignment diagnostics say something different: whether the data pick out
#' a single grouping of the samples at that `K`. When several chains settle on
#' genuinely different partitions, no single partition can be reported, however
#' well the model predicts. `lp_rhat` and `min_agreement` in `loo_table` carry
#' that information, and `identified` combines them.
#'
#' @section Why `elpd` comes from one chain:
#' `elpd` is the score of the best single chain, not of the pooled posterior.
#' PSIS-LOO weights each draw by `1 / p(y_i | theta)`, so draws that predict a
#' sample worst carry the most weight. If the chains occupy different modes,
#' pooling drags the estimate toward the worse mode, far enough that the pooled
#' score can fall below every chain that went into it, which no predictive
#' score for a mixture of those posteriors can legitimately do. The pooled
#' value is still reported as `elpd_pooled`, and a warning is issued whenever
#' it falls below `elpd_worst_chain`, but it should not be compared across `K`.
#'
#' A paired comparison between two values of `K` must therefore be built from
#' one coherent posterior at each.
#'
#' @param counts,covariates,scale_covariates,conc,alpha_shape,alpha_rate
#'   Passed to [eDNA_dmm()].
#' @param K_range Integer vector of `K` to evaluate. Default `2:5`.
#' @param chains,cores,method,iter,warmup,adapt_delta,seed Passed to
#'   [eDNA_dmm()].
#' @param agreement_threshold Smallest acceptable share of samples on which two
#'   chains agree after alignment. Default `0.9`.
#' @param rhat_threshold Largest acceptable `lp__` Rhat. Default `1.1`.
#' @param save_dir Directory to write each fit to as `fit_K<k>.rds`. When a
#'   file is already present it is loaded instead of refitted, which makes a
#'   long sweep resumable. `NULL` (default) saves nothing.
#' @param keep_fits Keep every fit in the returned object. Default `TRUE`.
#'   Fits are large; with `save_dir` set, `FALSE` releases each one after
#'   scoring it and keeps memory flat across the sweep.
#' @param verbose Print progress. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{`loo_table`}{One row per `K`: `elpd` and `se` for the best chain,
#'     `best_chain`, `elpd_pooled`, `elpd_worst_chain`, `p_loo`, `pareto_bad`,
#'     `n_obs`, `n_divergent`, `pct_permuted`, `lp_rhat`, `min_agreement`,
#'     `max_agreement`, `pi_rhat`, `beta_rhat`, `identified`, `minutes`.}
#'   \item{`loo_by_chain`}{One row per chain per `K`: `elpd`, `se`,
#'     `pareto_bad`, `lp_mean`. Plugs straight into
#'     [eDNA_dmm_k_diagnostics()] as `elpd_by_run`.}
#'   \item{`loo_compare`}{Output of [loo::loo_compare()], or `NULL` for a
#'     single `K`.}
#'   \item{`plot`}{ELPD against `K`, with the individual chains shown and any
#'     `K` whose chains disagreed left unfilled.}
#'   \item{`fits`}{Named list of fits, or of `NULL` when `keep_fits = FALSE`.}
#' }
#'
#' @seealso [eDNA_dmm()] for a single fit, [eDNA_dmm_k_diagnostics()] for the
#'   wider set of diagnostics used to choose `K`.
#'
#' @examples
#' \dontrun{
#' d   <- get_example_data()
#' res <- eDNA_loo(d$counts, d$covariates, K_range = 2:6)
#'
#' res$loo_table[, c("K", "elpd", "lp_rhat", "min_agreement", "identified")]
#' res$plot
#' max(res$loo_table$K[res$loo_table$identified])   # largest identified K
#' }
#'
#' @export
eDNA_loo <- function(
    counts,
    covariates          = NULL,
    K_range             = 2:5,
    scale_covariates    = TRUE,
    chains              = 4,
    cores               = NULL,
    method              = c("STEPHENS", "ECR-pivot", "ECR-iterative"),
    iter                = 4000,
    warmup              = 2000,
    adapt_delta         = 0.95,
    seed                = 13,
    conc                = 0.5,
    alpha_shape         = 5,
    alpha_rate          = 2,
    agreement_threshold = 0.9,
    rhat_threshold      = 1.1,
    save_dir            = NULL,
    keep_fits           = TRUE,
    verbose             = TRUE
) {
  method <- match.arg(method)
  if (!all(K_range == round(K_range)) || any(K_range < 2))
    rlang::abort("`K_range` must be a vector of integers all >= 2.")
  K_range <- as.integer(K_range)
  if (!is.null(save_dir)) dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  fits       <- stats::setNames(vector("list", length(K_range)), paste0("K", K_range))
  loo_list   <- stats::setNames(vector("list", length(K_range)), paste0("K", K_range))
  rows       <- list()
  chain_rows <- list()

  for (i in seq_along(K_range)) {
    k  <- K_range[i]
    t0 <- Sys.time()
    if (verbose) message(sprintf("\n===== K = %d =====", k))

    # A finished fit on disk is reused, so a sweep interrupted partway through
    # resumes where it stopped rather than starting over.
    path <- if (is.null(save_dir)) NULL else
      file.path(save_dir, sprintf("fit_K%d.rds", k))
    if (!is.null(path) && file.exists(path)) {
      if (verbose) message("  reusing ", basename(path))
      fit <- readRDS(path)
      if (is.null(fit$alignment))
        fit <- dmm_align_labels(fit, method = method, verbose = verbose)
    } else {
      fit <- eDNA_dmm(
        counts = counts, covariates = covariates, K = k,
        scale_covariates = scale_covariates,
        chains = chains, cores = cores, method = method,
        iter = iter, warmup = warmup, adapt_delta = adapt_delta,
        seed = seed, conc = conc, alpha_shape = alpha_shape,
        alpha_rate = alpha_rate, verbose = verbose)
      if (!is.null(path)) saveRDS(fit, path)
    }

    # ── Predictive score: pooled, and each chain on its own ──────────────────
    ll_all <- loo::extract_log_lik(fit$stan_fit, parameter_name = "log_lik")
    lo     <- loo::loo(ll_all)
    loo_list[[i]] <- lo

    ll_ch <- loo::extract_log_lik(fit$stan_fit, parameter_name = "log_lik",
                                  merge_chains = FALSE)
    lp_ch <- rstan::extract(fit$stan_fit, pars = "lp__", permuted = FALSE)
    n_ch  <- dim(ll_ch)[2]
    ch_elpd <- ch_se <- numeric(n_ch); ch_bad <- integer(n_ch)
    for (ch in seq_len(n_ch)) {
      lc <- loo::loo(ll_ch[, ch, , drop = FALSE])
      ch_elpd[ch] <- lc$estimates["elpd_loo", "Estimate"]
      ch_se[ch]   <- lc$estimates["elpd_loo", "SE"]
      ch_bad[ch]  <- sum(loo::pareto_k_values(lc) > 0.7)
      chain_rows[[length(chain_rows) + 1L]] <- data.frame(
        K = k, chain = ch, elpd = ch_elpd[ch], se = ch_se[ch],
        pareto_bad = ch_bad[ch], lp_mean = mean(lp_ch[, ch, 1]),
        stringsAsFactors = FALSE)
    }

    elpd_pooled <- lo$estimates["elpd_loo", "Estimate"]
    best_ch     <- which.max(ch_elpd)
    if (elpd_pooled < min(ch_elpd))
      warning(sprintf(
        paste0("K = %d: the pooled ELPD (%.1f) is below every individual chain ",
               "(worst %.1f). The chains are in different modes, so the pooled ",
               "posterior is not a valid predictive score. Use `elpd` or ",
               "`loo_by_chain` rather than `elpd_pooled`."),
        k, elpd_pooled, min(ch_elpd)), call. = FALSE)

    # ── Did the chains find the same solution? ───────────────────────────────
    al <- fit$alignment
    dg <- al$diagnostics
    grab <- function(q, col) {
      v <- dg[[col]][dg$quantity == q]
      if (length(v)) v[1] else NA_real_
    }
    ag  <- al$chain_agreement
    off <- if (is.null(ag)) NA_real_ else ag[lower.tri(ag)]

    sp    <- rstan::get_sampler_params(fit$stan_fit, inc_warmup = FALSE)
    n_div <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))

    identified <- isTRUE(al$lp_rhat < rhat_threshold) &&
      (is.null(ag) || min(off) >= agreement_threshold)

    rows[[i]] <- data.frame(
      K                = k,
      elpd             = ch_elpd[best_ch],
      se               = ch_se[best_ch],
      best_chain       = best_ch,
      elpd_pooled      = elpd_pooled,
      elpd_worst_chain = min(ch_elpd),
      p_loo            = lo$estimates["p_loo", "Estimate"],
      pareto_bad       = ch_bad[best_ch],
      n_obs            = length(loo::pareto_k_values(lo)),
      n_divergent      = n_div,
      pct_permuted     = al$pct_permuted,
      lp_rhat          = al$lp_rhat,
      min_agreement    = if (is.null(ag)) NA_real_ else min(off),
      max_agreement    = if (is.null(ag)) NA_real_ else max(off),
      pi_rhat          = grab("pi", "rhat"),
      beta_rhat        = grab("beta", "rhat"),
      identified       = identified,
      minutes          = as.numeric(difftime(Sys.time(), t0, units = "mins")),
      stringsAsFactors = FALSE)

    if (verbose)
      message(sprintf(
        "K = %d: ELPD %.1f (SE %.1f) | lp Rhat %.3f | agreement %.2f | %s",
        k, rows[[i]]$elpd, rows[[i]]$se, rows[[i]]$lp_rhat,
        rows[[i]]$min_agreement,
        if (identified) "one solution" else "SEVERAL SOLUTIONS"))

    # Fits are large. With save_dir set they are already on disk, so releasing
    # them keeps memory flat however many K the sweep covers.
    if (keep_fits) fits[[i]] <- fit else { rm(fit); gc(verbose = FALSE) }
  }

  loo_df   <- do.call(rbind, rows)
  chain_df <- do.call(rbind, chain_rows)

  p <- ggplot2::ggplot(loo_df, ggplot2::aes(x = .data$K, y = .data$elpd)) +
    ggplot2::geom_point(data = chain_df, ggplot2::aes(x = .data$K, y = .data$elpd),
                        colour = "grey45", size = 1.5, alpha = 0.7,
                        inherit.aes = FALSE) +
    ggplot2::geom_line(linewidth = 0.9, colour = "#1D3557") +
    ggplot2::geom_point(ggplot2::aes(fill = .data$identified), shape = 21,
                        size = 3.4, stroke = 0.9, colour = "#1D3557") +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#1D3557", `FALSE` = "white"),
      labels = c(`TRUE` = "chains agree", `FALSE` = "several solutions"),
      name = NULL) +
    ggplot2::scale_x_continuous(breaks = K_range) +
    ggplot2::labs(x = "Number of communities (K)", y = "LOO ELPD") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")

  list(loo_table    = loo_df,
       loo_by_chain = chain_df,
       # loo_compare needs at least two models to compare.
       loo_compare  = if (length(loo_list) >= 2)
                        do.call(loo::loo_compare, unname(loo_list)) else NULL,
       plot         = p,
       fits         = fits)
}


# =============================================================================
# get_example_data(): built-in example dataset
# =============================================================================

#' Load the built-in example dataset
#'
#' @description
#' A small, deliberately plain example survey: a site by taxon count table plus
#' the two numeric covariates that separate the communities. Twenty sites,
#' 33 taxa, four true communities.
#'
#' Use it to learn the expected input format for [eDNA_dmm()], to follow the
#' vignettes, or to check that an installation works.
#'
#' @section Data format:
#' \describe{
#'   \item{`counts`}{Integer matrix, 20 sites by 33 taxa. Rows are sites
#'     (`STN_001` to `STN_020`), columns are taxa (`Sp_1`, `Sp_2`, and so on).
#'     This is the only required input to [eDNA_dmm()].}
#'   \item{`covariates`}{Data frame, 20 rows by 2 numeric columns: `Depth` in
#'     metres and `Distance_shore`. Every column is numeric, so it can be
#'     passed straight to [eDNA_dmm()] with nothing dropped.}
#'   \item{`metadata`}{Data frame, 20 rows: `sample_id`, `TrueCommunity`,
#'     `Depth`, `Distance_shore`. The labelling columns, kept separate from the
#'     model covariates. Pass this to the plotting functions' `metadata`
#'     argument.}
#' }
#'
#' `covariates` and `metadata` are separate on purpose. [eDNA_dmm()] requires
#' every covariate column to be numeric, so an identifier such as `sample_id`,
#' or a ground-truth column such as `TrueCommunity`, cannot be part of
#' `covariates`; the model would try to fit it. The plotting functions want
#' exactly those labelling columns. Keeping the two apart means both calls work
#' without the user having to subset anything.
#'
#' @section Community structure:
#' Four communities separated by two covariates:
#' \itemize{
#'   \item **Community 1**: deep (80 m) and offshore (200 units)
#'   \item **Community 2**: surface (10 m) and offshore (200 units)
#'   \item **Community 3**: deep (80 m) and inshore (20 units)
#'   \item **Community 4**: surface (10 m) and inshore (20 units)
#' }
#'
#' @section Formatting your own data:
#' The `counts` matrix format is the required input to [eDNA_dmm()]:
#' - **Rows** are samples, one row per site
#' - **Columns** are taxa or ASVs; taxonomic annotation is not required
#' - **Values** are non-negative integer read counts, not normalised
#'
#' If your count table is in long format (sample, taxon, count), pivot it:
#' ```r
#' library(tidyr)
#' count_matrix <- pivot_wider(
#'   long_data,
#'   names_from  = taxon,
#'   values_from = reads,
#'   values_fill = 0
#' ) |>
#'   tibble::column_to_rownames("sample_id") |>
#'   as.matrix()
#' ```
#'
#' @return A named list with elements `counts`, `covariates` and `metadata`.
#'
#' @seealso [simulate_eDNA_survey()] to generate your own.
#'
#' @examples
#' d <- get_example_data()
#'
#' dim(d$counts)          # 20 sites by 33 taxa
#' d$counts[1:3, 1:5]
#' head(d$covariates)     # numeric only
#'
#' \dontrun{
#' fit <- eDNA_dmm(counts = d$counts, covariates = d$covariates, K = 4)
#' print(fit)
#' eDNA_dmm_structure(fit, metadata = d$metadata, facet_var = "TrueCommunity")
#' }
#'
#' @export
get_example_data <- function() {
  # `example_edna` is a lazy-loaded dataset: it lives in the *attached*
  # package environment, not in the package namespace, so an unqualified
  # reference here would resolve through the search path and could be
  # shadowed by any object of the same name sitting in the caller's global
  # environment (e.g. a leftover from an old R session/workspace). Using
  # `::` forces the lookup to the package's own copy.
  eDNAstructure::example_edna
}
