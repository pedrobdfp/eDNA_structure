# =============================================================================
# eDNA_loo(): Leave-one-out cross-validation for K selection
# =============================================================================

#' Compare DMM models across values of K using LOO cross-validation
#'
#' @description
#' Fits DMM models for a range of K values and compares them using
#' Leave-One-Out cross-validation (LOO-CV) via the `loo` package. Returns
#' both a comparison table and an elbow plot of LOO-ELPD vs K.
#'
#' LOO-ELPD (expected log predictive density) measures how well each model
#' predicts held-out observations. Higher values are better. The "elbow", 
#' the point where additional communities yield diminishing LOO-ELPD gains, 
#' is a useful heuristic for choosing K.
#'
#' @section Replicated data:
#' `eDNA_loo()` always uses the standard (summed) model; it has no `replication`
#' argument. This is deliberate. The replicate-aware model in [eDNA_dmm()] gives
#' each station its own composition parameter `theta_i`, so holding a station out
#' leaves that station's own parameter unidentified and station-level LOO is not
#' well defined. It is also considerably slower, which matters when sweeping many
#' values of K.
#'
#' The recommended workflow with replicated data is therefore: select K here on
#' the summed model, then refit at the chosen K with `replication` supplied to
#' [eDNA_dmm()] for the final parameter estimates. If you have replicates, sum
#' them per station for this step (e.g. `rowsum(rep_counts, replication)`).
#'
#' @param counts A count matrix (same as passed to [eDNA_dmm()]).
#' @param covariates A covariate matrix or data frame, or `NULL`.
#' @param K_range An integer vector of K values to evaluate. Default is `2:5`.
#'   Example: `K_range = 2:7`.
#' @param scale_covariates Logical. Default `TRUE`. See [eDNA_dmm()].
#' @param chains Integer. Number of chains per fit. Default `1`.
#' @param iter Integer. Iterations per chain. Default `4000`.
#' @param warmup Integer. Warmup iterations. Default `2000`.
#' @param adapt_delta Number. Default `0.95`.
#' @param seed Integer. Default `13`.
#' @param conc Number. Default `0.5`. See [eDNA_dmm()].
#' @param alpha_shape Number. Default `5`. See [eDNA_dmm()].
#' @param alpha_rate Number. Default `2`. See [eDNA_dmm()].
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{`loo_table`}{A data frame with LOO-ELPD estimates and SE per K.}
#'   \item{`loo_compare`}{The output of [loo::loo_compare()].}
#'   \item{`plot`}{A [ggplot2::ggplot()] elbow plot of LOO-ELPD vs K.}
#'   \item{`fits`}{A named list of `edna_dmm_fit` objects, one per K.}
#' }
#'
#' @seealso [eDNA_dmm()]
#' @export
eDNA_loo <- function(
    counts,
    covariates      = NULL,
    K_range         = 2:5,
    scale_covariates = TRUE,
    chains          = 1,
    iter            = 4000,
    warmup          = 2000,
    adapt_delta     = 0.95,
    seed            = 13,
    conc            = 0.5,
    alpha_shape     = 5,
    alpha_rate      = 2,
    verbose         = TRUE
) {
  if (!all(K_range == round(K_range)) || any(K_range < 2)) {
    rlang::abort("`K_range` must be a vector of integers all >= 2.")
  }
  K_range <- as.integer(K_range)

  fits      <- vector("list", length(K_range))
  names(fits) <- paste0("K", K_range)
  loo_list  <- vector("list", length(K_range))

  for (i in seq_along(K_range)) {
    k <- K_range[i]
    if (verbose) message(sprintf("\n===== Fitting K = %d =====", k))

    fits[[i]] <- eDNA_dmm(
      counts           = counts,
      covariates       = covariates,
      K                = k,
      scale_covariates = scale_covariates,
      chains           = chains,
      iter             = iter,
      warmup           = warmup,
      adapt_delta      = adapt_delta,
      seed             = seed,
      conc             = conc,
      alpha_shape      = alpha_shape,
      alpha_rate       = alpha_rate,
      verbose          = verbose
    )
    ll_mat       <- loo::extract_log_lik(fits[[i]]$stan_fit, parameter_name = "log_lik")
    loo_list[[i]] <- loo::loo(ll_mat)
    if (verbose) message(sprintf("K = %d: LOO-ELPD = %.1f (SE = %.1f)",
                                 k,
                                 loo_list[[i]]$estimates["elpd_loo", "Estimate"],
                                 loo_list[[i]]$estimates["elpd_loo", "SE"]))
  }

  names(loo_list) <- paste0("K", K_range)

  # ── LOO comparison ────────────────────────────────────────────────────────────
  loo_compare_result <- do.call(loo::loo_compare, unname(loo_list))

  loo_df <- data.frame(
    K    = K_range,
    elpd = vapply(loo_list, function(l) l$estimates["elpd_loo", "Estimate"], numeric(1)),
    se   = vapply(loo_list, function(l) l$estimates["elpd_loo", "SE"],       numeric(1))
  )

  # ── Elbow plot ────────────────────────────────────────────────────────────────
  p_loo <- ggplot2::ggplot(loo_df, ggplot2::aes(x = .data$K, y = .data$elpd)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$elpd - .data$se, ymax = .data$elpd + .data$se),
      alpha = 0.2, fill = "#457B9D"
    ) +
    ggplot2::geom_line(linewidth = 0.9, color = "#457B9D") +
    ggplot2::geom_point(size = 3.5, color = "#457B9D") +
    ggplot2::scale_x_continuous(breaks = K_range) +
    ggplot2::labs(
      x        = "K (number of communities)",
      y        = "LOO-ELPD (higher is better)",
      title    = "K selection - LOO cross-validation",
      subtitle = "Shaded band = \u00b11 SE  |  Elbow = point of diminishing returns"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey40")
    )

  list(
    loo_table   = loo_df,
    loo_compare = loo_compare_result,
    plot        = p_loo,
    fits        = fits
  )
}


# =============================================================================
# get_example_data() / get_example_replicates(): Built-in example datasets
# =============================================================================

#' Load the built-in example dataset
#'
#' @description
#' A small, deliberately plain example survey: a site-by-taxon count table plus
#' the two numeric covariates that separate the communities. Twenty sites,
#' 33 taxa, four true communities.
#'
#' Use it to learn the expected input format for [eDNA_dmm()], to follow the
#' vignettes, or to check that your installation works.
#'
#' @section Data format:
#' \describe{
#'   \item{`counts`}{Integer matrix, 20 sites x 33 taxa. Rows are sites
#'     (`STN_001`-`STN_020`), columns are taxa (`Sp_1`, `Sp_2`, ...).
#'     This is the only required input to [eDNA_dmm()].}
#'   \item{`covariates`}{Data frame, 20 rows x 2 numeric columns: `Depth` (m)
#'     and `Distance_shore`. Every column is numeric, so this can be passed
#'     straight to [eDNA_dmm()] with nothing dropped.}
#'   \item{`metadata`}{Data frame, 20 rows: `sample_id`, `TrueCommunity`,
#'     `Depth`, `Distance_shore`. The labelling columns, kept separate from
#'     the model covariates. Pass this to the plotting functions'
#'     `metadata` argument.}
#' }
#'
#' `covariates` and `metadata` are separate on purpose. `eDNA_dmm()` requires
#' every covariate column to be numeric, so an identifier column such as
#' `sample_id`, or a ground-truth column such as `TrueCommunity`, cannot be
#' part of `covariates` - the model would try to fit it. The plotting
#' functions want exactly those labelling columns. Keeping the two apart means
#' both calls work without the user having to subset anything.
#'
#' @section Community structure:
#' Four communities separated by two covariates:
#' \itemize{
#'   \item **Community 1**: deep (80 m) + offshore (200 units)
#'   \item **Community 2**: surface (10 m) + offshore (200 units)
#'   \item **Community 3**: deep (80 m) + inshore (20 units)
#'   \item **Community 4**: surface (10 m) + inshore (20 units)
#' }
#'
#' @section Formatting your own data:
#' The `counts` matrix format is the required input to [eDNA_dmm()]:
#' - **Rows** = samples (one row per site, or per replicate - see
#'   [get_example_replicates()])
#' - **Columns** = taxa, or ASVs; taxonomic annotation is not required
#' - **Values** = non-negative integer read counts, not normalised
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
#' @seealso [get_example_replicates()] for the replicated version of the same
#'   survey, [simulate_eDNA_survey()] to generate your own.
#'
#' @examples
#' d <- get_example_data()
#'
#' dim(d$counts)          # 20 sites x 33 taxa
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
  example_edna
}


#' Load the built-in replicated example dataset
#'
#' @description
#' The same survey as [get_example_data()], but with three replicates
#' per site instead of one pooled sample. Use it to see how replicated data
#' should be shaped for the hierarchical model.
#'
#' The only structural difference is that `counts` has one row per *replicate*
#' rather than per site, and a `replication` vector says which rows belong
#' together. Do not sum replicates before fitting - see [eDNA_dmm()].
#'
#' @section Data format:
#' \describe{
#'   \item{`counts`}{Integer matrix, 60 replicates x 33 taxa. Row names are
#'     `STN_001_B1`, `STN_001_B2`, ... - three replicates per site.}
#'   \item{`replication`}{Character vector of length 60 giving the site each
#'     row belongs to. Pass this as `replication` to [eDNA_dmm()].}
#'   \item{`covariates`}{Data frame, 20 rows x 2 numeric columns. Covariates
#'     are **station-level**: one row per site, not per replicate.}
#'   \item{`metadata`}{Data frame, 20 rows, station-level labelling columns.}
#' }
#'
#' @return A named list with elements `counts`, `replication`, `covariates`
#'   and `metadata`.
#'
#' @seealso [get_example_data()], [eDNA_dmm()]
#'
#' @examples
#' r <- get_example_replicates()
#'
#' dim(r$counts)            # 60 replicates x 33 taxa
#' r$counts[1:4, 1:5]       # STN_001_B1 .. STN_002_B1
#' head(r$replication, 6)    # which site each row came from
#' nrow(r$covariates)       # 20 - one row per site, not per replicate
#'
#' \dontrun{
#' fit <- eDNA_dmm(
#'   counts     = r$counts,
#'   replication = r$replication,
#'   covariates = r$covariates,
#'   K          = 4
#' )
#' fit$phi_mean   # replicate reproducibility
#' }
#'
#' @export
get_example_replicates <- function() {
  example_edna_reps
}
