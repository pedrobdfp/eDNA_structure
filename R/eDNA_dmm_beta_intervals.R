# =============================================================================
# eDNA_dmm_beta_intervals(): coefficient plot for softmax beta coefficients
# =============================================================================

#' Plot covariate effects as posterior intervals
#'
#' @description
#' Draws the softmax regression coefficients (beta) as a coefficient plot: one
#' row per community and covariate, a point at the posterior point estimate,
#' two nested credible intervals, and a dashed line at zero.
#'
#' This is a companion to [eDNA_dmm_beta()], which overlays prior and posterior
#' densities. Densities show how far the data moved a coefficient from its
#' prior, but the tail of a density is hard to read against zero. This function
#' answers the complementary question directly: does the interval cross zero,
#' and by how much.
#'
#' @section Reading the plot:
#' Each row is one coefficient, relative to the reference community:
#' - The **thick line** is the inner interval (50% by default).
#' - The **thin line** is the outer interval (90% by default).
#' - The **point** is filled when the outer interval excludes zero and hollow
#'   when it crosses zero.
#'
#' A coefficient whose outer interval lies entirely above zero means the
#' community becomes more likely as the covariate increases, relative to the
#' reference community; entirely below zero means less likely.
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param covariates_to_plot Character vector of covariate names to show, or
#'   `NULL` (default) for all of them.
#' @param show_intercept Logical. Include the intercept coefficients.
#'   Default `FALSE`.
#' @param intervals A length-2 numeric vector giving the inner (thick) and
#'   outer (thin) interval widths, each strictly between 0 and 1. Default
#'   `c(0.5, 0.9)`. Order does not matter; the wider one is drawn thin.
#' @param reference How the coefficients are identified for display.
#'   `"none"` (default) centres them across communities, so each
#'   value is that community's deviation from the average community and all
#'   `K` communities are drawn. An integer instead reads every
#'   coefficient as a contrast against that community, which is then omitted
#'   as the reference. `NULL` uses the fit's own `beta_reference` if it has one.
#'
#'   Changing this requires no refitting. The softmax is invariant to adding
#'   a constant across communities, so the choice is a normalisation applied
#'   to the finished draws, not an estimate.
#' @param point_est One of `"median"` (default) or `"mean"`: the posterior
#'   summary the point marks.
#' @param community_colors Named character vector of colors, or `NULL`
#'   (default) to use the package palette. Points and lines are colored by
#'   community; no legend is drawn because the y axis already names them.
#' @param color_by_community Logical. Color rows by community. Default `TRUE`.
#'   `FALSE` draws everything in a single dark color.
#' @param facet_scales Passed to [ggplot2::facet_wrap()]. Default `"free_x"`,
#'   so each covariate gets its own x range. Use `"fixed"` to compare
#'   magnitudes across covariates on a common scale.
#' @param point_size Size of the point estimate marker. Default `2.8`.
#' @param linewidth_inner,linewidth_outer Line widths of the inner and outer
#'   intervals. Defaults `1.6` and `0.6`.
#' @param base_size Base font size. Default `13`.
#' @param title,subtitle Plot title and subtitle, or `NULL`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`plot`}{A [ggplot2::ggplot()] object.}
#'   \item{`table`}{A data frame with one row per coefficient: `community`,
#'     `covariate`, `estimate`, the four interval bounds, and `excludes_zero`.}
#' }
#'
#' @seealso [eDNA_dmm_beta()] for the prior-versus-posterior density version.
#'
#' @examples
#' \dontrun{
#' fit <- eDNA_dmm(counts, covariates, K = 3)
#'
#' # Default 50% and 90% intervals
#' eDNA_dmm_beta_intervals(fit)$plot
#'
#' # 50% and 95% instead
#' eDNA_dmm_beta_intervals(fit, intervals = c(0.5, 0.95))$plot
#' }
#'
#' @export
eDNA_dmm_beta_intervals <- function(
    fit,
    covariates_to_plot = NULL,
    show_intercept     = FALSE,
    intervals          = c(0.5, 0.9),
    reference          = "none",
    point_est          = c("median", "mean"),
    community_colors   = NULL,
    color_by_community = TRUE,
    facet_scales       = "free_x",
    point_size         = 2.8,
    linewidth_inner    = 1.6,
    linewidth_outer    = 0.6,
    base_size          = 13,
    title              = NULL,
    subtitle           = NULL
) {
  check_stan_fit_object(fit, "eDNA_dmm_beta_intervals")
  point_est <- match.arg(point_est)

  K <- fit$K
  if (K < 2)
    rlang::abort("Coefficient plots require K >= 2 (at least one non-reference community).")
  if (length(fit$covariate_names) == 0)
    rlang::abort(c(
      "This model has no covariates (intercept-only).",
      i = "Coefficient plots only apply when covariates were used in `eDNA_dmm()`."
    ))

  if (!is.numeric(intervals) || length(intervals) != 2 ||
      any(intervals <= 0) || any(intervals >= 1))
    rlang::abort("`intervals` must be two numbers strictly between 0 and 1, e.g. c(0.5, 0.9).")

  # The wider interval is always the thin outer line, whichever order it came in.
  inner <- min(intervals)
  outer <- max(intervals)
  if (isTRUE(all.equal(inner, outer)))
    rlang::abort("`intervals` must give two different widths.")

  q_inner <- c((1 - inner) / 2, 1 - (1 - inner) / 2)
  q_outer <- c((1 - outer) / 2, 1 - (1 - outer) / 2)

  # ── Which coefficients ──────────────────────────────────────────────────────
  cov_labels <- c("intercept", fit$covariate_names)
  keep <- if (show_intercept) cov_labels else setdiff(cov_labels, "intercept")
  if (!is.null(covariates_to_plot)) {
    unknown <- setdiff(covariates_to_plot, cov_labels)
    if (length(unknown) > 0)
      rlang::abort(c(
        "Unknown covariate name(s) in `covariates_to_plot`.",
        i = paste0("Not in the model: ", paste(unknown, collapse = ", ")),
        i = paste0("Available: ", paste(cov_labels, collapse = ", "))
      ))
    keep <- intersect(keep, covariates_to_plot)
  }
  if (length(keep) == 0)
    rlang::abort(c(
      "No covariates to plot after filtering.",
      i = "Check `covariates_to_plot` and `show_intercept`."
    ))

  # ── Summarise the draws ─────────────────────────────────────────────────────
  # Goes through dmm_beta_draws() so that a fit relabelled by
  # dmm_relabel_communities() and any re-referencing are both honoured.
  # reference = "none" centres the coefficients across communities instead of
  # subtracting one of them. Because sum_j (beta_k - beta_j) = K * beta_k when
  # the coefficients are centred, each value is then the average of that
  # community's contrasts against all the others - a function of identified
  # contrasts only, depending on no choice of reference. Every community is
  # therefore estimated and none is dropped from the figure.
  centred <- identical(reference, "none")
  ref_use <- if (centred) NULL
             else if (!is.null(reference)) reference else fit$beta_reference
  d <- dmm_beta_draws(fit)
  if (!is.null(fit$community_order))
    d <- d[, fit$community_order, , drop = FALSE]
  if (centred) {
    for (j in seq_len(dim(d)[3])) d[, , j] <- d[, , j] - rowMeans(d[, , j])
  } else if (!is.null(ref_use)) {
    r <- d[, ref_use, , drop = FALSE]
    for (k in seq_len(dim(d)[2])) d[, k, ] <- d[, k, ] - r[, 1, ]
  }
  ref_idx <- if (centred) NA_integer_
             else if (!is.null(ref_use)) as.integer(ref_use) else K

  rows <- list()
  for (comm_i in seq_len(K)) {
    if (!is.na(ref_idx) && comm_i == ref_idx) next   # reference row is zero
    for (cov_name in keep) {
      draws <- d[, comm_i, cov_name]
      qi <- stats::quantile(draws, q_inner, names = FALSE)
      qo <- stats::quantile(draws, q_outer, names = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(
        community     = paste0("Community ", comm_i),
        covariate     = cov_name,
        estimate      = if (point_est == "median") stats::median(draws) else mean(draws),
        inner_lo      = qi[1], inner_hi = qi[2],
        outer_lo      = qo[1], outer_hi = qo[2],
        excludes_zero = (qo[1] > 0) || (qo[2] < 0),
        stringsAsFactors = FALSE
      )
    }
  }
  tbl <- do.call(rbind, rows)

  # Communities run top to bottom in order; covariates keep model order.
  comm_levels    <- setdiff(paste0("Community ", seq_len(K)),
                            if (is.na(ref_idx)) character(0)
                            else paste0("Community ", ref_idx))
  tbl$community  <- factor(tbl$community, levels = rev(comm_levels))
  tbl$covariate  <- factor(tbl$covariate, levels = keep)

  # ── Colors ──────────────────────────────────────────────────────────────────
  if (is.null(community_colors)) community_colors <- make_community_colors(K)
  pal <- if (isTRUE(color_by_community)) {
    community_colors[comm_levels]
  } else {
    stats::setNames(rep("#1D3557", length(comm_levels)), comm_levels)
  }

  # ── Plot ────────────────────────────────────────────────────────────────────
  # `fill` carries credibility and `colour` carries community, so shape 21 is
  # used throughout: a filled point is a coefficient whose outer interval
  # excludes zero, a hollow one is a coefficient that crosses it.
  p <- ggplot2::ggplot(
      tbl,
      ggplot2::aes(y = .data$community, x = .data$estimate,
                   colour = .data$community)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey40", linewidth = 0.6) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$outer_lo, xmax = .data$outer_hi),
      linewidth = linewidth_outer) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$inner_lo, xmax = .data$inner_hi),
      linewidth = linewidth_inner) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$excludes_zero),
      shape = 21, size = point_size, stroke = 0.9) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
    ggplot2::facet_wrap(~ .data$covariate, scales = facet_scales, nrow = 1) +
    ggplot2::labs(
      x = sprintf("Coefficient (%.0f%% and %.0f%% credible intervals)",
                  100 * inner, 100 * outer),
      y = NULL, title = title, subtitle = subtitle) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.background   = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold"),
      plot.subtitle      = ggplot2::element_text(color = "grey40",
                                                 size = base_size - 2)
    )

  tbl$community <- factor(as.character(tbl$community), levels = comm_levels)
  list(plot = p, table = tbl[order(tbl$covariate, tbl$community), ])
}
