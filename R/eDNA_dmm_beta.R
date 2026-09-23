# =============================================================================
# eDNA_dmm_beta(): Prior vs posterior plots for softmax beta coefficients
# =============================================================================

#' Plot covariate effects: prior vs posterior distributions of beta coefficients
#'
#' @description
#' Produces a prior-versus-posterior density plot for the softmax regression
#' coefficients (beta) that describe how environmental covariates drive
#' community membership. A posterior distribution that is pulled away from the
#' prior is evidence that a covariate genuinely predicts community membership.
#' A posterior that overlaps heavily with the prior means the data contain
#' little information about that covariate's effect.
#'
#' The function returns a [ggplot2::ggplot()] object.
#'
#' @section Joint vs separate layout:
#' The `layout` argument controls how communities are displayed:
#' - `"joint"` (default): all communities overlaid on the same panel per
#'   covariate, each community a different color. Good for comparing community
#'   effects at a glance.
#' - `"separate"`: one row per community (facet rows), one column per
#'   covariate (facet columns). Better when you want to inspect each community
#'   in detail without overlap.
#'
#' @section Reading the plot:
#' For each covariate:
#' - **Grey**: the prior distribution, Normal(0, `beta_prior_sd`). This is
#'   what the model assumed before seeing the data.
#' - **Colored**: the posterior distribution, what the model learned from the
#'   data. One color per community.
#'
#' By default (`reference = "none"`) every coefficient is that community's
#' *deviation from the average community* on that covariate, there is no
#' reference community and all `K` communities are drawn. Interpret:
#' - **Posterior >> 0**: this community becomes more likely, relative to the
#'   average community, at high values of the covariate.
#' - **Posterior << 0**: this community becomes less likely, relative to the
#'   average community, at high values of the covariate.
#' - **Posterior ~ prior**: the data do not constrain this coefficient; the
#'   covariate may not predict community membership.
#'
#' Pass an integer `reference` to instead read every coefficient as a contrast
#' against that one community (which is then omitted, having been re-expressed
#' as exactly zero); the sign interpretation above is the same, just measured
#' against that community instead of the average.
#'
#' @section Coefficient table:
#' The function also prints and invisibly returns a data frame of posterior
#' summaries (mean, 90% CI, probability of direction, ESS, reliability). To
#' access it:
#' ```r
#' result <- eDNA_dmm_beta(fit)
#' result$table
#' ```
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param layout A string: `"joint"` (default) or `"separate"`. Controls
#'   whether communities are overlaid or in separate facet rows. See the
#'   **Joint vs separate layout** section above.
#' @param covariates_to_plot A character vector of covariate names to include.
#'   Default is `NULL`, which plots all covariates (excluding the intercept).
#'   Example: `covariates_to_plot = c("depth", "latitude")`.
#' @param show_intercept Logical. If `FALSE` (the default), the intercept
#'   coefficient is excluded from the plot (it is rarely of direct interest).
#'   Set to `TRUE` to include it.
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
#' @param beta_prior_sd A positive number: the standard deviation of the
#'   Normal(0, `beta_prior_sd`) prior used on beta coefficients in the Stan
#'   model. Default is `1.0`. This **must match** the `to_vector(beta) ~ normal(0, sd)`
#'   line in the Stan model. The default Stan model uses `normal(0, 1.0)`. If
#'   you have changed the prior (e.g., to `normal(0, 2.5)`), update this argument.
#' @param n_prior_samples A positive integer: the number of draws from the
#'   prior distribution to use for the prior density. Default is `20000`.
#'   Higher values give smoother prior curves but make no difference to the
#'   posterior.
#' @param community_colors A named character vector mapping community names
#'   to hex colors. Pass `NULL` (the default) for the automatic HCL palette.
#' @param prior_color A string: the fill color for the prior density.
#'   Default is `"grey60"`.
#' @param prior_alpha A number in (0, 1]: transparency of the prior density.
#'   Default is `0.35`.
#' @param posterior_alpha A number in (0, 1]: transparency of the posterior
#'   density. Default is `0.55`.
#' @param show_annotations Logical. If `TRUE` (the default), annotate each
#'   panel with the posterior mean, 90% CI, and P(direction) value.
#' @param base_size A positive number: the base font size. Default is `13`.
#' @param title A string: the plot title. Default is auto-generated.
#' @param subtitle Plot subtitle. Default `NULL`, meaning no subtitle.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`plot`}{A [ggplot2::ggplot()] object.}
#'   \item{`table`}{A data frame with posterior summaries for each coefficient.
#'     Columns: `community`, `covariate`, `mean`, `median`, `ci_5`, `ci_95`,
#'     `prob_positive`, `prob_negative`, `ess`, `reliability`.}
#' }
#'
#' @seealso [eDNA_dmm()], [eDNA_dmm_structure()], [eDNA_dmm_nmds()]
#'
#' @examples
#' \dontrun{
#' data <- get_example_data()
#' fit  <- eDNA_dmm(data$counts, data$covariates, K = 3)
#'
#' # Default: joint layout, all covariates, prior SD = 1
#' result <- eDNA_dmm_beta(fit)
#' result$plot
#' result$table
#'
#' # Separate layout (each community in its own row)
#' eDNA_dmm_beta(fit, layout = "separate")$plot
#'
#' # Only plot specific covariates
#' eDNA_dmm_beta(fit, covariates_to_plot = c("depth"))$plot
#'
#' # If you changed the beta prior in Stan to normal(0, 2.5):
#' eDNA_dmm_beta(fit, beta_prior_sd = 2.5)$plot
#' }
#'
#' @export
eDNA_dmm_beta <- function(
    fit,
    layout              = "joint",
    covariates_to_plot  = NULL,
    show_intercept      = FALSE,
    beta_prior_sd       = 1.0,
    reference           = "none",
    n_prior_samples     = 20000,
    community_colors    = NULL,
    prior_color         = "grey60",
    prior_alpha         = 0.35,
    posterior_alpha     = 0.55,
    show_annotations    = NULL,
    base_size           = 13,
    title               = NULL,
    subtitle            = NULL
) {
  check_stan_fit_object(fit, "eDNA_dmm_beta")

  K <- fit$K
  if (K < 2) {
    rlang::abort("Beta plots require K >= 2 (at least one non-reference community).")
  }
  if (length(fit$covariate_names) == 0) {
    rlang::abort(
      c(
        "This model has no covariates (intercept-only).",
        i = "Beta plots only apply when covariates were used in `eDNA_dmm()`."
      )
    )
  }

  layout <- match.arg(layout, c("joint", "separate"))

  # -- Community colors ----------------------------------------------------------
  if (is.null(community_colors)) {
    community_colors <- make_community_colors(K)
  }

  # -- Beta draws (relabelled and re-referenced by dmm_beta_draws) --------------
  # Auto-annotate only when K=2 (one community, no overlap)
  if (is.null(show_annotations)) {
    show_annotations <- (K == 2)
  }
  

  cov_labels  <- c("intercept", fit$covariate_names)

  # -- Filter covariates ---------------------------------------------------------
  all_cov_labels <- if (show_intercept) cov_labels else cov_labels[cov_labels != "intercept"]

  if (!is.null(covariates_to_plot)) {
    bad <- setdiff(covariates_to_plot, fit$covariate_names)
    if (length(bad) > 0) {
      rlang::abort(
        c(
          paste0("Covariates not found in model: ", paste(bad, collapse = ", ")),
          i = paste0("Available covariates: ", paste(fit$covariate_names, collapse = ", "))
        )
      )
    }
    all_cov_labels <- intersect(all_cov_labels, covariates_to_plot)
  }

  if (length(all_cov_labels) == 0) {
    rlang::abort(
      c(
        "No covariates to plot after filtering.",
        i = "Check `covariates_to_plot` and `show_intercept` arguments."
      )
    )
  }

  # -- Build posterior long data frame -------------------------------------------
  # Draws come through dmm_beta_draws() rather than straight off the stan_fit,
  # so that a relabelled fit and any change of reference community are both
  # honoured. Reading the raw draws here would silently disagree with
  # eDNA_dmm_beta_intervals() on the same fit.
  # reference = "none" centres across communities instead of subtracting one,
  # so each value is the average of that community's contrasts against all the
  # others and every community is drawn. See eDNA_dmm_beta_intervals().
  centred <- identical(reference, "none")
  ref_use <- if (centred) NULL
             else if (!is.null(reference)) reference else fit$beta_reference
  d <- dmm_beta_draws(fit)
  if (!is.null(fit$community_order))
    d <- d[, fit$community_order, , drop = FALSE]
  if (centred) {
  } else if (!is.null(ref_use)) {
    r <- d[, ref_use, , drop = FALSE]
    for (k in seq_len(dim(d)[2])) d[, k, ] <- d[, k, ] - r[, 1, ]
  }
  ref_idx <- if (centred) NA_integer_
             else if (!is.null(ref_use)) as.integer(ref_use) else K

  posterior_rows <- vector("list", (K - 1) * length(all_cov_labels))
  idx <- 1L
  for (comm_i in seq_len(K)) {
    if (!is.na(ref_idx) && comm_i == ref_idx) next   # reference row is zero
    for (cov_name in all_cov_labels) {
      posterior_rows[[idx]] <- data.frame(
        community = paste0("Community ", comm_i),
        covariate = cov_name,
        value     = d[, comm_i, cov_name],
        source    = "Posterior",
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  posterior_df <- do.call(rbind, Filter(Negate(is.null), posterior_rows))

  # -- Build prior data frame -----------------------------------------------------
  prior_rows <- vector("list", length(all_cov_labels))
  for (j in seq_along(all_cov_labels)) {
    prior_rows[[j]] <- data.frame(
      community = "Prior",
      covariate = all_cov_labels[j],
      value     = stats::rnorm(n_prior_samples, 0, beta_prior_sd),
      source    = "Prior",
      stringsAsFactors = FALSE
    )
  }
  prior_df <- do.call(rbind, prior_rows)

  # -- Combine -------------------------------------------------------------------
  plot_df <- rbind(posterior_df, prior_df)

  # -- Annotation data frame -----------------------------------------------------
  if (show_annotations) {
    ann_rows <- vector("list", (K - 1) * length(all_cov_labels))
    idx <- 1L
    for (comm_i in seq_len(K)) {
      if (!is.na(ref_idx) && comm_i == ref_idx) next
      for (cov_name in all_cov_labels) {
        draws  <- d[, comm_i, cov_name]
        p_dir  <- max(mean(draws > 0), mean(draws < 0))
        ann_rows[[idx]] <- data.frame(
          community = paste0("Community ", comm_i),
          covariate = cov_name,
          label     = sprintf(
            "\u03b2 = %.2f\n[%.2f, %.2f]\nP(dir) = %.0f%%",
            mean(draws),
            stats::quantile(draws, 0.05),
            stats::quantile(draws, 0.95),
            100 * p_dir
          ),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
    ann_df <- do.call(rbind, Filter(Negate(is.null), ann_rows))
  }

  # -- Build plot ----------------------------------------------------------------
  # Color map: communities + "Prior" as grey
  comm_names    <- setdiff(paste0("Community ", seq_len(K)),
                           if (is.na(ref_idx)) character(0) else paste0("Community ", ref_idx))
  color_map     <- c(community_colors[comm_names], Prior = prior_color)
  alpha_map_vec <- c(
    stats::setNames(rep(posterior_alpha, length(comm_names)), comm_names),
    Prior = prior_alpha
  )

  title_str    <- title    %||% sprintf("Covariate effects on community membership  (K = %d)", K)
  subtitle_str <- subtitle

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$value, fill = .data$community, color = .data$community)
  ) +
    ggplot2::geom_density(
      ggplot2::aes(alpha = .data$community),
      linewidth = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      color      = "grey30",
      linewidth  = 0.4
    ) +
    ggplot2::scale_fill_manual(
      values = color_map,
      name   = NULL
    ) +
    ggplot2::scale_color_manual(
      values = color_map,
      name   = NULL
    ) +
    ggplot2::scale_alpha_manual(
      values = alpha_map_vec,
      guide  = "none"
    ) +
    ggplot2::labs(
      x        = "Coefficient (Z-scored covariate)",
      y        = "Density",
      title    = title_str,
      subtitle = subtitle_str
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      legend.position  = "right",
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = base_size - 2)
    )

  # -- Annotations ---------------------------------------------------------------
  if (show_annotations) {
    p <- p + ggplot2::geom_text(
      data        = ann_df,
      ggplot2::aes(x = Inf, y = Inf, label = .data$label),
      hjust       = 1.05,
      vjust       = 1.2,
      size        = 3,
      color       = "grey20",
      inherit.aes = FALSE
    )
  }

  # -- Faceting by layout --------------------------------------------------------
  if (layout == "joint") {
    # One panel per covariate; communities overlaid
    p <- p + ggplot2::facet_wrap(~ .data$covariate, scales = "free", nrow = 1)
  } else {
    # One row per community, one column per covariate
    p <- p + ggplot2::facet_grid(
      .data$community ~ .data$covariate,
      scales = "free"
    )
  }

  # -- Return --------------------------------------------------------------------
  # Filter beta_summary to the user-selected covariates
  beta_tbl <- fit$beta_summary
  if (!show_intercept) {
    beta_tbl <- beta_tbl[beta_tbl$covariate != "intercept", ]
  }
  if (!is.null(covariates_to_plot)) {
    beta_tbl <- beta_tbl[beta_tbl$covariate %in% covariates_to_plot, ]
  }

  list(plot = p, table = beta_tbl)
}
