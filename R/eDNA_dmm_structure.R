# =============================================================================
# eDNA_dmm_structure() — STRUCTURE-like community assignment bar plots
# =============================================================================

#' Plot posterior community assignments as STRUCTURE-like bar charts
#'
#' @description
#' Produces a stacked bar plot where each vertical bar is one sample, and the
#' bar is divided into colored segments whose heights represent the posterior
#' probability that the sample belongs to each community. This is the classic
#' "STRUCTURE plot" familiar from population genetics.
#'
#' The function returns a [ggplot2::ggplot()] object so you can further
#' customize it (e.g., adjust axis labels, add theme elements, combine with
#' other plots using `cowplot` or `patchwork`).
#'
#' @section Faceting by metadata:
#' The `facet_var` argument accepts a column name from a **metadata data frame**
#' you supply (`metadata` argument). Common choices:
#' - `facet_var = "year"` — panel per sampling year
#' - `facet_var = "depth_bin"` — panel per depth stratum
#' - `facet_var = "site"` — panel per geographic site
#'
#' Within each panel, samples are ordered by `sort_var` (default: latitude).
#'
#' @section Custom colors:
#' By default, community colors are evenly spaced HCL hues. To override with
#' your own palette, pass a named character vector to `community_colors`:
#' ```r
#' my_colors <- c("Community 1" = "#E63946", "Community 2" = "#457B9D")
#' eDNA_dmm_structure(fit, community_colors = my_colors)
#' ```
#' Names must match `"Community 1"`, `"Community 2"`, etc.
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param metadata An optional data frame with one row per sample containing
#'   additional variables for faceting and sorting. Must have a column whose
#'   values match the sample IDs in `fit$sample_info$sample_id`.
#'   Pass `NULL` (the default) to plot all samples in a single panel.
#' @param sample_id_col A string: the name of the column in `metadata` that
#'   matches sample IDs. Only used when `metadata` is not NULL. Default is
#'   `"sample_id"`.
#' @param facet_var A string: the name of a column in `metadata` to use for
#'   faceting (creating separate panels). Pass `NULL` (the default) for a
#'   single panel with all samples. Example: `"year"` or `"depth"`.
#' @param sort_var A string: the name of a column in `metadata` to sort
#'   samples within each panel. Default is `NULL`, which preserves the
#'   original row order. Example: `"latitude"` to sort south-to-north.
#' @param community_colors A named character vector mapping community names
#'   (e.g., `"Community 1"`) to color hex codes. Pass `NULL` (the default)
#'   to use the automatic HCL palette. Must include all K communities.
#' @param bar_width A positive number between 0 and 1: the width of each
#'   sample bar (passed to [ggplot2::geom_bar()]). Default is `0.9`. Increase
#'   toward 1 to remove gaps between bars; decrease for more spacing.
#' @param x_text Logical. If `TRUE`, show sample labels on the x-axis.
#'   Default is `FALSE` (labels are hidden), which is typical when samples
#'   are ordered continuously. Set to `TRUE` if you have few samples and want
#'   to label them.
#' @param base_size A positive number: the base font size (points) for the
#'   ggplot theme. Default is `11`. Increase for larger axis text.
#' @param title A string: the plot title. Default is an auto-generated title
#'   including K. Set to `""` for no title.
#' @param subtitle A string: the plot subtitle. Default is auto-generated.
#'   Set to `""` for no subtitle.
#' @param legend_position A string passed to [ggplot2::theme()]: where to
#'   place the legend. One of `"bottom"` (default), `"right"`, `"left"`,
#'   `"top"`, or `"none"`.
#'
#' @return A [ggplot2::ggplot()] object. Print it to display, or save with
#'   [ggplot2::ggsave()]. The object can be further modified with additional
#'   ggplot2 layers, scales, and theme adjustments.
#'
#' @seealso [eDNA_dmm()], [eDNA_dmm_nmds()], [eDNA_dmm_beta()]
#'
#' @examples
#' \dontrun{
#' data <- get_example_data()
#' fit  <- eDNA_dmm(data$counts, data$covariates, K = 3)
#'
#' # Basic plot — all samples in one panel
#' eDNA_dmm_structure(fit)
#'
#' # Faceted by year, sorted by latitude
#' eDNA_dmm_structure(
#'   fit,
#'   metadata  = data$covariates,
#'   facet_var = "year",
#'   sort_var  = "latitude"
#' )
#'
#' # Custom community colors
#' eDNA_dmm_structure(
#'   fit,
#'   community_colors = c("Community 1" = "#2196F3",
#'                        "Community 2" = "#FF5722",
#'                        "Community 3" = "#4CAF50")
#' )
#'
#' # Save with ggsave
#' p <- eDNA_dmm_structure(fit)
#' ggplot2::ggsave("structure_plot.png", p, width = 10, height = 4, dpi = 300)
#'
#' # Further ggplot2 customization
#' library(ggplot2)
#' p + theme(strip.text = element_text(size = 14, face = "bold"))
#' }
#'
#' @export
eDNA_dmm_structure <- function(
    fit,
    metadata          = NULL,
    sample_id_col     = "sample_id",
    facet_var         = NULL,
    sort_var          = NULL,
    community_colors  = NULL,
    bar_width         = 0.9,
    x_text            = FALSE,
    base_size         = 11,
    title             = NULL,
    subtitle          = NULL,
    legend_position   = "bottom"
) {
  check_stan_fit_object(fit, "eDNA_dmm_structure")

  K         <- fit$K
  si        <- fit$sample_info
  prob_cols <- paste0("prob_comm", seq_len(K))

  # ── Community colors ──────────────────────────────────────────────────────────
  if (is.null(community_colors)) {
    community_colors <- make_community_colors(K)
  } else {
    expected <- paste0("Community ", seq_len(K))
    missing_comms <- setdiff(expected, names(community_colors))
    if (length(missing_comms) > 0) {
      rlang::abort(
        c(
          "`community_colors` is missing names for some communities.",
          i = paste0("Missing: ", paste(missing_comms, collapse = ", ")),
          i = paste0("Expected names: ", paste(expected, collapse = ", "))
        )
      )
    }
  }

  # ── Merge metadata if provided ────────────────────────────────────────────────
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) {
      rlang::abort(
        c(
          "`metadata` must be a data frame.",
          i = paste0("You supplied an object of class: ", paste(class(metadata), collapse = ", "))
        )
      )
    }
    if (!sample_id_col %in% names(metadata)) {
      bad_cols <- setdiff(sample_id_col, names(metadata))
      near     <- names(metadata)[which.min(adist(sample_id_col, names(metadata), ignore.case = TRUE))]
      rlang::abort(
        c(
          paste0("Column '", sample_id_col, "' not found in `metadata`."),
          i = paste0("Available columns: ", paste(names(metadata), collapse = ", ")),
          i = paste0("Did you mean '", near, "'?")
        )
      )
    }
    si <- merge(si, metadata, by.x = "sample_id", by.y = sample_id_col, all.x = TRUE)
  }

  # ── Validate facet_var and sort_var ────────────────────────────────────────────
  if (!is.null(facet_var)) {
    if (!facet_var %in% names(si)) {
      rlang::abort(
        c(
          paste0("`facet_var = '", facet_var, "'` not found."),
          i = if (!is.null(metadata))
            paste0("Available columns in merged data: ", paste(names(si), collapse = ", "))
          else
            "Pass a `metadata` data frame containing this column."
        )
      )
    }
  }
  if (!is.null(sort_var)) {
    if (!sort_var %in% names(si)) {
      rlang::abort(
        c(
          paste0("`sort_var = '", sort_var, "'` not found."),
          i = if (!is.null(metadata))
            paste0("Available columns in merged data: ", paste(names(si), collapse = ", "))
          else
            "Pass a `metadata` data frame containing this column."
        )
      )
    }
  }

  # ── Sort samples ──────────────────────────────────────────────────────────────
  if (!is.null(sort_var)) {
    si <- si[order(si[[sort_var]]), ]
  }

  # ── Pivot to long format for ggplot ───────────────────────────────────────────
  si$x_pos <- seq_len(nrow(si))
  plot_long <- tidyr::pivot_longer(
    si,
    cols      = dplyr::all_of(prob_cols),
    names_to  = "community",
    values_to = "probability"
  )
  plot_long$community <- factor(
    plot_long$community,
    levels = prob_cols,
    labels = paste0("Community ", seq_len(K))
  )

  # ── Build plot ────────────────────────────────────────────────────────────────
  title_str    <- title    %||% sprintf("Posterior Community Assignments  (K = %d)", K)
  subtitle_str <- subtitle %||% sprintf(
    "%d samples  |  bar height = posterior membership probability", nrow(si)
  )

  p <- ggplot2::ggplot(plot_long,
                       ggplot2::aes(x = .data$x_pos,
                                    y = .data$probability,
                                    fill = .data$community)) +
    ggplot2::geom_bar(stat = "identity", width = bar_width) +
    ggplot2::scale_fill_manual(values = community_colors, name = NULL) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      expand = c(0, 0),
      breaks = c(0, 0.5, 1)
    ) +
    ggplot2::labs(
      x        = NULL,
      y        = "Membership probability",
      title    = title_str,
      subtitle = subtitle_str
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position  = legend_position,
      panel.spacing.x  = ggplot2::unit(0.3, "lines"),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = base_size - 2)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

  # ── X-axis labels ─────────────────────────────────────────────────────────────
  if (x_text) {
    p <- p + ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = base_size - 3),
      axis.ticks.x = ggplot2::element_line()
    )
  } else {
    p <- p + ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  # ── Faceting ──────────────────────────────────────────────────────────────────
  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_grid(
      reformulate(facet_var),
      scales = "free_x",
      space  = "free_x"
    )
  }

  p
}

# Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b
