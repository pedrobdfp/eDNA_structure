# =============================================================================
# eDNA_dmm_structure() — STRUCTURE-like community assignment bar plots
# =============================================================================

#' Plot posterior community assignments as STRUCTURE-like bar charts
#'
#' @description
#' Produces a stacked bar plot where each vertical bar is one sample, and the
#' bar is divided into colored segments whose heights represent the posterior
#' probability that the sample belongs to each community. It uses the same
#' visual format as [plot_true_compositions()] so the two plots can be
#' compared directly: one shows what species are present (observed
#' composition), the other shows what community they were assigned to
#' (posterior assignment probabilities).
#'
#' When `facet_row_var` is supplied, one sub-plot is built per row-level and
#' assembled vertically with cowplot — the only architecture that keeps panel
#' widths proportional to sample count via space = "free_x".
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param metadata An optional data frame with one row per sample.
#' @param sample_id_col Column in `metadata` matching sample IDs. Default `"sample_id"`.
#' @param facet_var Column for column facets (e.g. `"year"`). Default `NULL`.
#' @param facet_row_var Column for row facets (e.g. `"depth_bin"`). When
#'   supplied, one sub-plot is built per level and stacked with cowplot.
#'   Default `NULL`.
#' @param sort_var Column to sort samples within panels. Default `NULL`.
#' @param community_colors Named character vector mapping community names to
#'   colors. Default `NULL` (automatic HCL palette).
#' @param x_text Logical. Show sample labels on x-axis. Default `FALSE`.
#' @param base_size Base font size. Default `11`.
#' @param title Plot title. Default auto-generated.
#' @param subtitle Plot subtitle. Default auto-generated.
#' @param show_legend Logical. Show legend. Default `TRUE`.
#' @param legend_position Legend position. Default `"bottom"`.
#' @param legend_nrow Number of legend rows, passed to `guide_legend()`.
#'   Default `NULL` (ggplot chooses).
#' @param legend_ncol Number of legend columns, passed to `guide_legend()`.
#'   Default `NULL` (ggplot chooses).
#' @param legend_key_size Size (cm) of each legend color swatch. Default `0.4`.
#' @param bar_width Width of each bar, in units of one x slot. Default `0.9`,
#'   which leaves a thin white gutter so individual samples stay countable
#'   and matches [plot_true_compositions()]. Use `1` for a gapless block.
#' @param legend_rel_height Height of a top or bottom legend, as a fraction of
#'   the panel stack. Default `0.1`. Raise it when the legend has many rows:
#'   the legend is given exactly this much room, and rows that do not fit are
#'   clipped without warning.
#' @param legend_rel_width Width of a left or right legend, as a fraction of
#'   the panel stack. Default `0.2`.
#' @param strip_text_size Font size of the column-facet strip labels (e.g. the
#'   years). Default `NULL`, meaning `base_size`.
#' @param row_label_size Font size of the per-row label produced by
#'   `row_label_fn` (e.g. the depths), when `row_label_position = "title"`.
#'   Default `NULL`, meaning `base_size`.
#' @param ylab Shared y-axis label drawn once beside the panel stack, used only
#'   when `row_label_position = "title"`. Default `"Membership probability"`.
#' @param ylab_size Font size of that label. Default `NULL`, meaning
#'   `base_size`.
#' @param ylab_rel_width Width of the label strip, as a fraction of the panel
#'   stack. Default `0.03`.
#' @param axis_text_size Font size of the y-axis tick labels. Default `NULL`,
#'   meaning ggplot2's own scaling from `base_size`.
#' @param legend_text_size Font size of legend labels. Default `NULL`
#'   (inherits from `base_size`).
#' @param row_label_fn Function applied to each `facet_row_var` level to
#'   produce its display label — e.g. `function(x) paste0(x, " m")` for a
#'   depth column, or `function(x) paste0("Year: ", x)` for year.
#'   Default `as.character` (just prints the raw value).
#' @param panel_labels Logical. Label each row-facet panel with a lowercase
#'   letter. Default `TRUE`. Set `FALSE` when the plot is itself a panel of a
#'   larger composite figure, where the outer figure supplies the letters and
#'   an inner set would collide with them.
#' @param row_label_position `"title"` (bold label above each row's panel,
#'   default) or `"ylab"` (as the y-axis label, the old behavior).
#' @param vline_var Numeric column for vertical reference line. Default `NULL`.
#' @param vline_value Threshold value in `vline_var` space. Default `NULL`.
#' @param vline_color Line color. Default `"black"`.
#' @param vline_linetype Line type. Default `"dashed"`.
#' @param vline_linewidth Line width. Default `0.7`.
#'
#' @return A ggplot2 object (no `facet_row_var`) or a cowplot grid object.
#'
#' @seealso [eDNA_dmm()], [plot_true_compositions()], [eDNA_dmm_nmds()], [eDNA_dmm_beta()]
#' @export
eDNA_dmm_structure <- function(
    fit,
    metadata            = NULL,
    sample_id_col       = "sample_id",
    facet_var           = NULL,
    facet_row_var       = NULL,
    sort_var            = NULL,
    community_colors    = NULL,
    x_text              = FALSE,
    base_size           = 11,
    title               = NULL,
    subtitle            = NULL,
    show_legend         = TRUE,
    legend_position     = "bottom",
    legend_nrow         = NULL,
    legend_ncol         = NULL,
    legend_key_size     = 0.4,
    legend_text_size    = NULL,
    legend_rel_height   = 0.1,
    legend_rel_width    = 0.2,
    bar_width           = 0.9,
    strip_text_size     = NULL,
    row_label_size      = NULL,
    ylab                = "Membership probability",
    ylab_size           = NULL,
    ylab_rel_width      = 0.03,
    axis_text_size      = NULL,
    row_label_fn        = as.character,
    row_label_position  = c("title", "ylab"),
    panel_labels        = TRUE,
    vline_var           = NULL,
    vline_value         = NULL,
    vline_color         = "black",
    vline_linetype      = "dashed",
    vline_linewidth     = 0.7
) {
  check_stan_fit_object(fit, "eDNA_dmm_structure")
  row_label_position <- match.arg(row_label_position)

  K         <- fit$K
  si        <- fit$sample_info
  prob_cols <- paste0("prob_comm", seq_len(K))

  # ── Community colors ──────────────────────────────────────────────────────
  if (is.null(community_colors)) {
    community_colors <- make_community_colors(K)
  } else {
    expected      <- paste0("Community ", seq_len(K))
    missing_comms <- setdiff(expected, names(community_colors))
    if (length(missing_comms) > 0)
      rlang::abort(c("`community_colors` missing names.",
                     i = paste0("Missing: ", paste(missing_comms, collapse = ", "))))
  }

  # ── Merge metadata ────────────────────────────────────────────────────────
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata))
      rlang::abort("`metadata` must be a data frame.")
    if (!sample_id_col %in% names(metadata)) {
      near <- names(metadata)[which.min(adist(sample_id_col, names(metadata),
                                              ignore.case = TRUE))]
      rlang::abort(c(paste0("Column '", sample_id_col, "' not found in `metadata`."),
                     i = paste0("Did you mean '", near, "'?")))
    }
    si <- merge(si, metadata, by.x = "sample_id", by.y = sample_id_col, all.x = TRUE)
  }

  # ── Validate ──────────────────────────────────────────────────────────────
  for (v in c(facet_var, facet_row_var, sort_var, vline_var)) {
    if (!is.null(v) && !v %in% names(si))
      rlang::abort(paste0("'", v, "' not found in merged data. Available: ",
                          paste(names(si), collapse = ", ")))
  }
  if (!is.null(vline_var) && !is.numeric(si[[vline_var]]))
    rlang::abort("`vline_var` must be numeric.")

  # ── Titles ────────────────────────────────────────────────────────────────
  title_str    <- title    %||% sprintf("Posterior Community Assignments  (K = %d)", K)
  subtitle_str <- subtitle %||% sprintf(
    "%d samples  |  bar height = posterior membership probability", nrow(si))

  # ── Core single-panel builder ─────────────────────────────────────────────
  # Uses sample_id as discrete x within each subset — gapless bars guaranteed.
  build_panel <- function(dat, row_label = NULL, show_legend = FALSE,
                          show_ylab = TRUE) {
    if (!is.null(sort_var)) dat <- dat[order(dat[[sort_var]]), ]

    # Ordered factor from this subset only — key to gapless bars
    dat$x_label <- factor(dat$sample_id, levels = unique(dat$sample_id))

    # Vline position within facet_var panels
    vline_df <- NULL
    if (!is.null(vline_var) && !is.null(vline_value) && vline_var %in% names(dat)) {
      if (!is.null(facet_var) && facet_var %in% names(dat)) {
        vline_df <- dat |>
          dplyr::select(dplyr::all_of(c(facet_var, vline_var, "x_label"))) |>
          dplyr::distinct() |>
          dplyr::group_by(dplyr::across(dplyr::all_of(facet_var))) |>
          dplyr::summarise(
            # The position must be LOCAL to the panel. facet_grid(scales =
            # "free_x") drops unused levels and re-indexes each panel from 1, so
            # a global level index would land outside the panel and stretch its
            # x range — which silently destroys the proportional panel widths
            # that space = "free_x" is supposed to give.
            cape_x = match(
              .data$x_label[which.min(abs(.data[[vline_var]] - vline_value))],
              sort(unique(.data$x_label))
            ) + 0.5,
            .groups = "drop"
          )
      } else {
        closest  <- which.min(abs(dat[[vline_var]] - vline_value))
        vline_df <- data.frame(cape_x = as.integer(dat$x_label[closest]) + 0.5)
      }
    }

    plot_long <- tidyr::pivot_longer(dat, cols = dplyr::all_of(prob_cols),
                                     names_to = "community", values_to = "probability")
    plot_long$community <- factor(plot_long$community, levels = prob_cols,
                                  labels = paste0("Community ", seq_len(K)))
    # ensure x_label factor levels are preserved after pivot
    plot_long$x_label <- factor(plot_long$x_label, levels = levels(dat$x_label))

    panel_title <- if (!is.null(row_label) && row_label_position == "title") row_label else NULL
    # In "ylab" mode the row label IS the y-axis title. In "title" mode the axis
    # title would otherwise repeat on every stacked row panel, so callers ask for
    # it on one panel only.
    panel_ylab  <- if (!is.null(row_label) && row_label_position == "ylab") {
      row_label
    } else if (show_ylab) {
      ylab
    } else {
      NULL
    }

    p <- ggplot2::ggplot(
      plot_long,
      ggplot2::aes(x = .data$x_label, y = .data$probability, fill = .data$community)
    ) +
      # bar_width controls the gutter between samples. Below 1 leaves white
      # striping so a single station can be picked out of a run of hundreds;
      # 1 butts the bars together into a continuous block.
      ggplot2::geom_bar(stat = "identity", position = "stack", color = NA,
                        width = bar_width) +
      ggplot2::scale_fill_manual(values = community_colors, name = NULL) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0),
                                  breaks = c(0, 0.5, 1)) +
      ggplot2::labs(x = NULL, y = panel_ylab, title = panel_title) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = legend_nrow, ncol = legend_ncol)) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        legend.position  = if (show_legend) legend_position else "none",
        legend.key.size  = ggplot2::unit(legend_key_size, "cm"),
        legend.text      = if (!is.null(legend_text_size)) ggplot2::element_text(size = legend_text_size) else ggplot2::element_text(),
        axis.text.y      = ggplot2::element_text(size = axis_text_size %||% base_size),
        panel.spacing.x  = ggplot2::unit(0.3, "lines"),
        strip.background = ggplot2::element_blank(),
        strip.text       = ggplot2::element_text(face = "bold", size = strip_text_size %||% base_size),
        plot.title       = ggplot2::element_text(face = "bold", size = row_label_size %||% base_size, hjust = 0),
        plot.subtitle    = ggplot2::element_text(color = "grey40", size = base_size - 2)
      )

    if (!x_text) {
      p <- p + ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                              axis.ticks.x = ggplot2::element_blank())
    }

    if (!is.null(vline_df)) {
      p <- p + ggplot2::geom_vline(
        data = vline_df,
        ggplot2::aes(xintercept = .data$cape_x),
        color     = vline_color, linetype  = vline_linetype,
        linewidth = vline_linewidth, inherit.aes = FALSE
      )
    }

    if (!is.null(facet_var) && facet_var %in% names(dat)) {
      p <- p + ggplot2::facet_grid(reformulate(facet_var),
                                   scales = "free_x", space = "free_x")
    }

    p
  }

  # ── No row facet: single ggplot ───────────────────────────────────────────
  if (is.null(facet_row_var)) {
    p <- build_panel(si, show_legend = show_legend)
    return(p + ggplot2::labs(title = title_str, subtitle = subtitle_str))
  }

  # ── Row facet: cowplot assembly ───────────────────────────────────────────
  row_levels <- if (is.factor(si[[facet_row_var]])) levels(si[[facet_row_var]]) else
    sort(unique(si[[facet_row_var]]))

  # In "ylab" mode each panel's y-axis title IS its row label. In "title" mode
  # the row label moves to the panel title, so no panel carries an axis title —
  # a shared one is drawn for the whole stack below.
  panels <- Filter(Negate(is.null), lapply(row_levels, function(lv) {
    dat <- si[si[[facet_row_var]] == lv, ]
    if (nrow(dat) == 0) return(NULL)
    build_panel(dat, row_label = row_label_fn(lv), show_legend = FALSE,
                show_ylab = row_label_position == "ylab")
  }))

  # An empty title AND subtitle means the caller is embedding this plot as a
  # panel of a larger figure, where the heading lives in the figure legend. In
  # that case the title row is dropped entirely rather than drawn empty — an
  # empty row still consumes height and shows up as a gap above the panel.
  has_heading <- nzchar(title_str) || nzchar(subtitle_str)

  title_grob <- if (has_heading) {
    cowplot::ggdraw() +
      cowplot::draw_label(title_str, fontface = "bold", size = base_size + 1,
                          x = 0.02, hjust = 0) +
      cowplot::draw_label(subtitle_str, size = base_size - 1, color = "grey40",
                          x = 0.02, y = 0.25, hjust = 0)
  } else NULL

  # Helper: stack the title row on top only when there is one to draw.
  with_title <- function(...) {
    parts   <- list(...)
    heights <- attr(parts, "heights")
    if (is.null(title_grob)) return(parts[[1]])
    cowplot::plot_grid(plotlist = c(list(title_grob), parts), ncol = 1,
                       rel_heights = c(0.06, rep(1, length(parts))))
  }

  stacked <- cowplot::plot_grid(plotlist = panels, ncol = 1,
                                labels = if (isTRUE(panel_labels)) letters[seq_along(panels)] else NULL,
                                label_size = 14)

  # One shared, rotated y-axis label for the stack. Putting the axis title on
  # the panels themselves either repeats it down every row or gets clipped by
  # cowplot's alignment, so it is drawn once here instead.
  if (row_label_position == "title") {
    ylab_grob <- cowplot::ggdraw() +
      cowplot::draw_label(ylab, angle = 90, size = ylab_size %||% base_size)
    stacked <- cowplot::plot_grid(ylab_grob, stacked, nrow = 1,
                                  rel_widths = c(ylab_rel_width, 1))
  }

  if (!show_legend) {
    return(with_title(stacked))
  }

  legend_plot <- build_panel(si, show_legend = TRUE) +
    ggplot2::theme(legend.position = legend_position) +
    ggplot2::guides(fill = ggplot2::guide_legend())

  if (legend_position == "right") {
    legend_grob <- cowplot::get_plot_component(legend_plot, "guide-box-right", return_all = TRUE)
    inner <- cowplot::plot_grid(stacked, legend_grob, nrow = 1, rel_widths = c(1, legend_rel_width))
    with_title(inner)
  } else if (legend_position == "left") {
    legend_grob <- cowplot::get_plot_component(legend_plot, "guide-box-left", return_all = TRUE)
    inner <- cowplot::plot_grid(legend_grob, stacked, nrow = 1, rel_widths = c(legend_rel_width, 1))
    with_title(inner)
  } else {
    guide_box_name <- if (legend_position == "top") "guide-box-top" else "guide-box-bottom"
    legend_grob <- cowplot::get_plot_component(legend_plot, guide_box_name, return_all = TRUE)
    if (legend_position == "top") {
      if (is.null(title_grob)) cowplot::plot_grid(legend_grob, stacked, ncol = 1, rel_heights = c(legend_rel_height, 1)) else cowplot::plot_grid(title_grob, legend_grob, stacked, ncol = 1, rel_heights = c(0.06, legend_rel_height, 1))
    } else {
      if (is.null(title_grob)) cowplot::plot_grid(stacked, legend_grob, ncol = 1, rel_heights = c(1, legend_rel_height)) else cowplot::plot_grid(title_grob, stacked, legend_grob, ncol = 1, rel_heights = c(0.06, 1, legend_rel_height))
    }
  }
}
