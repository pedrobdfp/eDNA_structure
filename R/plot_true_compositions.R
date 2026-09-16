# =============================================================================
# plot_true_compositions(): visualize observed species composition per sample
# =============================================================================

#' Plot observed species composition per sample (pre-fit diagnostic)
#'
#' @description
#' Produces a stacked bar plot of **observed species frequencies**: one bar
#' per sample, colored by species. This is a pre-fitting diagnostic that shows
#' the raw community signal in your data before any model is applied. It uses
#' the same visual format as [eDNA_dmm_structure()] so the two plots can be
#' compared directly: the top panel shows what species are present, the bottom
#' shows what community they were assigned to.
#'
#' When `facet_row_var` is supplied, one sub-plot is built per row-level and
#' they are assembled vertically with cowplot. This mirrors the original
#' make_bar_figure() approach that gives each depth stratum its own
#' proportional-width panels, the only architecture that keeps panel widths
#' proportional to sample count via space = "free_x".
#'
#' @param counts A numeric matrix or data frame of read counts (samples × taxa).
#' @param metadata An optional data frame with one row per sample.
#' @param sample_id_col Column in `metadata` matching `rownames(counts)`.
#'   Default `"sample_id"`.
#' @param facet_var Column for column facets (e.g. `"year"`). Default `NULL`.
#' @param facet_row_var Column for row facets (e.g. `"depth_bin"`). When
#'   supplied, one sub-plot is built per level and stacked with cowplot.
#'   Default `NULL`.
#' @param sort_var Column to sort samples within panels. Default `NULL`.
#' @param top_n Number of most abundant taxa shown individually; rest collapsed
#'   to `"Other"`. Ranked by mean observed frequency across samples. Ignored
#'   when `taxa_include` is supplied. Default `20`.
#' @param taxa_include Character vector naming exactly which taxa to show
#'   individually, everything else pooled into `"Other"`. Supply the result of
#'   [dmm_taxon_order()] here to make this figure share a taxon set, and
#'   therefore a colour palette, with [eDNA_dmm_compositions()]. Names absent
#'   from `counts` are dropped with a warning. Default `NULL`, meaning rank by
#'   `top_n`.
#' @param base_size Base font size. Default `11`.
#' @param title Plot title. Default auto-generated.
#' @param subtitle Plot subtitle. Default `NULL`, meaning no subtitle.
#' @param show_legend Logical. Show the taxon legend. Default `FALSE` - a
#'   full taxon legend is usually too large to be useful on this plot.
#' @param legend_position Legend position. Default `"bottom"`.
#' @param panel_labels Logical. Label each row-facet panel with a lowercase
#'   letter. Default `TRUE`. Set `FALSE` when the plot is itself a panel of a
#'   larger composite figure, where the outer figure supplies the letters and
#'   an inner set would collide with them.
#' @param vline_var Numeric column for vertical reference line. Default `NULL`.
#' @param vline_value Threshold value in `vline_var` space. Default `NULL`.
#' @param vline_color Line color. Default `"black"`.
#' @param vline_linetype Line type. Default `"dashed"`.
#' @param vline_linewidth Line width. Default `0.7`.
#' @param legend_nrow Number of legend rows, passed to `guide_legend()`.
#'   Default `NULL` (ggplot chooses).
#' @param legend_ncol Number of legend columns, passed to `guide_legend()`.
#'   Default `NULL` (ggplot chooses).
#' @param legend_key_size Size (cm) of each legend color swatch. Default `0.4`.
#' @param n_hues Number of hue families in the taxon palette. Taxa are sorted
#'   and laid out hue block by hue block, so a legend of `n_hues` columns puts
#'   roughly one colour family per column. Default `7`.
#' @param bar_width Width of each bar, in units of one x slot. Default `0.9`,
#'   which leaves a thin white gutter so individual samples stay countable.
#'   Use `1` for a gapless STRUCTURE-style block.
#' @param ylab Shared y-axis label drawn once beside the panel stack, used only
#'   when `row_label_position = "title"`. Default `"Proportion"`.
#' @param strip_text_size Font size of the column-facet strip labels (e.g. the
#'   years). Default `NULL`, meaning `base_size`.
#' @param row_label_size Font size of the per-row label produced by
#'   `row_label_fn` (e.g. the depths), when `row_label_position = "title"`.
#'   Default `NULL`, meaning `base_size`.
#' @param ylab_size Font size of the shared y-axis label. Default `NULL`,
#'   meaning `base_size`.
#' @param ylab_rel_width Width of that label strip, as a fraction of the panel
#'   stack. Default `0.03`.
#' @param legend_rel_height Height of a top or bottom legend, as a fraction of
#'   the panel stack. Default `0.1`. Raise it when the legend has many rows:
#'   the legend is given exactly this much room, and rows that do not fit are
#'   clipped without warning.
#' @param legend_rel_width Width of a left or right legend, as a fraction of
#'   the panel stack. Default `0.2`.
#' @param legend_text_size Font size of legend labels. Default `NULL`
#'   (inherits from `base_size`).
#' @param row_label_fn Function applied to each `facet_row_var` level to
#'   produce its display label, e.g. `function(x) paste0(x, " m")` for a
#'   depth column, or `function(x) paste0("Year: ", x)` for year. This is
#'   what makes the label text generalize across datasets: it's tied to
#'   whatever `facet_row_var` means for your data.
#'   Default `as.character` (just prints the raw value).
#' @param row_label_position `"title"` (bold label above each row's panel,
#'   default) or `"ylab"` (as the y-axis label, the old behavior).
#' @return A ggplot2 object (no `facet_row_var`) or a cowplot grid object.
#'
#' @seealso [eDNA_dmm_structure()], [simulate_eDNA_survey()]
#' @export
plot_true_compositions <- function(
    counts,
    metadata           = NULL,
    sample_id_col      = "sample_id",
    facet_var          = NULL,
    facet_row_var      = NULL,
    sort_var           = NULL,
    top_n              = 20,
    taxa_include       = NULL,
    base_size          = 11,
    title              = NULL,
    subtitle           = NULL,
    show_legend        = FALSE,
    legend_position    = "bottom",
    legend_nrow        = NULL,
    legend_ncol        = NULL,
    legend_key_size    = 0.4,
    legend_text_size   = NULL,
    legend_rel_height  = 0.1,
    legend_rel_width   = 0.2,
    n_hues             = 7,
    bar_width          = 0.9,
    ylab               = "Proportion",
    strip_text_size    = NULL,
    row_label_size     = NULL,
    ylab_size          = NULL,
    ylab_rel_width     = 0.03,
    row_label_fn       = as.character,
    row_label_position = c("title", "ylab"),
    panel_labels       = TRUE,
    vline_var        = NULL,
    vline_value      = NULL,
    vline_color      = "black",
    vline_linetype   = "dashed",
    vline_linewidth  = 0.7
) {
  row_label_position <- match.arg(row_label_position)
  
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  N          <- nrow(counts)
  S          <- ncol(counts)
  taxa_names <- if (!is.null(colnames(counts))) colnames(counts) else paste0("Sp_", seq_len(S))
  sample_ids <- if (!is.null(rownames(counts))) rownames(counts) else paste0("S", seq_len(N))
  
  row_tots <- rowSums(counts)
  row_tots[row_tots == 0] <- 1
  freq_mat <- counts / row_tots
  
  # An explicit set wins over the top_n ranking. The palette is a deterministic
  # function of the sorted taxon set, so passing the same set that
  # eDNA_dmm_compositions() used is what makes the two figures share colours.
  if (!is.null(taxa_include)) {
    missing_taxa <- setdiff(taxa_include, taxa_names)
    if (length(missing_taxa) > 0)
      rlang::warn(paste0("`taxa_include` names not present in `counts`, dropped: ",
                         paste(missing_taxa, collapse = ", ")))
    top_taxa <- intersect(taxa_include, taxa_names)
    if (length(top_taxa) == 0)
      rlang::abort("`taxa_include` matched none of the columns of `counts`.")
  } else {
    mean_freq <- colMeans(freq_mat)
    top_taxa  <- names(sort(mean_freq, decreasing = TRUE))[seq_len(min(top_n, S))]
  }
  other_taxa <- setdiff(taxa_names, top_taxa)
  
  plot_df <- as.data.frame(freq_mat)
  colnames(plot_df) <- taxa_names
  plot_df$sample_id <- sample_ids
  
  if (length(other_taxa) > 0) {
    plot_df$Other <- rowSums(plot_df[, other_taxa, drop = FALSE])
    plot_df       <- plot_df[, c("sample_id", top_taxa, "Other")]
  } else {
    plot_df <- plot_df[, c("sample_id", top_taxa)]
  }
  
  si <- plot_df
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata))
      rlang::abort("`metadata` must be a data frame.")
    metadata <- as.data.frame(metadata)
    metadata[[".edna_sample_id"]] <- resolve_metadata_ids(
      metadata, sample_id_col, plot_df$sample_id
    )
    # Drop columns that already exist in `si` so the merge cannot emit .x/.y pairs
    dup <- setdiff(intersect(names(metadata), names(si)), ".edna_sample_id")
    if (length(dup))
      metadata <- metadata[, setdiff(names(metadata), dup), drop = FALSE]
    si <- merge(si, metadata, by.x = "sample_id", by.y = ".edna_sample_id",
                all.x = TRUE, sort = FALSE)
  }
  
  named_taxa <- sort(top_taxa)
  n_named    <- length(named_taxa)
  tax_colors <- make_taxa_colors(named_taxa, n_hues = n_hues,
                                 include_other = length(other_taxa) > 0)
  
  value_cols   <- c(top_taxa, if (length(other_taxa) > 0) "Other")
  taxon_levels <- c(sort(top_taxa), if (length(other_taxa) > 0) "Other")
  
  title_str    <- title    %||% "Observed species composition"
  subtitle_str <- subtitle %||% ""
  
  build_panel <- function(dat, row_label = NULL, show_legend = FALSE,
                          show_title = FALSE, show_ylab = TRUE) {
    if (!is.null(sort_var) && sort_var %in% names(dat)) {
      dat <- dat[order(dat[[sort_var]]), ]
    }
    dat$x_label <- factor(dat$sample_id, levels = unique(dat$sample_id))
    
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
            # x range: which silently destroys the proportional panel widths
            # that space = "free_x" is supposed to give.
            cape_x = match(
              .data$x_label[which.min(abs(.data[[vline_var]] - vline_value))],
              sort(unique(.data$x_label))
            ) + 0.5,
            .groups = "drop"
          )
      } else {
        closest <- which.min(abs(dat[[vline_var]] - vline_value))
        vline_df <- data.frame(cape_x = as.integer(dat$x_label[closest]) + 0.5)
      }
    }
    
    plot_long <- tidyr::pivot_longer(dat, cols = dplyr::all_of(value_cols),
                                     names_to = "taxon", values_to = "frequency")
    plot_long$taxon <- factor(plot_long$taxon, levels = taxon_levels)
    
    panel_title <- if (!is.null(row_label) && row_label_position == "title") row_label else NULL
    # In "ylab" mode the row label IS the y-axis title. In "title" mode the axis
    # title would otherwise repeat on every stacked row panel, so callers ask for
    # it on one panel only.
    panel_ylab  <- if (!is.null(row_label) && row_label_position == "ylab") {
      row_label
    } else if (show_ylab) {
      "Proportion"
    } else {
      NULL
    }
    
    p <- ggplot2::ggplot(
      plot_long,
      ggplot2::aes(x = .data$x_label, y = .data$frequency, fill = .data$taxon)
    ) +
      # bar_width controls the gutter between samples. A value below 1 leaves
      # 10% gap between every pair of samples, which reads as white striping.
      ggplot2::geom_bar(stat = "identity", position = "stack", color = NA,
                        width = bar_width) +
      ggplot2::scale_fill_manual(values = tax_colors, name = NULL, drop = FALSE) +
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
        axis.text.x      = ggplot2::element_blank(),
        axis.ticks.x     = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text       = ggplot2::element_text(face = "bold", size = strip_text_size %||% base_size),
        panel.spacing.x  = ggplot2::unit(0.3, "lines"),
        plot.title       = ggplot2::element_text(face = "bold", size = row_label_size %||% base_size, hjust = 0),
        plot.subtitle    = ggplot2::element_text(color = "grey40", size = base_size - 2)
      )
    
    if (!is.null(vline_df)) {
      p <- p + ggplot2::geom_vline(
        data        = vline_df,
        ggplot2::aes(xintercept = .data$cape_x),
        color       = vline_color,
        linetype    = vline_linetype,
        linewidth   = vline_linewidth,
        inherit.aes = FALSE
      )
    }
    
    if (!is.null(facet_var) && facet_var %in% names(dat)) {
      p <- p + ggplot2::facet_grid(reformulate(facet_var), scales = "free_x", space = "free_x")
    }
    
    p
  }
  
  if (is.null(facet_row_var)) {
    p <- build_panel(si, show_legend = show_legend, show_title = TRUE)
    p <- p + ggplot2::labs(title = title_str, subtitle = subtitle_str)
    return(p)
  }
  
  if (!facet_row_var %in% names(si))
    rlang::abort(paste0("`facet_row_var = '", facet_row_var, "'` not found."))
  
  row_levels <- if (is.factor(si[[facet_row_var]])) {
    levels(si[[facet_row_var]])
  } else {
    sort(unique(si[[facet_row_var]]))
  }
  
  # In "ylab" mode each panel's y-axis title IS its row label. In "title" mode
  # the row label moves to the panel title, so no panel carries an axis title: 
  # a shared one is drawn for the whole stack below.
  panels <- lapply(seq_along(row_levels), function(i) {
    lv  <- row_levels[i]
    dat <- si[si[[facet_row_var]] == lv, ]
    if (nrow(dat) == 0) return(NULL)
    build_panel(dat, row_label = row_label_fn(lv), show_legend = FALSE,
                show_title = FALSE,
                show_ylab = row_label_position == "ylab")
  })
  panels <- Filter(Negate(is.null), panels)
  
  # With both title and subtitle empty: the manuscript case, where that text
  # belongs in the figure legend: the heading row is dropped entirely rather
  # than reserved as blank space above the panels.
  has_heading <- nzchar(title_str) || nzchar(subtitle_str)
  title_grob <- if (has_heading) {
    cowplot::ggdraw() +
      cowplot::draw_label(title_str, fontface = "bold", size = base_size + 1, x = 0.02, hjust = 0) +
      cowplot::draw_label(subtitle_str, size = base_size - 1, color = "grey40", x = 0.02, y = 0.25, hjust = 0)
  } else NULL

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

  # Assemble rows top to bottom, prepending the heading only when there is one.
  with_title <- function(parts, heights) {
    if (is.null(title_grob))
      return(cowplot::plot_grid(plotlist = parts, ncol = 1, rel_heights = heights))
    cowplot::plot_grid(plotlist = c(list(title_grob), parts), ncol = 1,
                       rel_heights = c(0.06, heights))
  }

  if (show_legend) {
    legend_plot <- build_panel(si, show_legend = TRUE) +
      ggplot2::theme(legend.position = legend_position)

    if (legend_position %in% c("right", "left")) {
      side <- paste0("guide-box-", legend_position)
      legend_grob <- cowplot::get_plot_component(legend_plot, side, return_all = TRUE)
      inner <- if (legend_position == "right")
        cowplot::plot_grid(stacked, legend_grob, nrow = 1, rel_widths = c(1, legend_rel_width))
      else
        cowplot::plot_grid(legend_grob, stacked, nrow = 1, rel_widths = c(legend_rel_width, 1))
      with_title(list(inner), 1)
    } else {
      guide_box_name <- if (legend_position == "top") "guide-box-top" else "guide-box-bottom"
      legend_grob <- cowplot::get_plot_component(legend_plot, guide_box_name, return_all = TRUE)
      # legend_rel_height is a fraction of the panel stack, so a legend with
      # many rows needs it raised or cowplot clips the bottom rows silently.
      if (legend_position == "top")
        with_title(list(legend_grob, stacked), c(legend_rel_height, 1))
      else
        with_title(list(stacked, legend_grob), c(1, legend_rel_height))
    }
  } else {
    with_title(list(stacked), 1)
  }
}
