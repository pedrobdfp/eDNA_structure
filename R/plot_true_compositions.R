# =============================================================================
# plot_true_compositions() — visualize observed species composition per sample
# =============================================================================

#' Plot observed species composition per sample (pre-fit diagnostic)
#'
#' @description
#' Produces a stacked bar plot of **observed species frequencies** — one bar
#' per sample, colored by species. This is a pre-fitting diagnostic that shows
#' the raw community signal in your data before any model is applied. It uses
#' the same visual format as [eDNA_dmm_structure()] so the two plots can be
#' compared directly: the top panel shows what species are present, the bottom
#' shows what community they were assigned to.
#'
#' This function is most useful with **simulated data** (where the true
#' community of each sample is known) or with any dataset where you want to
#' visually inspect raw compositional structure.
#'
#' @param counts A numeric matrix or data frame of read counts (samples × taxa),
#'   the same object passed to [eDNA_dmm()].
#' @param metadata An optional data frame with one row per sample for faceting.
#'   Pass `NULL` (default) for a single panel.
#' @param sample_id_col A string: column in `metadata` matching `rownames(counts)`.
#'   Default `"sample_id"`.
#' @param facet_var A string: column in `metadata` to facet by (e.g.,
#'   `"TrueCommunity"` for simulated data). Default `NULL`.
#' @param sort_var A string: column in `metadata` to sort samples within panels.
#'   Default `NULL`.
#' @param top_n An integer: number of most abundant taxa to show individually;
#'   the rest are collapsed into `"Other"`. Default `20`.
#' @param bar_width A number in (0, 1]: bar width. Default `0.95`.
#' @param base_size A positive number: base font size. Default `11`.
#' @param title A string: plot title. Default auto-generated.
#' @param subtitle A string: plot subtitle. Default auto-generated.
#' @param legend_position A string: `"none"` (default — there are too many
#'   species for a useful legend), `"right"`, `"bottom"`, etc.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @seealso [eDNA_dmm_structure()], [simulate_eDNA_survey()]
#'
#' @examples
#' \dontrun{
#' sim <- simulate_eDNA_survey(n_communities = 3, n_species = 20, seed = 1)
#'
#' # All samples, no faceting
#' plot_true_compositions(sim$counts)
#'
#' # Faceted by true community (only possible when you know true labels)
#' plot_true_compositions(
#'   sim$counts,
#'   metadata  = sim$covariates,
#'   facet_var = "TrueCommunity",
#'   sort_var  = "Depth"
#' )
#' }
#'
#' @export
plot_true_compositions <- function(
    counts,
    metadata         = NULL,
    sample_id_col    = "sample_id",
    facet_var        = NULL,
    sort_var         = NULL,
    top_n            = 20,
    bar_width        = 0.95,
    base_size        = 11,
    title            = NULL,
    subtitle         = NULL,
    legend_position  = "none"
) {
  # ── Coerce to matrix ───────────────────────────────────────────────────────
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  N <- nrow(counts)
  S <- ncol(counts)

  taxa_names  <- if (!is.null(colnames(counts))) colnames(counts) else paste0("Sp_", seq_len(S))
  sample_ids  <- if (!is.null(rownames(counts))) rownames(counts) else paste0("S", seq_len(N))

  # ── Relative frequencies ──────────────────────────────────────────────────
  row_tots <- rowSums(counts)
  row_tots[row_tots == 0] <- 1
  freq_mat <- counts / row_tots

  # ── Top N taxa ─────────────────────────────────────────────────────────────
  mean_freq   <- colMeans(freq_mat)
  top_taxa    <- names(sort(mean_freq, decreasing = TRUE))[seq_len(min(top_n, S))]
  other_taxa  <- setdiff(taxa_names, top_taxa)

  # ── Long data frame ────────────────────────────────────────────────────────
  plot_df <- as.data.frame(freq_mat)
  colnames(plot_df) <- taxa_names
  plot_df$sample_id <- sample_ids

  if (length(other_taxa) > 0) {
    plot_df$Other <- rowSums(plot_df[, other_taxa, drop = FALSE])
    plot_df       <- plot_df[, c("sample_id", top_taxa, "Other")]
  } else {
    plot_df <- plot_df[, c("sample_id", top_taxa)]
  }

  # ── Merge metadata ─────────────────────────────────────────────────────────
  si <- plot_df
  if (!is.null(metadata)) {
    if (!sample_id_col %in% names(metadata)) {
      rlang::abort(paste0("Column '", sample_id_col, "' not found in `metadata`. ",
                          "Available: ", paste(names(metadata), collapse = ", ")))
    }
    si <- merge(si, metadata, by.x = "sample_id", by.y = sample_id_col, all.x = TRUE)
  }

  # ── Sort and x positions ──────────────────────────────────────────────────
  if (!is.null(sort_var) && sort_var %in% names(si)) {
    si <- si[order(si[[sort_var]]), ]
  }

  if (!is.null(facet_var) && facet_var %in% names(si)) {
    si <- si |>
      dplyr::group_by(dplyr::across(dplyr::all_of(facet_var))) |>
      dplyr::mutate(x_pos = dplyr::row_number()) |>
      dplyr::ungroup()
  } else {
    si$x_pos <- seq_len(nrow(si))
  }

  # ── Pivot long ─────────────────────────────────────────────────────────────
  value_cols <- c(top_taxa, if (length(other_taxa) > 0) "Other")
  plot_long  <- tidyr::pivot_longer(
    si,
    cols      = dplyr::all_of(value_cols),
    names_to  = "taxon",
    values_to = "frequency"
  )
  all_taxa_ordered <- c(top_taxa, if (length(other_taxa) > 0) "Other")
  plot_long$taxon  <- factor(plot_long$taxon, levels = all_taxa_ordered)

  # ── Colors ─────────────────────────────────────────────────────────────────
  n_colors   <- length(top_taxa)
  hues       <- seq(15, 375, length.out = n_colors + 1)[seq_len(n_colors)]
  tax_colors <- c(
    stats::setNames(grDevices::hcl(h = hues, c = 65, l = 55), top_taxa),
    if (length(other_taxa) > 0) c(Other = "grey75") else NULL
  )

  # ── Titles ─────────────────────────────────────────────────────────────────
  title_str    <- title    %||% "Observed species composition"
  subtitle_str <- subtitle %||% sprintf(
    "%d samples  |  top %d taxa shown individually", N, min(top_n, S)
  )

  # ── Plot ───────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    plot_long,
    ggplot2::aes(x = .data$x_pos, y = .data$frequency, fill = .data$taxon)
  ) +
    ggplot2::geom_bar(stat = "identity", width = bar_width, position = "stack",
                      color = NA) +
    ggplot2::scale_fill_manual(values = tax_colors, name = NULL) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0),
                                breaks = c(0, 0.5, 1)) +
    ggplot2::labs(
      x        = NULL,
      y        = "Species frequency",
      title    = title_str,
      subtitle = subtitle_str
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position  = legend_position,
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.spacing.x  = ggplot2::unit(0.3, "lines"),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = base_size - 2)
    )

  if (!is.null(facet_var) && facet_var %in% names(si)) {
    p <- p + ggplot2::facet_grid(
      reformulate(facet_var),
      scales = "free_x",
      space  = "free_x"
    )
  }

  p
}
