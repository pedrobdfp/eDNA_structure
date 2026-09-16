# =============================================================================
# eDNA_dmm_compositions(): posterior community composition bar plots
# dmm_taxon_order(): shared taxon ranking driving both composition plots
# =============================================================================

#' Rank taxa by their weight in a fitted model
#'
#' @description
#' Returns taxa ordered by total posterior weight across communities,
#' `colSums(pi_mean)`. A taxon that dominates a single community ranks highly
#' even if it is rare overall, which is the ordering that matters when the
#' question is what distinguishes communities from one another.
#'
#' This is deliberately NOT the same ranking as mean observed frequency across
#' samples, which rewards being widespread instead. Use this function to drive
#' BOTH [eDNA_dmm_compositions()] and [plot_true_compositions()] when the two
#' figures sit side by side, because the taxon palette is a deterministic
#' function of the sorted taxon set: same set in, same colours out. Selecting
#' the two sets independently gives two different palettes, and shared taxa
#' then appear in different colours in the two panels.
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param top_n Number of taxa to return. `NULL` (default) returns all of them,
#'   ranked.
#' @return A character vector of taxon names, highest weight first.
#' @seealso [eDNA_dmm_compositions()], [plot_true_compositions()]
#' @export
dmm_taxon_order <- function(fit, top_n = NULL) {
  check_stan_fit_object(fit, "dmm_taxon_order")
  pi_mean    <- fit$pi_mean
  taxa_names <- colnames(pi_mean) %||% paste0("Sp_", seq_len(ncol(pi_mean)))
  colnames(pi_mean) <- taxa_names
  ranked <- names(sort(colSums(pi_mean), decreasing = TRUE))
  if (is.null(top_n)) ranked else ranked[seq_len(min(top_n, length(ranked)))]
}

#' Plot posterior mean species composition for each community
#'
#' @description
#' Produces a stacked bar plot of the **posterior mean species composition**
#' for each fitted community (the pi matrix). One bar per community, colored
#' by species using the same structured palette as [plot_true_compositions()],
#' so species colors are directly comparable between the two plots.
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param top_n Number of most abundant taxa to show individually; the rest
#'   are collapsed into `"Other"`. Default `20`.
#' @param base_size Base font size. Default `13`.
#' @param title Plot title. Default auto-generated.
#' @param subtitle Plot subtitle. Default `NULL`, meaning no subtitle.
#' @param legend_position Legend position. Default `"right"`.
#' @param bar_width Bar width. Default `0.7`.
#'
#' @return A [ggplot2::ggplot()] object.
#' @seealso [plot_true_compositions()], [eDNA_dmm_structure()]
#' @export
eDNA_dmm_compositions <- function(
    fit,
    top_n            = 20,
    base_size        = 13,
    title            = NULL,
    subtitle         = NULL,
    legend_position  = "right",
    bar_width        = 0.7
) {
  check_stan_fit_object(fit, "eDNA_dmm_compositions")
  
  K          <- fit$K
  pi_mean    <- fit$pi_mean          # K x S matrix, posterior mean compositions
  taxa_names <- colnames(pi_mean)
  
  if (is.null(taxa_names))
    taxa_names <- paste0("Sp_", seq_len(ncol(pi_mean)))
  
  # ── Which taxa to show ─────────────────────────────────────────────────────
  # Ranked by total posterior weight across communities. To make a companion
  # plot_true_compositions() figure use the SAME colours, pass this same set to
  # its `taxa_include` argument: see dmm_taxon_order().
  top_taxa   <- dmm_taxon_order(fit, top_n)
  other_taxa <- setdiff(taxa_names, top_taxa)
  
  # ── Structured color palette: identical to plot_true_compositions() ───────
  named_taxa <- sort(top_taxa)
  n_named    <- length(named_taxa)
  tax_colors <- make_taxa_colors(named_taxa,
                                 include_other = length(other_taxa) > 0)
  
  taxon_levels <- c(sort(top_taxa), if (length(other_taxa) > 0) "Other")
  
  # ── Build long data frame ──────────────────────────────────────────────────
  comp_df <- as.data.frame(pi_mean)
  colnames(comp_df) <- taxa_names
  comp_df$Community <- as.character(seq_len(K))
  
  if (length(other_taxa) > 0) {
    comp_df$Other <- rowSums(comp_df[, other_taxa, drop = FALSE])
  }
  
  plot_long <- tidyr::pivot_longer(
    comp_df,
    cols      = dplyr::all_of(c(top_taxa, if (length(other_taxa) > 0) "Other")),
    names_to  = "taxon",
    values_to = "proportion"
  )
  plot_long$taxon     <- factor(plot_long$taxon, levels = taxon_levels)
  plot_long$Community <- factor(plot_long$Community, levels = as.character(seq_len(K)))
  
  # ── Titles ─────────────────────────────────────────────────────────────────
  title_str    <- title    %||% sprintf("Posterior community compositions  (K = %d)", K)
  subtitle_str <- subtitle
  
  # ── Plot ───────────────────────────────────────────────────────────────────
  ggplot2::ggplot(
    plot_long,
    ggplot2::aes(x = .data$Community, y = .data$proportion, fill = .data$taxon)
  ) +
    ggplot2::geom_bar(stat = "identity", width = bar_width, color = NA) +
    ggplot2::scale_fill_manual(values = tax_colors, name = "Taxon", drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(
      x        = "Community",
      y        = "Posterior mean proportion",
      title    = title_str,
      subtitle = subtitle_str
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position  = legend_position,
      strip.background = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = base_size - 2)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend())
}
