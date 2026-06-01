#' eDNA_structure: Dirichlet-Multinomial Mixture Models for eDNA Metabarcoding
#'
#' @description
#' `eDNA_structure` fits Bayesian Dirichlet-Multinomial Mixture (DMM) models
#' to eDNA read count data from metabarcoding surveys. The model simultaneously
#' estimates:
#'
#' - How many latent ecological communities are present
#' - The taxonomic composition (relative frequencies) of each community
#' - How environmental covariates (depth, latitude, etc.) drive community membership
#'
#' Input can be raw ASV tables, taxonomically-annotated taxon tables, or
#' any sample × feature count matrix.
#'
#' @section Core workflow:
#' 1. **Fit the model**: [eDNA_dmm()]
#' 2. **Visualize assignments**: [eDNA_dmm_structure()]
#' 3. **Ordinate communities**: [eDNA_dmm_nmds()]
#' 4. **Inspect covariate effects**: [eDNA_dmm_beta()]
#' 5. **Get example data**: [get_example_data()]
#'
#' @section Getting started:
#' ```r
#' # Load example data and run a quick K=2 fit
#' data <- get_example_data()
#' fit  <- eDNA_dmm(counts = data$counts, covariates = data$covariates, K = 2)
#' eDNA_dmm_structure(fit)
#' ```
#'
#' @section Package philosophy:
#' This package is designed to be approachable for users unfamiliar with
#' Bayesian mixture models. All functions expose every relevant argument
#' explicitly (no `...` pass-through), include informative error messages
#' for common mistakes, and produce [ggplot2][ggplot2::ggplot2-package]
#' objects that can be freely customized downstream.
#'
#' @docType package
#' @name eDNA_structure-package
#' @aliases eDNA_structure
"_PACKAGE"

#' @importFrom rstan stan extract as.matrix.stanfit
#' @importFrom ggplot2 ggplot aes geom_bar geom_point geom_line geom_hline
#'   geom_vline geom_density geom_linerange stat_ellipse
#'   scale_fill_manual scale_color_manual scale_x_continuous scale_y_continuous
#'   scale_size_continuous facet_wrap facet_grid labs theme theme_bw
#'   element_text element_blank element_rect margin guides guide_legend
#'   position_dodge expansion coord_fixed
#' @importFrom dplyr select filter mutate arrange group_by summarise ungroup
#'   bind_cols bind_rows left_join distinct pull n across starts_with
#'   rename slice_max
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom vegan vegdist metaMDS scores
#' @importFrom posterior ess_bulk as_draws_matrix
#' @importFrom loo loo extract_log_lik
#' @importFrom scales percent percent_format pretty_breaks
#' @importFrom rlang .data abort warn inform
#' @importFrom methods is
NULL
