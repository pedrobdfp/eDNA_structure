#' eDNA_structure: Dirichlet-Multinomial Mixture Models for eDNA Metabarcoding
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

.onLoad <- function(libname, pkgname) {
  rstan::rstan_options(auto_write = TRUE)

  modules <- paste0("stan_fit4", names(stanmodels), "_mod")

  for (m in modules) {
    tryCatch(
      Rcpp::loadModule(m, what = FALSE),
      error = function(e) {
        packageStartupMessage(
          "Note: Stan module '", m, "' could not be loaded automatically (",
          conditionMessage(e), "). ",
          "Model-fitting functions in this package may not work; ",
          "all other functions are unaffected."
        )
      }
    )
  }
}
