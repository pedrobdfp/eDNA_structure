# =============================================================================
# Lazy, cached Stan model compilation
# =============================================================================
# Stan models are compiled on first use (not at package install time), then
# cached on disk so subsequent R sessions reuse the compiled model instantly.
# This avoids the Rcpp Module / RCPP_MODULE Windows DLL export issue entirely:
# there is no precompiled module baked into the package's compiled code.

.eDNA_stanmodels_cache <- new.env(parent = emptyenv())

#' @keywords internal
.get_dmm_stanmodel <- function() {
  if (is.null(.eDNA_stanmodels_cache$dmm)) {

    cache_dir  <- tools::R_user_dir("eDNAstructure", which = "cache")
    cache_path <- file.path(cache_dir, "dmm_model.rds")

    if (file.exists(cache_path)) {
      model <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    } else {
      model <- NULL
    }

    if (is.null(model)) {
      stan_file <- system.file("stan", "dmm.stan", package = "eDNAstructure")
      if (!file.exists(stan_file)) {
        rlang::abort("Could not find dmm.stan in the installed package.")
      }

      message("Compiling Stan model (first use on this machine, ~30-60s)...")
      message("This will be cached and instant on future calls.")

      model <- rstan::stan_model(file = stan_file, model_name = "dmm")

      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      tryCatch(
        saveRDS(model, cache_path),
        error = function(e) {
          packageStartupMessage(
            "Note: could not cache compiled Stan model to disk (",
            conditionMessage(e), "). Model will be recompiled next session."
          )
        }
      )
    }

    .eDNA_stanmodels_cache$dmm <- model
  }

  .eDNA_stanmodels_cache$dmm
}

#' Clear the cached compiled Stan model
#'
#' @description
#' Forces recompilation of the Stan model on next use. Useful after updating
#' the package, or if the cached model becomes stale/corrupted.
#'
#' @export
eDNA_clear_stan_cache <- function() {
  cache_dir  <- tools::R_user_dir("eDNAstructure", which = "cache")
  cache_path <- file.path(cache_dir, "dmm_model.rds")
  if (file.exists(cache_path)) {
    file.remove(cache_path)
    message("Cleared cached Stan model: ", cache_path)
  } else {
    message("No cached Stan model found.")
  }
  rm(list = ls(.eDNA_stanmodels_cache), envir = .eDNA_stanmodels_cache)
  invisible(NULL)
}
