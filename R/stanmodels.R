# =============================================================================
# Lazy, cached Stan model compilation
# =============================================================================
# Stan models are compiled on first use (not at package install time), then
# cached on disk so subsequent R sessions reuse the compiled model instantly.
# This avoids the Rcpp Module / RCPP_MODULE Windows DLL export issue entirely:
# there is no precompiled module baked into the package's compiled code.
#
# The package ships one Stan program, "dmm" (inst/stan/dmm.stan). It is
# compiled the first time it is used and cached thereafter.

.eDNA_stanmodels_cache <- new.env(parent = emptyenv())

# Stan programs shipped in inst/stan, without the .stan extension.
.eDNA_stan_model_names <- c("dmm")

#' @keywords internal
.stanmodel_cache_path <- function(name) {
  file.path(tools::R_user_dir("eDNAstructure", which = "cache"),
            paste0(name, "_model.rds"))
}

# Sidecar file recording the md5 of the .stan source the cached .rds was
# compiled from. Compared against the currently-installed source on every
# load, so a package update (which can rewrite the Stan code, e.g. changing
# the beta parameterization) invalidates a stale cached binary instead of it
# being silently reused forever.
#' @keywords internal
.stanmodel_hash_path <- function(name) {
  file.path(tools::R_user_dir("eDNAstructure", which = "cache"),
            paste0(name, "_model.hash"))
}

#' Internal accessor for a compiled Stan model
#'
#' Returns the compiled `stanmodel` for `name`, loading it from the in-session
#' cache, then the per-machine disk cache (only if it still matches the
#' installed .stan source), and compiling from source otherwise.
#'
#' @param name Model name, one of `.eDNA_stan_model_names`.
#' @return A [rstan::stanmodel-class] object.
#' @keywords internal
.get_stanmodel <- function(name = "dmm") {
  name <- match.arg(name, .eDNA_stan_model_names)

  if (!is.null(.eDNA_stanmodels_cache[[name]])) {
    return(.eDNA_stanmodels_cache[[name]])
  }

  stan_file <- system.file("stan", paste0(name, ".stan"),
                           package = "eDNAstructure")
  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    rlang::abort(paste0("Could not find ", name, ".stan in the installed package."))
  }
  current_hash <- unname(tools::md5sum(stan_file))

  cache_path <- .stanmodel_cache_path(name)
  hash_path  <- .stanmodel_hash_path(name)

  cache_is_current <- file.exists(cache_path) && file.exists(hash_path) &&
    isTRUE(tryCatch(readLines(hash_path, n = 1L) == current_hash,
                     error = function(e) FALSE))

  model <- if (cache_is_current) {
    tryCatch(readRDS(cache_path), error = function(e) NULL)
  } else {
    NULL
  }

  if (is.null(model)) {
    if (file.exists(cache_path) && !cache_is_current) {
      message("Stan source for '", name,
              "' has changed since it was last compiled; recompiling.")
    }

    message("Compiling Stan model '", name,
            "' (first use on this machine, ~30-60s)...")
    message("This will be cached and instant on future calls.")

    model <- rstan::stan_model(file = stan_file, model_name = name)

    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    tryCatch({
      saveRDS(model, cache_path)
      writeLines(current_hash, hash_path)
    }, error = function(e) {
      packageStartupMessage(
        "Note: could not cache compiled Stan model to disk (",
        conditionMessage(e), "). Model will be recompiled next session."
      )
    })
  }

  .eDNA_stanmodels_cache[[name]] <- model
  model
}

# Back-compatible wrapper for the standard model.
#' @keywords internal
.get_dmm_stanmodel <- function() {
  .get_stanmodel("dmm")
}

#' Clear the cached compiled Stan models
#'
#' @description
#' Forces recompilation of the Stan model on next use. Useful after updating
#' the package, or if a cached model becomes stale or corrupted.
#'
#' @return Invisibly `NULL`, called for its side effect.
#'
#' @export
eDNA_clear_stan_cache <- function() {
  found <- FALSE
  for (name in .eDNA_stan_model_names) {
    cache_path <- .stanmodel_cache_path(name)
    hash_path  <- .stanmodel_hash_path(name)
    if (file.exists(cache_path)) {
      file.remove(cache_path)
      message("Cleared cached Stan model: ", cache_path)
      found <- TRUE
    }
    if (file.exists(hash_path)) file.remove(hash_path)
  }
  if (!found) message("No cached Stan models found.")
  rm(list = ls(.eDNA_stanmodels_cache), envir = .eDNA_stanmodels_cache)
  invisible(NULL)
}
