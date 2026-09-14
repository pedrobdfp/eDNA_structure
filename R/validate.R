# =============================================================================
# Input validation helpers
# All user-facing checks live here. The philosophy: catch mistakes early,
# explain what went wrong AND how to fix it.
# =============================================================================

#' @keywords internal
validate_counts <- function(counts, call = rlang::caller_env()) {
  # ── Type check ──────────────────────────────────────────────────────────────
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    rlang::abort(
      c(
        "`counts` must be a numeric matrix or data frame.",
        i = paste0("You supplied an object of class: ", paste(class(counts), collapse = ", ")),
        i = "Expected format: rows = samples, columns = taxa (or ASVs).",
        i = "Example: a matrix where counts[i, j] = reads of taxon j in sample i."
      ),
      call = call
    )
  }

  if (is.data.frame(counts)) {
    # Check for non-numeric columns (e.g. a stray Sample ID column)
    non_num <- names(counts)[!vapply(counts, is.numeric, logical(1))]
    if (length(non_num) > 0) {
      rlang::abort(
        c(
          "`counts` data frame contains non-numeric columns.",
          i = paste0("Non-numeric columns found: ", paste(non_num, collapse = ", ")),
          i = "All columns must be numeric read counts.",
          i = "If one column is a sample ID, move it to `rownames(counts)` first.",
          i = 'Example: `rownames(counts) <- counts$SampleID; counts$SampleID <- NULL`'
        ),
        call = call
      )
    }
    counts <- as.matrix(counts)
  }

  # ── Dimension check ─────────────────────────────────────────────────────────
  if (nrow(counts) < 3) {
    rlang::abort(
      c(
        paste0("`counts` has only ", nrow(counts), " row(s); at least 3 samples are required."),
        i = "The DMM cannot meaningfully partition community structure with fewer than 3 samples.",
        i = "Each row of `counts` should be one sample (station, replicate, etc.)."
      ),
      call = call
    )
  }

  if (ncol(counts) < 2) {
    rlang::abort(
      c(
        paste0("`counts` has only ", ncol(counts), " column(s); at least 2 taxa are required."),
        i = "Each column of `counts` should be one taxon or ASV."
      ),
      call = call
    )
  }

  # ── Value checks ─────────────────────────────────────────────────────────────
  if (any(counts < 0, na.rm = TRUE)) {
    rlang::abort(
      c(
        "`counts` contains negative values.",
        i = "Read counts must be zero or positive integers.",
        i = paste0("Negative values found at: ",
                   paste(which(counts < 0, arr.ind = TRUE)[1:min(3, sum(counts < 0)), ,
                               drop = FALSE] |>
                           apply(1, function(r) paste0("[", r[1], ",", r[2], "]")),
                         collapse = ", "))
      ),
      call = call
    )
  }

  if (any(is.na(counts))) {
    n_na <- sum(is.na(counts))
    rlang::abort(
      c(
        paste0("`counts` contains ", n_na, " NA value(s)."),
        i = "All count values must be non-missing.",
        i = "Replace NAs with 0 if a taxon was simply not detected: `counts[is.na(counts)] <- 0`"
      ),
      call = call
    )
  }

  # ── Zero-row samples ─────────────────────────────────────────────────────────
  row_sums <- rowSums(counts)
  empty_samples <- which(row_sums == 0)
  if (length(empty_samples) > 0) {
    sample_names <- if (!is.null(rownames(counts))) rownames(counts)[empty_samples] else
      paste0("row ", empty_samples)
    rlang::abort(
      c(
        paste0(length(empty_samples), " sample(s) have zero total reads and cannot be modelled."),
        i = paste0("Empty sample(s): ", paste(sample_names[seq_len(min(5, length(sample_names)))],
                                              collapse = ", ")),
        i = "Remove these rows before fitting: `counts <- counts[rowSums(counts) > 0, ]`"
      ),
      call = call
    )
  }

  # ── Warn about suspiciously low read counts ───────────────────────────────────
  low_samples <- which(row_sums < 100)
  if (length(low_samples) > 0) {
    rlang::warn(
      c(
        paste0(length(low_samples), " sample(s) have fewer than 100 total reads."),
        i = "Very low read counts produce unreliable composition estimates.",
        i = "Consider filtering: `counts <- counts[rowSums(counts) >= 100, ]`"
      )
    )
  }

  # ── All-zero taxa ─────────────────────────────────────────────────────────────
  col_sums <- colSums(counts)
  zero_taxa <- which(col_sums == 0)
  if (length(zero_taxa) > 0) {
    taxa_names <- if (!is.null(colnames(counts))) colnames(counts)[zero_taxa] else
      paste0("column ", zero_taxa)
    rlang::warn(
      c(
        paste0(length(zero_taxa), " taxon/taxa column(s) have zero reads across all samples."),
        i = paste0("Zero-count column(s): ",
                   paste(taxa_names[seq_len(min(5, length(taxa_names)))], collapse = ", ")),
        i = "These are dropped automatically before fitting.",
        i = "To suppress this warning, remove them yourself: `counts <- counts[, colSums(counts) > 0]`"
      )
    )
  }

  # ── Storage mode ──────────────────────────────────────────────────────────────
  if (!is.integer(counts)) {
    if (!all(counts == floor(counts), na.rm = TRUE)) {
      rlang::warn(
        c(
          "`counts` contains non-integer values. These will be coerced to integers via `round()`.",
          i = "Stan requires integer count data.",
          i = "If your data are already integers stored as doubles, this is expected and safe.",
          i = "If your data are normalized read proportions, you need raw counts."
        )
      )
    }
    storage.mode(counts) <- "integer"
  }

  counts
}

#' @keywords internal
validate_covariates <- function(covariates, counts, scale_covariates,
                                n_expected = NULL, call = rlang::caller_env()) {
  # `n_expected` overrides nrow(counts) for the replicate-aware model, where
  # covariates are station-level but `counts` has one row per replicate.
  N <- n_expected %||% nrow(counts)

  # ── NULL means intercept-only model ──────────────────────────────────────────
  if (is.null(covariates)) {
    return(matrix(numeric(0), nrow = N, ncol = 0))
  }

  # ── Type ──────────────────────────────────────────────────────────────────────
  if (!is.matrix(covariates) && !is.data.frame(covariates)) {
    rlang::abort(
      c(
        "`covariates` must be a numeric matrix or data frame (or NULL for an intercept-only model).",
        i = paste0("You supplied an object of class: ", paste(class(covariates), collapse = ", ")),
        i = "Expected format: rows = samples (same order as `counts`), columns = covariates.",
        i = "Example: a data frame with columns `depth` and `latitude`, one row per sample."
      ),
      call = call
    )
  }

  if (is.data.frame(covariates)) {
    non_num <- names(covariates)[!vapply(covariates, is.numeric, logical(1))]
    if (length(non_num) > 0) {
      rlang::abort(
        c(
          "`covariates` data frame contains non-numeric columns.",
          i = paste0("Non-numeric columns: ", paste(non_num, collapse = ", ")),
          i = "All covariate columns must be numeric.",
          i = "For categorical covariates, create indicator (dummy) columns manually."
        ),
        call = call
      )
    }
    covariates <- as.matrix(covariates)
  }

  # ── Dimension alignment ───────────────────────────────────────────────────────
  if (nrow(covariates) != N) {
    rlang::abort(
      c(
        paste0("`covariates` has ", nrow(covariates), " rows but ", N, " were expected."),
        i = "`covariates` and the samples being modelled must have the same number of rows, one per sample, in the same order.",
        i = "Check that you haven't filtered one object without filtering the other."
      ),
      call = call
    )
  }

  # ── Missing values ────────────────────────────────────────────────────────────
  if (any(is.na(covariates))) {
    n_na <- sum(is.na(covariates))
    na_cols <- colnames(covariates)[apply(covariates, 2, anyNA)]
    rlang::abort(
      c(
        paste0("`covariates` contains ", n_na, " NA value(s)."),
        i = paste0("Columns with NAs: ", paste(na_cols, collapse = ", ")),
        i = "Stan cannot handle missing covariate values.",
        i = "Options: (1) impute NAs, (2) drop rows with NAs from both `counts` and `covariates`,",
        i = "         (3) set `covariates = NULL` to fit an intercept-only model."
      ),
      call = call
    )
  }

  # ── Near-zero variance covariates ─────────────────────────────────────────────
  col_vars <- apply(covariates, 2, var)
  zero_var <- which(col_vars < .Machine$double.eps * 100)
  if (length(zero_var) > 0) {
    zero_names <- if (!is.null(colnames(covariates))) colnames(covariates)[zero_var] else
      paste0("column ", zero_var)
    rlang::abort(
      c(
        paste0("Covariate(s) have (near-)zero variance and cannot be used: ",
               paste(zero_names, collapse = ", ")),
        i = "A constant covariate provides no information about community membership.",
        i = "Remove it from `covariates` or check that your data are correct."
      ),
      call = call
    )
  }

  # ── Scaling ───────────────────────────────────────────────────────────────────
  if (scale_covariates) {
    scaled <- scale(covariates)
    attr(scaled, "scale_center") <- attr(scaled, "scaled:center")
    attr(scaled, "scale_scale")  <- attr(scaled, "scaled:scale")
    attr(scaled, "scaled:center") <- NULL
    attr(scaled, "scaled:scale")  <- NULL
    covariates <- scaled
  } else {
    # Warn if covariates look unscaled and the user said not to scale
    col_ranges <- apply(covariates, 2, function(x) diff(range(x)))
    big_range  <- col_ranges > 100
    if (any(big_range)) {
      rlang::warn(
        c(
          "Some covariates have a large numeric range and `scale_covariates = FALSE`.",
          i = paste0("Wide-range columns: ",
                     paste(colnames(covariates)[big_range], collapse = ", ")),
          i = "Unscaled covariates can cause slow mixing and poor convergence.",
          i = "Strongly recommended: set `scale_covariates = TRUE` (the default)."
        )
      )
    }
  }

  covariates
}

# =============================================================================
# Replicate structure (station_id)
# =============================================================================

#' Validate a station_id vector against a replicate-level count matrix
#'
#' Returns the information `eDNA_dmm()` needs to build ragged replicate data,
#' plus a `fallback` flag. `fallback = TRUE` means no station has more than one
#' replicate, in which case the replicate model's two dispersion parameters
#' (`alpha`, `phi`) are not separately identified — and the caller should simply
#' use the standard model, which is the correct model for unreplicated data.
#'
#' @return A list with `index` (integer station index per row), `levels`,
#'   `n_stations`, `reps` (replicates per station) and `fallback`.
#' @keywords internal
validate_station_id <- function(station_id, counts, call = rlang::caller_env()) {
  R <- nrow(counts)

  if (is.matrix(station_id) || is.data.frame(station_id)) {
    if (ncol(as.data.frame(station_id)) != 1) {
      rlang::abort(
        c(
          "`station_id` must be a vector, not a multi-column object.",
          i = "Supply one station label per row of `counts`, e.g. `metadata$station`."
        ),
        call = call
      )
    }
    station_id <- as.data.frame(station_id)[[1]]
  }

  if (!is.atomic(station_id)) {
    rlang::abort(
      c(
        "`station_id` must be an atomic vector (character, factor, or numeric).",
        i = paste0("You supplied an object of class: ",
                   paste(class(station_id), collapse = ", "))
      ),
      call = call
    )
  }

  if (length(station_id) != R) {
    rlang::abort(
      c(
        paste0("`station_id` has length ", length(station_id),
               " but `counts` has ", R, " rows."),
        i = "There must be exactly one station label per row of `counts`.",
        i = "With replicates, each row of `counts` is one replicate, and rows from the same station repeat that station's label."
      ),
      call = call
    )
  }

  if (anyNA(station_id)) {
    rlang::abort(
      c(
        paste0("`station_id` contains ", sum(is.na(station_id)), " NA value(s)."),
        i = "Every replicate must be assigned to a station.",
        i = "Drop those rows from `counts` and `station_id` together, or label them."
      ),
      call = call
    )
  }

  station_chr <- as.character(station_id)
  # Levels in order of first appearance, so station ordering is stable and
  # predictable (row 1's station is station 1).
  levels_ <- unique(station_chr)
  index   <- match(station_chr, levels_)
  n_stations <- length(levels_)

  if (n_stations < 2) {
    rlang::abort(
      c(
        paste0("`station_id` identifies only ", n_stations, " station(s)."),
        i = "The mixture model needs at least 2 stations (ideally many more).",
        i = "Check that `station_id` labels stations, not replicates: all replicates of one station share a label."
      ),
      call = call
    )
  }

  reps <- table(factor(station_chr, levels = levels_))
  max_reps <- max(reps)

  # No replication anywhere: alpha and phi are not separately identified, and
  # the replicate model would buy the user nothing. Signal a fallback rather
  # than warning or erroring — unreplicated data stays fully supported.
  if (max_reps == 1) {
    return(list(
      index      = index,
      levels     = levels_,
      n_stations = n_stations,
      reps       = reps,
      fallback   = TRUE
    ))
  }

  n_replicated <- sum(reps > 1)
  if (n_replicated < 3) {
    rlang::inform(
      c(
        paste0("Only ", n_replicated, " station(s) have more than one replicate."),
        i = "The replicate-level dispersion `phi` is estimated from replicated stations only, so it will lean heavily on its prior.",
        i = "Singleton stations are still fitted normally (partial pooling); this is a note, not a problem."
      )
    )
  }

  list(
    index      = index,
    levels     = levels_,
    n_stations = n_stations,
    reps       = reps,
    fallback   = FALSE
  )
}

#' Collapse replicate-level covariates to one row per station
#'
#' Accepts covariates supplied either per station (already `n_stations` rows) or
#' per replicate (one row per row of `counts`). In the latter case the rows are
#' collapsed to one per station, erroring if any covariate varies within a
#' station — the model has no replicate-level covariate term, so such a covariate
#' could not be used.
#'
#' @keywords internal
collapse_station_covariates <- function(covariates, station, call = rlang::caller_env()) {
  if (is.null(covariates)) return(NULL)

  n_rows <- if (is.data.frame(covariates)) nrow(covariates) else NROW(covariates)

  # Already station-level: nothing to do.
  if (n_rows == station$n_stations) return(covariates)

  if (n_rows != length(station$index)) {
    rlang::abort(
      c(
        paste0("`covariates` has ", n_rows, " rows, which matches neither the number of stations (",
               station$n_stations, ") nor the number of replicate rows (",
               length(station$index), ")."),
        i = "Supply covariates either one row per station (in order of first appearance in `station_id`), or one row per row of `counts`.",
        i = "Covariates are station-level: the model has no replicate-level covariate term."
      ),
      call = call
    )
  }

  # First row of each station, in station-index order.
  first_row <- match(seq_len(station$n_stations), station$index)
  collapsed <- if (is.data.frame(covariates)) {
    covariates[first_row, , drop = FALSE]
  } else {
    covariates[first_row, , drop = FALSE]
  }

  # Constancy check: a covariate that varies within a station cannot be a
  # station-level covariate, and silently taking the first value would be wrong.
  varying <- character(0)
  col_names <- colnames(covariates)
  if (is.null(col_names)) col_names <- paste0("column ", seq_len(NCOL(covariates)))

  for (j in seq_len(NCOL(covariates))) {
    col <- if (is.data.frame(covariates)) covariates[[j]] else covariates[, j]
    ok <- vapply(split(col, station$index), function(v) {
      if (length(v) < 2) return(TRUE)
      if (is.numeric(v)) {
        isTRUE(all.equal(v, rep(v[1], length(v)), tolerance = 1e-8))
      } else {
        all(v == v[1])
      }
    }, logical(1))
    if (!all(ok)) varying <- c(varying, col_names[j])
  }

  if (length(varying) > 0) {
    rlang::abort(
      c(
        paste0("Covariate(s) vary between replicates of the same station: ",
               paste(varying, collapse = ", ")),
        i = "Covariates in this model act on the station's community membership, so they must be constant within a station.",
        i = "Either aggregate them yourself (e.g. take the station mean) or drop them from `covariates`."
      ),
      call = call
    )
  }

  if (is.data.frame(collapsed)) rownames(collapsed) <- NULL
  collapsed
}

#' @keywords internal
validate_K <- function(K, N, call = rlang::caller_env()) {
  if (!is.numeric(K) || length(K) != 1 || K != round(K)) {
    rlang::abort(
      c(
        "`K` must be a single positive integer (the number of communities to fit).",
        i = paste0("You supplied: K = ", paste(K, collapse = ", ")),
        i = "Example: `K = 3` fits a model with 3 latent communities."
      ),
      call = call
    )
  }
  K <- as.integer(K)
  if (K < 2) {
    rlang::abort(
      c(
        "`K` must be at least 2.",
        i = "K = 1 is a single-community model (no mixture), which is not meaningful here.",
        i = "Start with K = 2 and increase if needed."
      ),
      call = call
    )
  }
  if (K >= N) {
    rlang::abort(
      c(
        paste0("`K` (", K, ") must be less than the number of samples N (", N, ")."),
        i = "You cannot have more communities than samples.",
        i = paste0("With ", N, " samples, maximum K is ", N - 1, ".",
                   " Typical values are K = 2 to 6.")
      ),
      call = call
    )
  }
  if (K > 10) {
    rlang::warn(
      c(
        paste0("K = ", K, " is unusually large."),
        i = "Models with many communities are prone to label switching and slow convergence.",
        i = "Consider fitting K = 2 through K = 6 and using LOO-CV to select."
      )
    )
  }
  K
}

#' @keywords internal
validate_covariate_names <- function(covariates, covariate_names, call = rlang::caller_env()) {
  if (is.null(covariate_names)) return(invisible(NULL))

  current_names <- colnames(covariates)
  if (is.null(current_names)) {
    rlang::warn(
      c(
        "`covariates` has no column names. The `covariate_names` argument will be used as labels.",
        i = "To avoid this, set colnames on your covariates matrix: `colnames(covariates) <- c('depth', 'lat')`"
      )
    )
    return(invisible(NULL))
  }

  bad <- setdiff(covariate_names, current_names)
  if (length(bad) > 0) {
    # Suggest near-matches (typos)
    suggestions <- vapply(bad, function(b) {
      dists   <- adist(b, current_names, ignore.case = TRUE)
      nearest <- current_names[which.min(dists)]
      if (min(dists) <= 2) paste0("Did you mean '", nearest, "'?") else ""
    }, character(1))
    sugg_msgs <- suggestions[nchar(suggestions) > 0]

    rlang::abort(
      c(
        paste0("Covariate name(s) not found in `covariates`: ",
               paste(bad, collapse = ", ")),
        i = paste0("Available names: ", paste(current_names, collapse = ", ")),
        if (length(sugg_msgs) > 0) i = paste(sugg_msgs, collapse = "; ") else NULL
      ),
      call = call
    )
  }
}

#' @keywords internal
check_stan_fit_object <- function(x, fn_name = "this function", call = rlang::caller_env()) {
  if (!inherits(x, "edna_dmm_fit")) {
    rlang::abort(
      c(
        paste0("`fit` must be an `edna_dmm_fit` object (the output of `eDNA_dmm()`)."),
        i = paste0("You supplied an object of class: ", paste(class(x), collapse = ", ")),
        i = paste0("Run `fit <- eDNA_dmm(counts, covariates, K = 2)` first, then pass `fit` to `", fn_name, "`.")
      ),
      call = call
    )
  }
}

#' @keywords internal
make_community_colors <- function(k) {
  hues <- seq(15, 375, length.out = k + 1)[seq_len(k)]
  setNames(grDevices::hcl(h = hues, c = 80, l = 55),
           paste0("Community ", seq_len(k)))
}

#' Structured taxon palette: hue families, shaded within each
#'
#' Taxa are sorted and then laid out hue block by hue block, so neighbouring
#' names share a hue and differ in lightness. That is what lets a legend of
#' thirty-odd taxa be read as a handful of colour families rather than as
#' thirty arbitrary swatches.
#'
#' Counts are spread across ALL `n_hues` hues before any shading is added.
#' The earlier version built a shades-by-hues grid and read it column-major,
#' which meant a taxon count that the first few hues could absorb left the
#' last hue unused: 30 taxa became six hues of five shades and the seventh
#' hue, the pink one, never appeared in any figure.
#'
#' @param taxa Character vector of taxon names. `"Other"` is ignored here and
#'   appended as grey by `include_other`.
#' @param n_hues Number of hue families. Default `7`.
#' @param include_other Append `Other = "grey70"`. Default `TRUE`.
#' @return A named character vector of colours.
#' @keywords internal
make_taxa_colors <- function(taxa, n_hues = 7, include_other = TRUE) {
  other <- if (isTRUE(include_other)) c("Other" = "grey70") else NULL
  taxa_sorted <- sort(taxa[taxa != "Other"])
  n <- length(taxa_sorted)
  if (n == 0) return(other)

  n_hues    <- min(n_hues, n)
  base_hues <- seq(15, 375, length.out = n_hues + 1)[seq_len(n_hues)]
  # Even split of n taxa over n_hues hues, e.g. 30 over 7 gives 5,4,5,4,4,4,4.
  per_hue   <- diff(round(seq(0, n, length.out = n_hues + 1)))

  cols <- unlist(lapply(seq_len(n_hues), function(j) {
    k <- per_hue[[j]]
    if (k == 0) return(character(0))
    lum <- if (k == 1) 58 else seq(75, 40, length.out = k)
    grDevices::hcl(h = base_hues[[j]], c = 80, l = lum)
  }))

  c(stats::setNames(cols[seq_len(n)], taxa_sorted), other)
}
