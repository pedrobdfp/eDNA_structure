# =============================================================================
# dmm_beta_draws(): softmax coefficients against a chosen reference community
# =============================================================================


#' Posterior draws of the softmax coefficients, against any reference community
#'
#' @description
#' Returns the full `K x (P+1)` coefficient matrix for every posterior draw,
#' by default in the **centred** parameterization [eDNA_dmm()] actually fits
#' (each community's coefficient is its deviation from the average community,
#' nothing is fixed at zero), optionally re-expressed as a contrast against a
#' chosen reference community instead.
#'
#' The softmax is invariant to adding a constant to all communities within a
#' covariate, so centred and reference-coded views are the same model in
#' different coordinates. Passing `reference` re-expresses every coefficient
#' as community `k` relative to that one community, which then becomes exactly
#' zero (contrast coding). This is sometimes an easier reading than the
#' default deviation-from-average view, e.g. contrasting every community
#' against the most distinct one.
#'
#' Changing the reference requires **no refitting**: subtracting community
#' `j`'s row from every row moves the zero and leaves the model untouched.
#' Doing it on the draws rather than on the summaries means the intervals
#' carry their uncertainty across correctly.
#'
#' @param fit An `edna_dmm_fit` object from [eDNA_dmm()].
#' @param reference Integer community index to re-express the coefficients
#'   against (that community's row becomes exactly zero), or `NULL` (default)
#'   to keep the fitted centred parameterization, in which every community
#'   has a real, generally nonzero coefficient.
#' @param n_draws Maximum number of posterior draws to return, or `NULL`
#'   (default) for all of them.
#'
#' @return A numeric array with dimensions `draws x K x (P+1)`. The second
#'   dimension is named `"Community 1" ... "Community K"`, the third
#'   `c("intercept", fit$covariate_names)`. With `reference` set, that
#'   community's slice is exactly zero; otherwise every slice generally has a
#'   nonzero value (the centred parameterization has no fixed-zero row).
#'
#' @seealso [eDNA_dmm_beta()] and [eDNA_dmm_beta_intervals()], which both take
#'   a `reference` argument and use this internally.
#'
#' @examples
#' \dontrun{
#' fit <- eDNA_dmm(counts, covariates, K = 4)
#'
#' # Coefficients against community 2 instead of the fitted reference
#' d <- dmm_beta_draws(fit, reference = 2)
#' dim(d)                       # draws x 4 x (P+1)
#' colMeans(d[, , "depth"])     # posterior mean depth effect per community
#'
#' # The plotting functions take the same argument
#' eDNA_dmm_beta_intervals(fit, reference = 2)$plot
#' }
#'
#' @export
dmm_beta_draws <- function(fit, reference = NULL, n_draws = NULL) {
  check_stan_fit_object(fit, "dmm_beta_draws")
  K <- fit$K
  P <- length(fit$covariate_names)
  if (P == 0)
    rlang::abort(c(
      "This model has no covariates (intercept-only).",
      i = "There are no covariate coefficients to re-reference."
    ))
  cov_labels <- c("intercept", fit$covariate_names)

  if (!is.null(fit$beta_draws_aligned)) {
    # A fit that has been through dmm_align_labels() carries draws whose
    # community labels mean the same thing in every draw. The raw draws in
    # `stan_fit` do NOT: with several chains each one numbered its communities
    # arbitrarily, so rebuilding from `stan_fit` would pool draws that are
    # about different communities. That does not error, it just silently
    # inflates every interval and pulls the estimates toward zero. Prefer the
    # aligned draws whenever they exist.
    out <- fit$beta_draws_aligned
    if (!identical(dimnames(out)[[3]], cov_labels))
      dimnames(out) <- list(NULL, paste0("Community ", seq_len(K)), cov_labels)
    if (!is.null(n_draws) && n_draws < dim(out)[1])
      out <- out[seq_len(n_draws), , , drop = FALSE]
  } else {
    beta_mat <- as.matrix(fit$stan_fit, pars = "beta")
    n <- nrow(beta_mat)
    if (!is.null(n_draws)) n <- min(n, n_draws)

    # How many community rows does this model carry? A reference-coded model
    # stores K - 1 and the reference row is restored below as explicit zeros;
    # a sum-to-zero model stores all K and needs no restoration. Reading it off
    # the parameter names rather than assuming keeps both codings working
    # through the same path.
    n_rows <- suppressWarnings(max(as.integer(
      sub("^beta\\[([0-9]+),.*$", "\\1", colnames(beta_mat)))))
    if (!is.finite(n_rows)) n_rows <- K - 1L

    # Rebuild the full K x (P+1) matrix per draw, with the fitted reference's
    # row restored as explicit zeros so every community is represented.
    out <- array(0, dim = c(n, K, P + 1),
                 dimnames = list(NULL, paste0("Community ", seq_len(K)), cov_labels))
    for (comm_i in seq_len(min(n_rows, K))) {
      for (cov_j in seq_len(P + 1)) {
        pname <- sprintf("beta[%d,%d]", comm_i, cov_j)
        if (!pname %in% colnames(beta_mat)) next
        out[, comm_i, cov_j] <- beta_mat[seq_len(n), pname]
      }
    }
  }

  if (!is.null(reference)) {
    reference <- as.integer(reference)
    if (length(reference) != 1 || is.na(reference) ||
        reference < 1 || reference > K)
      rlang::abort(sprintf("`reference` must be a single integer in 1..%d.", K))
    ref <- out[, reference, , drop = FALSE]
    for (k in seq_len(K)) out[, k, ] <- out[, k, ] - ref[, 1, ]
  }
  out
}
