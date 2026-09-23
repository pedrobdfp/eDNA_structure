# =============================================================================
# Generate the eDNAstructure example dataset
#
# Run once from the package root:
#   source("data-raw/generate_example_data.R")
#
# Produces example_edna: a small, deliberately plain, unreplicated survey (one
# row of counts per site). Trimmed to what a user actually needs in order to
# fit a model. The full simulation return object (contributor tables,
# per-organism shedding, raw metabarcoding long format, composition
# attributes) is deliberately NOT shipped: it is an implementation detail of
# simulate_eDNA_survey(), not an example of what eDNA data look like. Users
# who want it call simulate_eDNA_survey() themselves.
# =============================================================================

devtools::load_all(".")

set.seed(42)

# Community means, chosen so that:
#  - every community sits several SDs from 0 (no physically impossible
#    negative Depth after sampling noise)
#  - the between-community gap is a moderate multiple of the SD (~5), enough
#    for a clean visual separation in plot_true_compositions() without being
#    so extreme that the covariates become (quasi-)perfectly separating in
#    the softmax regression, which starves beta of gradient and produces poor
#    chain mixing (see DIFFERENCES_FROM_PUBLIC-style debugging notes: the
#    original 80/10 x 200/20 design was 14-36 SDs apart at sd = 5, which is
#    why beta would not converge at K = 4 no matter how long chains ran).
community_covariate_means <- matrix(
  c(70, 140,
    40, 140,
    70, 100,
    40, 100),
  nrow = 4, byrow = TRUE,
  dimnames = list(NULL, c("Depth", "Distance_shore"))
)

covariate_sds <- matrix(
  c(6, 8),
  nrow = 4, ncol = 2, byrow = TRUE,
  dimnames = list(NULL, c("Depth", "Distance_shore"))
)

# Shared simulation settings.
sim_args <- list(
  n_communities             = 4,
  n_species                 = 40,
  samples_per_community     = 5,
  community_covariate_means = community_covariate_means,
  covariate_sds             = covariate_sds,
  mean_read_depth           = 10000,
  seq_reps                  = 1,
  spillover                 = 0.15,
  shedding_error            = 0.3,
  decay_rate                = 0.1,
  seed                      = 42,
  n_dominant_range          = c(2, 3),
  n_unique_low_freq         = 4,
  community_groups          = c(1, 2, 1, 2),
  n_group_shared            = 4,
  group_shared_presence     = 1.0,
  dominant_freq             = c(0.08, 0.14),
  low_unique_freq           = c(0.02, 0.04),
  group_shared_freq         = c(0.02, 0.04),
  shared_freq               = c(0.02, 0.04),
  shared_presence           = 0,
  bio_reps                  = 1
)

sim <- do.call(simulate_eDNA_survey, sim_args)

counts <- sim$counts
storage.mode(counts) <- "integer"

# Covariates: numeric only, so the object can be passed straight to eDNA_dmm()
# without the user first having to drop ID or ground-truth columns. Rounded
# to 3 significant figures: these are illustrative depths/distances, not
# measurements, and the extra decimals only add visual noise in the
# quickstart/tutorial output.
covariates <- data.frame(
  Depth          = signif(as.numeric(sim$covariates$Depth), 3),
  Distance_shore = signif(as.numeric(sim$covariates$Distance_shore), 3),
  row.names      = rownames(counts)
)

# Metadata: the labelling columns, kept separate from the model covariates.
# This is what the plotting functions take as their `metadata` argument.
metadata <- data.frame(
  sample_id        = rownames(counts),
  TrueCommunity    = as.integer(sim$metadata$TrueCommunity),
  Depth            = covariates$Depth,
  Distance_shore   = covariates$Distance_shore,
  stringsAsFactors = FALSE
)

example_edna <- list(
  counts     = counts,
  covariates = covariates,
  metadata   = metadata
)

# -----------------------------------------------------------------------------
# Report and save
# -----------------------------------------------------------------------------

cat("example_edna\n")
cat("  counts    :", paste(dim(example_edna$counts), collapse = " x "), "\n")
cat("  read depth:", paste(range(rowSums(example_edna$counts)), collapse = "-"), "\n")
cat("  Depth     :", paste(range(covariates$Depth), collapse = "-"), "\n")
cat("  Shore     :", paste(range(covariates$Distance_shore), collapse = "-"), "\n")

save(example_edna, file = "data/example_edna.rda", compress = "xz")
message("Saved: data/example_edna.rda")

# -----------------------------------------------------------------------------
# Flat CSV mirror, for the "Getting Started" vignette's "one file" walkthrough
# -----------------------------------------------------------------------------

raw_csv <- cbind(
  data.frame(
    sample_id     = rownames(counts),
    Depth         = covariates$Depth,
    Distance_shore = covariates$Distance_shore,
    TrueCommunity = metadata$TrueCommunity
  ),
  as.data.frame(counts)
)
write.csv(raw_csv, "inst/extdata/example_edna_raw.csv", row.names = FALSE, quote = FALSE)
message("Saved: inst/extdata/example_edna_raw.csv")
