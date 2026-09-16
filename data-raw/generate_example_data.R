# =============================================================================
# Generate the eDNAstructure example datasets
#
# Run once from the package root:
#   source("data-raw/generate_example_data.R")
#
# Produces two small, deliberately plain objects:
#
#   example_edna       - an unreplicated survey: one row of counts per site
#   example_edna_reps  - the same survey with 3 replicates per site
#
# Both are trimmed to what a user actually needs in order to fit a model. The
# full simulation return object (contributor tables, per-organism shedding,
# raw metabarcoding long format, composition attributes) is deliberately NOT
# shipped: it is an implementation detail of simulate_eDNA_survey(), not an
# example of what eDNA data look like. Users who want it call
# simulate_eDNA_survey() themselves.
# =============================================================================

devtools::load_all(".")

set.seed(2026)

community_covariate_means <- matrix(
  c(80, 200,
    10, 200,
    80,  20,
    10,  20),
  nrow = 4, byrow = TRUE,
  dimnames = list(NULL, c("Depth", "Distance_shore"))
)

# Shared simulation settings. The seed and parameters match the original
# example dataset, so `counts` is unchanged from previous versions.
sim_args <- list(
  n_communities             = 4,
  n_species                 = 40,
  samples_per_community     = 5,
  community_covariate_means = community_covariate_means,
  mean_read_depth           = 10000,
  seq_reps                  = 1,
  spillover                 = 0.15,
  shedding_error            = 0.3,
  decay_rate                = 0.1,
  seed                      = 2026,
  n_dominant_range          = c(2, 3),
  n_unique_low_freq         = 4,
  community_groups          = c(1, 2, 1, 2),
  n_group_shared            = 4,
  group_shared_presence     = 1.0,
  dominant_freq             = c(0.08, 0.14),
  low_unique_freq           = c(0.02, 0.04),
  group_shared_freq         = c(0.02, 0.04),
  shared_freq               = c(0.02, 0.04),
  shared_presence           = 0
)

# -----------------------------------------------------------------------------
# 1. Unreplicated example: one row per site
# -----------------------------------------------------------------------------

sim <- do.call(simulate_eDNA_survey, c(sim_args, list(bio_reps = 1)))

counts <- sim$counts
storage.mode(counts) <- "integer"

# Covariates: numeric only, so the object can be passed straight to eDNA_dmm()
# without the user first having to drop ID or ground-truth columns.
covariates <- data.frame(
  Depth          = as.numeric(sim$covariates$Depth),
  Distance_shore = as.numeric(sim$covariates$Distance_shore),
  row.names      = rownames(counts)
)

# Metadata: the labelling columns, kept separate from the model covariates.
# This is what the plotting functions take as their `metadata` argument.
metadata <- data.frame(
  sample_id        = rownames(counts),
  TrueCommunity    = as.integer(sim$metadata$TrueCommunity),
  Depth            = as.numeric(sim$covariates$Depth),
  Distance_shore   = as.numeric(sim$covariates$Distance_shore),
  stringsAsFactors = FALSE
)

example_edna <- list(
  counts     = counts,
  covariates = covariates,
  metadata   = metadata
)

# -----------------------------------------------------------------------------
# 2. Replicated example: three replicates per site
# -----------------------------------------------------------------------------

sim_rep <- do.call(simulate_eDNA_survey, c(sim_args, list(bio_reps = 3)))

rep_long         <- sim_rep$metab_df
rep_long$station <- paste0("STN_", sprintf("%03d", rep_long$SampleID))
rep_long$rep_key <- paste0(rep_long$station, "_B", rep_long$BioRep)
rep_long$Species <- paste0("Sp_", rep_long$Species)

rep_wide <- tidyr::pivot_wider(
  rep_long[, c("rep_key", "station", "Species", "Counts")],
  names_from  = "Species",
  values_from = "Counts",
  values_fill = 0L
)
rep_wide <- rep_wide[order(rep_wide$station, rep_wide$rep_key), ]

rep_counts <- as.matrix(rep_wide[, !(names(rep_wide) %in% c("rep_key", "station"))])
rownames(rep_counts)     <- rep_wide$rep_key
storage.mode(rep_counts) <- "integer"

rep_replication <- rep_wide$station

# Station-level covariates: one row per station, in order of first appearance
station_levels <- unique(rep_replication)
idx <- match(station_levels,
             sim_rep$metadata$sample_id)

rep_covariates <- data.frame(
  Depth          = as.numeric(sim_rep$covariates$Depth)[idx],
  Distance_shore = as.numeric(sim_rep$covariates$Distance_shore)[idx],
  row.names      = station_levels
)

rep_metadata <- data.frame(
  sample_id        = station_levels,
  TrueCommunity    = as.integer(sim_rep$metadata$TrueCommunity)[idx],
  Depth            = as.numeric(sim_rep$covariates$Depth)[idx],
  Distance_shore   = as.numeric(sim_rep$covariates$Distance_shore)[idx],
  stringsAsFactors = FALSE
)

example_edna_reps <- list(
  counts     = rep_counts,
  replication = rep_replication,
  covariates = rep_covariates,
  metadata   = rep_metadata
)

# -----------------------------------------------------------------------------
# Report and save
# -----------------------------------------------------------------------------

cat("example_edna\n")
cat("  counts    :", paste(dim(example_edna$counts), collapse = " x "), "\n")
cat("  read depth:", paste(range(rowSums(example_edna$counts)), collapse = "-"), "\n")

cat("example_edna_reps\n")
cat("  counts    :", paste(dim(example_edna_reps$counts), collapse = " x "), "\n")
cat("  stations  :", length(unique(example_edna_reps$replication)), "\n")
cat("  reps/stn  :", paste(range(table(example_edna_reps$replication)), collapse = "-"), "\n")

save(example_edna,      file = "data/example_edna.rda",      compress = "xz")
save(example_edna_reps, file = "data/example_edna_reps.rda", compress = "xz")
message("Saved: data/example_edna.rda, data/example_edna_reps.rda")
