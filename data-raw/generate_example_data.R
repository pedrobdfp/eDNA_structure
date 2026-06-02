# =============================================================================
# Generate the eDNA_structure example dataset using the simulation pipeline
# Run once from the package root: source("data-raw/generate_example_data.R")
# =============================================================================

devtools::load_all()

set.seed(2026)

community_covariate_means <- matrix(
  c(80, 200,
    10, 200,
    80,  20,
    10,  20),
  nrow = 4, byrow = TRUE,
  dimnames = list(NULL, c("Depth", "Distance_shore"))
)

example_edna <- simulate_eDNA_survey(
  n_communities             = 4,
  n_species                 = 40,
  samples_per_community     = 5,
  community_covariate_means = community_covariate_means,
  mean_read_depth           = 10000,
  bio_reps                  = 1,
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

cat("Count matrix:", paste(dim(example_edna$counts), collapse = " x "), "\n")
cat("Read depth range:", paste(range(rowSums(example_edna$counts)), collapse = "-"), "\n")

usethis::use_data(example_edna, overwrite = TRUE, compress = "xz")
message("Saved: data/example_edna.rda")
