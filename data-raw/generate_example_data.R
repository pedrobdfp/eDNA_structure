# =============================================================================
# Generate the eDNA_structure example dataset
# Run this script ONCE to create data/example_edna.rda
# This file is not part of the package itself; it's a data-generation script.
# =============================================================================
#
# The dataset simulates a coastal eDNA metabarcoding survey:
#   - 60 sampling stations along a latitudinal gradient (39N–41.5N)
#   - 3 depth strata: 50 m, 150 m, 300 m
#   - 3 sampling years: 2019, 2021, 2023
#   - 25 taxa (mix of fish, invertebrates — plausible NE Pacific fauna)
#   - TRUE underlying structure: K=3 communities driven by depth and latitude
#
# This is simulated data for illustration only. Taxon names are real species
# but read counts are synthetic.
# =============================================================================

set.seed(2025)

# ── 1. Sample metadata ────────────────────────────────────────────────────────
n_samples   <- 60
latitudes   <- seq(39.1, 41.4, length.out = 20)
depths      <- c(50, 150, 300)
years       <- c(2019, 2021, 2023)

# Create a balanced design: 20 latitude positions x 3 depths
cov_raw <- expand.grid(
  latitude = latitudes,
  depth    = depths,
  stringsAsFactors = FALSE
)[1:n_samples, ]

# Assign years roughly evenly across the spatial gradient
cov_raw$year <- rep(years, length.out = n_samples)
cov_raw$depth_bin <- ifelse(cov_raw$depth <= 100, "shallow", "deep")
cov_raw$sample_id <- paste0("STN_", sprintf("%03d", seq_len(n_samples)))
rownames(cov_raw)  <- cov_raw$sample_id

# ── 2. Taxon names (25 plausible NE Pacific taxa) ─────────────────────────────
taxa <- c(
  "Engraulis mordax",       "Sardinops sagax",        "Merluccius productus",
  "Sebastes entomelas",     "Sebastes melanops",      "Scomber japonicus",
  "Trachurus symmetricus",  "Citharichthys sordidus", "Paralichthys californicus",
  "Anoplopoma fimbria",     "Theragra chalcogramma",  "Gadus macrocephalus",
  "Stenobrachius leucopsarus", "Diaphus theta",        "Lampanyctus ritteri",
  "Clupea pallasii",        "Thunnus alalunga",       "Mola mola",
  "Orcinus orca",           "Tursiops truncatus",     "Zalophus californianus",
  "Octopus bimaculatus",    "Dosidicus gigas",        "Loligo opalescens",
  "Pleuronectes vetulus"
)
S <- length(taxa)  # 25

# ── 3. True community compositions (K=3) ──────────────────────────────────────
# Community 1: shallow, southern — small pelagic fish dominant
pi1 <- c(0.25, 0.20, 0.10, 0.05, 0.05, 0.08, 0.07, 0.02, 0.02, 0.01,
          0.01, 0.01, 0.02, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.005,
          0.005, 0.005, 0.005, 0.005, 0.005)
pi1 <- pi1 / sum(pi1)

# Community 2: deep, northern — groundfish and mesopelagic fish
pi2 <- c(0.02, 0.02, 0.15, 0.12, 0.10, 0.02, 0.03, 0.05, 0.04, 0.08,
          0.10, 0.07, 0.05, 0.04, 0.04, 0.01, 0.005, 0.005, 0.005, 0.005,
          0.005, 0.01, 0.01, 0.01, 0.02)
pi2 <- pi2 / sum(pi2)

# Community 3: mid-depth, transitional — mixed mesopelagic and squid
pi3 <- c(0.05, 0.05, 0.08, 0.05, 0.05, 0.04, 0.05, 0.03, 0.03, 0.04,
          0.04, 0.03, 0.08, 0.07, 0.06, 0.02, 0.02, 0.01, 0.01, 0.01,
          0.005, 0.04, 0.06, 0.07, 0.03)
pi3 <- pi3 / sum(pi3)

pi_true <- rbind(pi1, pi2, pi3)
rownames(pi_true) <- paste0("Community_", 1:3)

# ── 4. Membership probabilities (softmax driven by depth and latitude) ─────────
lat_z   <- (cov_raw$latitude - mean(cov_raw$latitude)) / sd(cov_raw$latitude)
depth_z <- (cov_raw$depth    - mean(cov_raw$depth))    / sd(cov_raw$depth)

# True betas:
# Community 1: positive lat (southern), negative depth (shallow)
# Community 2: negative lat (northern), positive depth (deep)
# Community 3: reference

beta1_intercept <- 0.3;  beta1_lat <- -1.5; beta1_depth <- -2.0
beta2_intercept <- -0.2; beta2_lat <-  1.5; beta2_depth <-  2.0

eta <- cbind(
  beta1_intercept + beta1_lat * lat_z + beta1_depth * depth_z,
  beta2_intercept + beta2_lat * lat_z + beta2_depth * depth_z,
  rep(0, n_samples)
)

softmax_fn <- function(x) {
  x <- x - max(x)
  exp(x) / sum(exp(x))
}

gamma_mat <- t(apply(eta, 1, softmax_fn))

# ── 5. Simulate count data ────────────────────────────────────────────────────
# Dirichlet-Multinomial with alpha = 3 (moderate overdispersion)
alpha_true <- 3.0
n_reads    <- sample(2000:8000, n_samples, replace = TRUE)

# Use GDD (generative draws) per sample
# For each sample: draw community assignment then draw DirMult counts
counts_mat <- matrix(0L, nrow = n_samples, ncol = S,
                     dimnames = list(cov_raw$sample_id, taxa))

rdirichlet <- function(alpha_vec) {
  x <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)
  x / sum(x)
}

for (i in seq_len(n_samples)) {
  z_i     <- sample(1:3, 1, prob = gamma_mat[i, ])
  alpha_k <- alpha_true * pi_true[z_i, ]
  theta_i <- rdirichlet(alpha_k)
  counts_mat[i, ] <- as.integer(round(rmultinom(1, n_reads[i], theta_i)))
}

# ── 6. Assemble output ────────────────────────────────────────────────────────
covariates_df <- cov_raw[, c("sample_id", "latitude", "depth", "year", "depth_bin")]

example_edna <- list(
  counts     = counts_mat,
  covariates = covariates_df,
  # Include true parameters so users can validate their fits
  true_params = list(
    K        = 3,
    pi       = pi_true,
    gamma    = gamma_mat,
    beta_lat   = c(beta1_lat, beta2_lat),
    beta_depth = c(beta1_depth, beta2_depth)
  )
)

# Save to package data directory
usethis::use_data(example_edna, overwrite = TRUE, compress = "xz")
message("Saved: data/example_edna.rda")
