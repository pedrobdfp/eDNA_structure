# eDNAstructure

**Dirichlet-Multinomial Mixture Models for eDNA Metabarcoding Community Structure**

`eDNAstructure` is an R package for fitting Bayesian Dirichlet-Multinomial Mixture (DMM) models to environmental DNA (eDNA) read count data from metabarcoding surveys. Given a sample × taxon count matrix and optional environmental covariates, the model identifies latent ecological communities, estimates their taxonomic compositions, and quantifies how environmental gradients (depth, latitude, temperature, etc.) drive community membership — all within a fully Bayesian framework with principled uncertainty quantification.

---

## Installation

### Step 1 — Install a C++ toolchain

Stan compiles models to C++ and requires a toolchain on your machine:

- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS**: Run `xcode-select --install` in Terminal
- **Linux**: Install `build-essential` (Ubuntu/Debian) or equivalent

### Step 2 — Install rstan

```r
install.packages("rstan")
```

Verify it works before proceeding:

```r
library(rstan)
example(stan_model, package = "rstan", run.dontrun = TRUE)
```

If you see sampling output without errors, Stan is ready. Full guide: <https://mc-stan.org/rstan/articles/rstan.html>

### Step 3 — Install eDNAstructure

```r
install.packages("remotes")
remotes::install_github("pedrobdfp/eDNA_structure")
```

### Dependencies

Installed automatically:

| Package | Purpose |
|---------|---------|
| `rstan` (≥ 2.21) | Bayesian inference via Stan |
| `ggplot2` (≥ 3.4) | All visualizations |
| `dplyr`, `tidyr` | Data manipulation |
| `vegan` (≥ 2.6) | NMDS ordination |
| `posterior` (≥ 1.4) | MCMC diagnostics (ESS, Rhat) |
| `loo` (≥ 2.6) | Leave-one-out cross-validation |
| `scales` | Axis formatting |

---

## Quick start

```r
library(eDNAstructure)

# Load the built-in example dataset (or simulate your own — see below)
data <- get_example_data()

# Fit a K=4 DMM with depth and distance from shore as covariates
fit <- eDNA_dmm(
  counts     = data$counts,
  covariates = data$covariates[, c("Depth", "Distance_shore")],
  K          = 4
)

print(fit)
summary(fit)

# Structure plot — one bar per sample, colored by community membership probability
eDNA_dmm_structure(
  fit,
  metadata  = data$covariates,
  facet_var = "TrueCommunity",
  sort_var  = "Depth"
)

# NMDS ordination colored by community assignment
eDNA_dmm_nmds(fit)$plot

# Prior vs posterior distributions for covariate effects
eDNA_dmm_beta(fit)$plot
```

For a complete walkthrough — simulation, data formatting, K selection, all plot options, and parameter recovery — see the **[full tutorial](vignettes/tutorial.Rmd)**.

---

## Input data format

### Count matrix

The primary input to `eDNA_dmm()` is a **sample × taxon** matrix of non-negative integer read counts:

- **Rows** = samples (stations, replicates, individuals, etc.)
- **Columns** = taxa or ASVs — taxonomic annotation is not required
- **Values** = raw integer read counts (do not normalize)

```r
data$counts[1:3, 1:5]
#          Sp_1  Sp_2  Sp_3  Sp_4  Sp_5
# STN_001   412   310   121    73     0
# STN_002   389   275    98    61    14
# STN_003    52    41   487   312   208
```

If your data are in long format, convert them first:

```r
library(tidyr)
count_matrix <- long_df |>
  pivot_wider(names_from = taxon, values_from = reads, values_fill = 0) |>
  tibble::column_to_rownames("sample_id") |>
  as.matrix()
```

### Covariate data frame

A **sample × covariate** data frame in the same row order as the count matrix. All covariates are Z-score standardized internally by default.

```r
head(data$covariates)
#   sample_id  TrueCommunity  Depth  Distance_shore
# 1   STN_001              1     82             198
# 2   STN_002              1     79             204
# 3   STN_003              2     11             197
```

---

## Functions

### `eDNA_dmm()` — Fit the DMM

The core function. Fits a Dirichlet-Multinomial Mixture model via Stan and returns an `edna_dmm_fit` object.

```r
fit <- eDNA_dmm(
  counts           = my_counts,   # sample × taxon integer count matrix
  covariates       = my_covs,     # sample × covariate data frame, or NULL
  K                = 4,           # number of latent communities to fit
  scale_covariates = TRUE,        # Z-score standardize covariates (strongly recommended)
  chains           = 1,           # number of MCMC chains (see note on label switching below)
  iter             = 4000,        # total iterations per chain (including warmup)
  warmup           = 2000,        # warmup iterations to discard
  adapt_delta      = 0.95,        # HMC target acceptance rate; increase to 0.99 if divergences
  max_treedepth    = 12,          # increase to 14–15 if "max treedepth exceeded" warnings
  seed             = 13,          # random seed for reproducibility
  conc             = 0.5,         # Dirichlet prior concentration: < 1 = sparse communities
  alpha_shape      = 5,           # Gamma prior shape for overdispersion parameter alpha
  alpha_rate       = 2,           # Gamma prior rate  (prior mean = shape/rate = 2.5)
  verbose          = TRUE         # print sampling progress
)
```

The returned `edna_dmm_fit` object contains:

| Element | Description |
|---------|-------------|
| `sample_info` | Data frame: posterior membership probabilities and MAP assignment per sample |
| `pi_mean` | Matrix [K × S]: posterior mean community compositions |
| `beta_summary` | Data frame: covariate coefficient summaries with ESS and reliability |
| `alpha_mean` | Scalar: posterior mean overdispersion |
| `stan_fit` | Raw `rstan::stanfit` object for advanced diagnostics |

> **On single chains:** Mixture models suffer from label switching across chains — "Community 1" in chain A may correspond to "Community 2" in chain B, making multi-chain Rhat diagnostics meaningless. A single long chain sidesteps this entirely. Use within-chain ESS (reported by `summary()`) as your convergence criterion.

---

### `eDNA_dmm_structure()` — Structure bar plot

Produces a STRUCTURE-style plot: one vertical bar per sample, divided into colored segments by posterior community membership probability. The same format used in population genetics structure plots.

```r
p <- eDNA_dmm_structure(
  fit,
  metadata         = my_metadata,    # data frame with additional sample variables
  sample_id_col    = "sample_id",    # column in metadata matching sample IDs in fit
  facet_var        = "TrueCommunity", # facet panels by this variable (e.g. site, year, depth)
  sort_var         = "Depth",        # sort samples within each panel by this variable
  community_colors = NULL,           # named hex vector (e.g. c("Community 1" = "#E63946"))
                                     # or NULL for automatic HCL palette
  bar_width        = 0.9,            # bar width (0–1); 1 = no gaps between bars
  x_text           = FALSE,          # show sample ID labels on x-axis?
  base_size        = 11,             # base font size in points
  title            = NULL,           # plot title; NULL = auto-generated
  subtitle         = NULL,           # plot subtitle; NULL = auto-generated
  legend_position  = "bottom"        # "bottom", "right", "left", "top", or "none"
)
```

Returns a `ggplot2` object — save with `ggsave()` or extend with additional ggplot2 layers.

---

### `eDNA_dmm_nmds()` — NMDS ordination

Runs NMDS on Bray-Curtis dissimilarities (with the eDNA index transformation by default) and plots samples colored by their MAP community assignment. Point size scales with assignment certainty.

```r
result <- eDNA_dmm_nmds(
  fit,
  k                = 2,          # NMDS dimensions (2 or 3); increase if stress > 0.2
  nmds_axes        = c(1, 2),    # which two axes to display; e.g. c(1,3) for axes 1 and 3
  distance         = "bray",     # dissimilarity metric passed to vegan::vegdist()
  use_edna_index   = TRUE,       # apply eDNA index transform before computing distances
  trymax           = 100,        # maximum random NMDS starts (more = less risk of local optima)
  seed             = 42,         # random seed for NMDS
  show_ellipse     = TRUE,       # draw 95% confidence ellipse per community?
  ellipse_type     = "t",        # ellipse type: "t" (robust) or "norm" (normal-based)
  community_colors = NULL,       # named hex vector or NULL for automatic palette
  size_range       = c(1.5, 5),  # point size range mapping 50% to 100% certainty
  alpha            = 0.85,       # point transparency
  base_size        = 13,
  title            = NULL,
  subtitle         = NULL,
  legend_position  = "right"
)

result$plot   # ggplot2 object
result$nmds   # vegan::metaMDS object (access species scores, stress, etc.)
```

---

### `eDNA_dmm_beta()` — Covariate effects

Overlays the prior and posterior distributions for each softmax regression coefficient. A posterior pulled away from the prior is evidence that the covariate genuinely predicts community membership.

```r
result <- eDNA_dmm_beta(
  fit,
  layout             = "joint",    # "joint": communities overlaid per covariate panel
                                   # "separate": one row per community, one column per covariate
  covariates_to_plot = NULL,       # character vector of covariate names to include, or NULL for all
  show_intercept     = FALSE,      # include the intercept term?
  beta_prior_sd      = 1.0,        # prior SD — must match the Stan model (default: Normal(0,1))
  n_prior_samples    = 4000,       # prior draws for the density curve (more = smoother)
  community_colors   = NULL,       # named hex vector or NULL
  prior_color        = "grey60",   # fill color for the prior density
  prior_alpha        = 0.35,       # prior density transparency
  posterior_alpha    = 0.55,       # posterior density transparency
  show_annotations   = NULL,       # NULL = auto (shown for K=2 only); TRUE or FALSE to override
  base_size          = 13,
  title              = NULL,
  subtitle           = NULL
)

result$plot    # ggplot2 object
result$table   # data frame: mean, 90% CI, P(direction), ESS, reliability per coefficient
```

---

### `eDNA_loo()` — K selection via LOO cross-validation

Fits models across a range of K values and compares them using Leave-One-Out cross-validation. Returns an elbow plot and a comparison table to guide K selection.

```r
loo_result <- eDNA_loo(
  counts           = my_counts,
  covariates       = my_covs,
  K_range          = 2:5,        # integer vector of K values to evaluate
  scale_covariates = TRUE,
  chains           = 1,
  iter             = 4000,
  warmup           = 2000,
  adapt_delta      = 0.95,
  seed             = 13,
  conc             = 0.5,
  alpha_shape      = 5,
  alpha_rate       = 2,
  verbose          = TRUE
)

loo_result$plot        # LOO-ELPD elbow plot (higher = better; look for the elbow)
loo_result$loo_table   # data frame: K, LOO-ELPD, SE
loo_result$loo_compare # loo::loo_compare() output
loo_result$fits        # named list of edna_dmm_fit objects, one per K
```

---

### `plot_true_compositions()` — Raw species composition

Visualizes the observed species frequencies per sample as stacked bars — the same layout as `eDNA_dmm_structure()`, allowing direct before/after comparison. Most useful before fitting to inspect the raw community signal, and with simulated data where true community labels are known.

```r
p <- plot_true_compositions(
  counts,
  metadata        = my_metadata,   # data frame for faceting and sorting
  sample_id_col   = "sample_id",
  facet_var       = "TrueCommunity", # facet by known or hypothesized grouping
  sort_var        = "Depth",
  top_n           = 20,             # show top N taxa individually; rest collapsed to "Other"
  bar_width       = 0.95,
  base_size       = 11,
  title           = NULL,
  subtitle        = NULL,
  legend_position = "none"          # default none — too many taxa for a useful legend
)
```

Returns a `ggplot2` object.

---

### `get_example_data()` — Built-in example dataset

Returns the built-in simulated dataset: 20 samples × 40 taxa across 4 communities separated by depth and distance from shore. Generated by `simulate_eDNA_survey()` with known ground truth, so fitted parameters can be compared to the true values.

```r
data <- get_example_data()
# data$counts                — integer matrix [20 × 40]
# data$covariates            — data frame: sample_id, TrueCommunity, Depth, Distance_shore
# data$community_compositions — true composition matrix [4 × 40]
# data$metab_df              — raw simulated metabarcoding reads
# data$sample_metadata       — full simulation metadata
```

---

### Simulation pipeline

`eDNAstructure` includes a mechanistic simulation pipeline for generating eDNA datasets with known community structure. Use it for method validation, power analysis, or teaching.

```r
# Full pipeline in one call
sim <- simulate_eDNA_survey(
  n_communities         = 4,
  n_species             = 40,
  samples_per_community = 5,
  seed                  = 2026
)
# sim$counts      — ready for eDNA_dmm()
# sim$covariates  — ready for eDNA_dmm()

# Or run each step individually for full control:
community_mat   <- generate_community_compositions(n_communities = 4, n_species = 40, seed = 1)
contrib_obj     <- generate_contributors(community_mat, samples_per_community = 5)
eDNA_obj        <- generate_eDNA(contrib_obj$contributors_list)
metab_df        <- simulate_metabarcoding(eDNA_obj, mean_read_depth = 10000)
sample_metadata <- generate_sample_covariates(contrib_obj$contributors_list,
                                               community_covariates = my_cov_means)
```

| Function | Purpose |
|----------|---------|
| `simulate_eDNA_survey()` | Full pipeline in one call |
| `generate_community_compositions()` | K community frequency vectors over S species |
| `generate_contributors()` | Organisms shedding eDNA per sample |
| `generate_eDNA()` | Shedding, exponential decay, bottle sub-sampling |
| `simulate_metabarcoding()` | Amplification bias and multinomial read counts |
| `generate_sample_covariates()` | Environmental metadata drawn from community-specific distributions |

---

## The model

For sample *i*, the DMM marginalizes over a latent community assignment *z*_i:

1. **Compositions**: π_k ~ Dirichlet(conc · **1**_S) for k = 1…K
2. **Membership**: P(*z*_i = k) = softmax(β_0k + β_1k · x_1i + … + β_Pk · x_Pi), community K = reference
3. **Counts**: **x**_i | *z*_i = k ~ DirichletMultinomial(N_i, α · π_k)

The global overdispersion α absorbs both technical (PCR, sequencing) and ecological compositional variance. Marginalizing over *z*_i makes inference exact and avoids discrete sampling.

---

## Frequently asked questions

**Why only one chain?**
Label switching: "Community 1" in chain A may map to "Community 2" in chain B. Multi-chain Rhat values are pathological even when each chain converges perfectly. One long chain avoids this. Check within-chain ESS instead (printed by `summary()`).

**I have divergent transitions. What do I do?**
Increase `adapt_delta` toward `0.99`. If they persist, try lower K or verify your count matrix has no all-zero rows or columns.

**Can I use raw ASVs instead of taxonomy-collapsed counts?**
Yes. The model treats each column as a compositional unit and does not use taxonomy. ASVs give finer resolution; taxa collapse dimensionality and often converge faster.

**How do I include year as a covariate?**
Pass it as a numeric column. But if you have only a few discrete years, the linearity assumption may be too strong — consider fitting without year and testing it post-hoc via multinomial regression on the posterior assignments.

**The first run takes forever — is something wrong?**
No. Stan compiles the model to C++ on the first call (1–2 minutes). All subsequent calls in the same session skip compilation entirely. Adding `rstan_options(auto_write = TRUE)` to your script (or `.Rprofile`) caches the compiled model across sessions.

---

## Citation

If you use this package in published research, please cite:

> Brandão-Dias et al. (year). Multinomial mixture models from environmental DNA reveal
> depth stability and dynamic surface turnover of marine vertebrate communities. *Under review.*

> Brandão-Dias et al. (year). eDNAstructure: Dirichlet-Multinomial Mixture Models for
> eDNA Metabarcoding Community Structure. R package version 0.1.0.
> https://github.com/pedrobdfp/eDNA_structure

Please also cite Stan:

> Carpenter B. et al. (2017). Stan: A probabilistic programming language.
> *Journal of Statistical Software*, 76(1).

---

## License

MIT © eDNAstructure authors
