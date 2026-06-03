# eDNAstructure

**Dirichlet-Multinomial Mixture Models for eDNA Metabarcoding Community Structure**

`eDNA_structure` is an R package for fitting Bayesian mixture models to
environmental DNA (eDNA) read count data. The model identifies latent
ecological communities in your metabarcoding data, estimates their taxonomic
compositions, and quantifies how environmental variables (depth, latitude,
temperature, etc.) predict community membership.

---

## What it does

Given a **sample × taxon read count matrix** and optional environmental
covariates, `eDNAstructure`:

1. Fits a **Dirichlet-Multinomial Mixture (DMM) model** using Bayesian
   inference (Stan under the hood)
2. Returns posterior community membership probabilities for every sample
3. Estimates the taxonomic composition of each community
4. Quantifies covariate effects through softmax regression with full
   uncertainty propagation
5. Produces publication-ready `ggplot2` visualizations

---

## Installation

### Prerequisites

`eDNAstructure` requires R ≥ 4.1.0 and a working Stan/rstan installation.
If you have not installed rstan before, follow these steps first:

**1. Install a C++ toolchain:**
- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS**: Run `xcode-select --install` in Terminal
- **Linux**: Install `build-essential` (Ubuntu/Debian) or equivalent

**2. Install rstan:**
```r
install.packages("rstan")
```

**3. Verify Stan works:**
```r
library(rstan)
example(stan_model, package = "rstan", run.dontrun = TRUE)
```
If you see sampling output without errors, Stan is ready.

Full rstan installation guide: <https://mc-stan.org/rstan/articles/rstan.html>

### Install eDNA_structure

Several chunks of bleack text will appear. They are related to the stan code pre-compilling.

```r
# Install the remotes package if needed
install.packages("remotes")

# Install eDNA_structure from GitHub
remotes::install_github("pedrobdfp/eDNA_structure")
```

### Dependencies

Installed automatically with the package:

| Package | Purpose |
|---------|---------|
| `rstan` (≥ 2.21) | Bayesian inference via Stan |
| `ggplot2` (≥ 3.4) | All visualizations |
| `dplyr`, `tidyr` | Data manipulation |
| `vegan` (≥ 2.6) | NMDS ordination and community ecology distances |
| `posterior` (≥ 1.4) | MCMC diagnostic statistics (ESS, Rhat) |
| `loo` (≥ 2.6) | Leave-one-out cross-validation for K selection |
| `scales` | Axis label formatting |

---

## Quick start

```r
library(eDNAstructure)

# Load the built-in example dataset
data <- get_example_data()
str(data)

# Fit a K=3 DMM with depth and latitude as covariates
fit <- eDNA_dmm(
  counts     = data$counts,
  covariates = data$covariates[, c("latitude", "depth")],
  K          = 3
)

# Print a summary
print(fit)
summary(fit)

# Structure plot
eDNA_dmm_structure(
  fit,
  metadata  = data$covariates,
  facet_var = "depth_bin",
  sort_var  = "latitude"
)

# NMDS ordination
eDNA_dmm_nmds(fit)$plot

# Covariate effects (prior vs posterior)
eDNA_dmm_beta(fit)$plot
```

---

## Input data format

### Count matrix

Your primary input should be a **sample × taxon count matrix**:

- **Rows** = samples (stations, replicates, individuals, etc.)
- **Columns** = taxa or ASVs — taxonomic annotation is not required
- **Values** = non-negative integer read counts

```r
# Example: what the count matrix looks like
head(data$counts[, 1:4])
#          Engraulis mordax  Sardinops sagax  Merluccius productus  Sebastes entomelas
# STN_001              412              310                   121                  73
# STN_002              389              275                    98                  61
# STN_003               52               41                   487                 312
```

If your data are in long format (sample, taxon, count columns), convert them:

```r
library(tidyr)
count_matrix <- long_df |>
  pivot_wider(names_from = taxon, values_from = reads, values_fill = 0) |>
  tibble::column_to_rownames("sample_id") |>
  as.matrix()
```

### Covariate data frame

```r
head(data$covariates)
#   sample_id  latitude  depth  year  depth_bin
# 1   STN_001     39.10     50  2019    shallow
# 2   STN_002     39.22     50  2021    shallow
# 3   STN_003     39.34    150  2023       deep
```

Covariates are Z-score standardized internally by default.
The covariate data frame must be in the same row order as the count matrix.

---

## Core functions

### `eDNA_dmm()` — Fit the model

```r
fit <- eDNA_dmm(
  counts           = my_counts,    # sample × taxon count matrix
  covariates       = my_covs,      # sample × covariate data frame (or NULL)
  K                = 3,            # number of communities
  scale_covariates = TRUE,         # Z-score covariates (strongly recommended)
  chains           = 1,            # single chain (see ?eDNA_dmm)
  iter             = 4000,         # total MCMC iterations
  warmup           = 2000,         # warmup iterations to discard
  adapt_delta      = 0.95,         # HMC step-size target (increase if divergences)
  max_treedepth    = 12,           # increase if "max treedepth" warnings
  seed             = 13,           # random seed for reproducibility
  conc             = 0.5,          # Dirichlet concentration (< 1 = sparse)
  alpha_shape      = 5,            # Gamma prior shape for overdispersion
  alpha_rate       = 2,            # Gamma prior rate (mean = shape/rate = 2.5)
  verbose          = TRUE          # print Stan sampling progress
)
```

Returns an `edna_dmm_fit` object with elements including:
- `sample_info`: posterior community membership per sample
- `pi_mean`: posterior mean community compositions (K × S matrix)
- `beta_summary`: covariate coefficient summaries with ESS and reliability
- `stan_fit`: the raw `rstan::stanfit` object for advanced diagnostics

### `eDNA_dmm_structure()` — Structure bar plot

```r
p <- eDNA_dmm_structure(
  fit,
  metadata         = my_metadata,   # data frame with extra sample variables
  sample_id_col    = "sample_id",   # column in metadata matching fit sample IDs
  facet_var        = "year",        # create panels by this variable
  sort_var         = "latitude",    # sort samples within panels
  community_colors = NULL,          # auto-generated; or named hex vector
  bar_width        = 0.9,           # width of each sample bar (0–1)
  x_text           = FALSE,         # show sample labels on x-axis?
  base_size        = 11,            # base font size
  title            = NULL,          # custom title (NULL = auto)
  subtitle         = NULL,          # custom subtitle (NULL = auto)
  legend_position  = "bottom"       # "bottom", "right", "left", "top", "none"
)
```

Returns a `ggplot2` object.

### `eDNA_dmm_nmds()` — NMDS ordination

```r
result <- eDNA_dmm_nmds(
  fit,
  k                = 2,            # NMDS dimensions (2 or 3)
  nmds_axes        = c(1, 2),      # which axes to plot
  distance         = "bray",       # dissimilarity index (see ?vegan::vegdist)
  use_edna_index   = TRUE,         # apply eDNA index transformation before ordination
  trymax           = 100,          # max random NMDS starts
  seed             = 42,           # random seed for NMDS
  show_ellipse     = TRUE,         # draw 95% confidence ellipses?
  ellipse_type     = "t",          # "t" or "norm"
  community_colors = NULL,         # named hex vector or NULL for auto
  size_range       = c(1.5, 5),    # point size range (50% to 100% certainty)
  alpha            = 0.85,         # point transparency
  base_size        = 13,           # base font size
  title            = NULL,         # custom title
  subtitle         = NULL,         # custom subtitle
  legend_position  = "right"       # legend placement
)

result$plot   # ggplot2 object
result$nmds   # vegan metaMDS object (for species scores etc.)
```

### `eDNA_dmm_beta()` — Covariate coefficient plots

```r
result <- eDNA_dmm_beta(
  fit,
  layout             = "joint",     # "joint" (communities overlaid) or "separate" (faceted rows)
  covariates_to_plot = NULL,        # character vector of covariate names, or NULL for all
  show_intercept     = FALSE,       # include the intercept coefficient?
  beta_prior_sd      = 1.0,         # SD of the Normal(0, sd) prior — must match Stan model
  n_prior_samples    = 4000,        # draws from prior for smooth density curve
  community_colors   = NULL,        # named hex vector or NULL
  prior_color        = "grey60",    # fill color for prior density
  prior_alpha        = 0.35,        # transparency of prior density
  posterior_alpha    = 0.55,        # transparency of posterior density
  show_annotations   = TRUE,        # annotate panels with mean, CI, P(direction)
  base_size          = 13,          # base font size
  title              = NULL,        # custom title
  subtitle           = NULL         # custom subtitle
)

result$plot    # ggplot2 object
result$table   # data frame with coefficient summaries
```

### `eDNA_loo()` — LOO-CV for K selection

```r
loo_result <- eDNA_loo(
  counts      = my_counts,
  covariates  = my_covs,
  K_range     = 2:5,          # fit models for K = 2, 3, 4, 5
  iter        = 4000,
  warmup      = 2000,
  seed        = 13
)

loo_result$plot        # elbow plot (LOO-ELPD vs K)
loo_result$loo_table   # data frame with ELPD and SE per K
loo_result$fits        # list of edna_dmm_fit objects, one per K
```

### `get_example_data()` — Built-in example dataset

```r
data <- get_example_data()
# data$counts      — 60 × 25 count matrix
# data$covariates  — data frame: sample_id, latitude, depth, year, depth_bin
# data$true_params — true generating parameters (K=3, pi, gamma)
```

---

## The model

The generative model for sample *i*:

1. **Community compositions**: π_k ~ Dirichlet(conc · **1**_S) for each community *k*
2. **Membership probability**: P(*z*_i = *k*) = softmax(β_0k + β_1k · cov1_i + … + β_Pk · covP_i), with community *K* as reference
3. **Observed counts**: **x**_i | *z*_i = *k* ~ DirichletMultinomial(N_i, α · π_k)

*z*_i is marginalized out, making inference exact. The global overdispersion
parameter α absorbs both technical (PCR, sequencing) and ecological variance.

---

## Frequently Asked Questions

**Why use one chain instead of many?**

Mixture models have a symmetry: "Community 1" in one chain may correspond to
"Community 2" in another (label switching). Multi-chain Rhat diagnostics
appear pathological even when each chain converges perfectly. A single long
chain avoids this. Use within-chain ESS as your convergence metric instead
(reported in `summary()`).

**My model has divergent transitions. What do I do?**

Try `adapt_delta = 0.99`. If divergences persist, try fewer communities (lower
K) or check that your count matrix has no all-zero rows or columns.

**Can I use raw ASVs instead of taxonomically-collapsed data?**

Yes. The model treats each column as a "type" — it does not use taxonomy at
all. Fitting on ASVs gives finer resolution; fitting on taxa reduces
dimensionality and often converges faster.

**My taxa have unfamiliar names / no taxonomy. Is that a problem?**

No. Column names are used only as labels in plots.

**How do I include year as a covariate?**

Year can be included as a numeric column in your covariate data frame. However,
if you only have a few discrete time points (e.g., 3 years), treating year as
continuous assumes a linear trend. Consider fitting the model without year and
then testing year effects post-hoc (e.g., multinomial regression of community
assignment on year).

---

## Citation

If you use this package in published research, please cite:

> Brandão-Dias et al. (year). Multinomial mixture models from environmental DNA reveal 
> depth stability and dynamic surface turnover of marine vertebrate communities. Under review
 
> Brandão-Dias et al. (year). eDNAstructure: Dirichlet-Multinomial Mixture Models for
> eDNA Metabarcoding Community Structure. R package version 0.1.0.
> https://github.com/pedrobdfp/eDNA_structure

Please also cite Stan:

> Carpenter B. et al. (2017). Stan: A probabilistic programming language.
> *Journal of Statistical Software*, 76(1).

---

## License

MIT © eDNAstructure authors
