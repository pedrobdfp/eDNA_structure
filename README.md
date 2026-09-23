# eDNAstructure

**Dirichlet-Multinomial Mixture Models for eDNA Metabarcoding Community Structure**

`eDNAstructure` is an R package for fitting Bayesian Dirichlet-Multinomial Mixture (DMM) models to environmental DNA (eDNA) read count data from metabarcoding surveys. Given a sample x taxon count matrix and optional environmental covariates, the model identifies latent ecological communities, estimates their taxonomic compositions, and quantifies how environmental gradients correlate with community membership, all within a fully Bayesian framework with principled uncertainty quantification.

---

## Installation

### Step 1: Install a C++ toolchain

Stan compiles models to C++ and requires a toolchain on your machine:

- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS**: Run `xcode-select --install` in Terminal
- **Linux**: Install `build-essential` (Ubuntu/Debian) or equivalent

### Step 2: Install rstan

```r
install.packages("rstan")
```

Verify it works before proceeding:

```r
library(rstan)
example(stan_model, package = "rstan", run.dontrun = TRUE)
```

If you see sampling output without errors, Stan is ready. Stan can be a little annoying to install, but they have an excellent guide:
Full guide: <https://mc-stan.org/rstan/articles/rstan.html>

### Step 3: Install eDNAstructure

```r
install.packages("remotes")
remotes::install_github("pedrobdfp/eDNA_structure", upgrade = "never")
```

During installation, a large amount of black text may appear; this is the Stan model compiling to C++. It only happens once. Every subsequent call to `eDNA_dmm()` goes straight to sampling with no compilation output.

### Dependencies

Installed automatically:

| Package | Purpose |
|---------|---------|
| `rstan` (>= 2.21) | Bayesian inference via Stan |
| `ggplot2` (>= 3.4) | All visualizations |
| `dplyr`, `tidyr` | Data manipulation |
| `vegan` (>= 2.6) | NMDS ordination |
| `posterior` (>= 1.4) | MCMC diagnostics (ESS, Rhat) |
| `loo` (>= 2.6) | Leave-one-out cross-validation |
| `scales` | Axis formatting |

---

## Quick start
A quick look at the core package functions
```r
library(eDNAstructure)

# The built-in example: a site-by-taxon count table plus two numeric covariates
data <- get_example_data()

data$counts[1:3, 1:5]   # 20 sites x 32 taxa: this is all the model needs
head(data$covariates)   # Depth, Distance_shore (numeric only)

# Look at the raw (observed) composition before fitting anything
plot_true_compositions(
  data$counts,
  metadata  = data$metadata,      # labelling columns live here
  facet_var = "TrueCommunity",
  show_legend        = TRUE
)

# Fit the model for a given number of communities (K)
# This is the main function of the package, which fits the DMM model to the data.
fit <- eDNA_dmm(
  counts     = data$counts,
  covariates = data$covariates,   # goes in as-is, nothing to subset
  K          = 4
)

# Some useful summaries
print(fit)
summary(fit)

# Plot the results
# First the structure plot. Its default is to plot the exact same way as plot_true_compositions()
# Note that the community labels are arbitrary
eDNA_dmm_structure(fit, metadata = data$metadata,
                   facet_var = "TrueCommunity", sort_var = "Depth")			   
eDNA_dmm_compositions(fit)   # what each community looks like in species space
eDNA_dmm_nmds(fit)$plot      # ordination, colored by assignment
eDNA_dmm_beta_intervals(fit)$plot  # covariate effects, with credible intervals
```

### Choosing K

In the example dataset, the true number of communities we have simulated is k=4. However, this will be unknown for any real dataset. So, the first step in inference is figuring out how many communities there are.
```r
# This function will call eDNA_dmm() and fit the data with several different Ks,
# and then compare the fits with LOO.
# This is computationally demanding and may take several minutes
loo_result <- eDNA_loo(data$counts, data$covariates, K_range = 2:5)
# You will see some warnings associated with k=3 and k=5. 
# This may happen when K is blatantly incorrect, which may cause the model to not converge.
loo_result$plot        # ELPD against K, unfilled where chains disagreed

# eDNA_dmm_k_diagnostics() is the fuller picture behind that one elbow plot:
# six panels (predictive fit, marginal gain, reliability, convergence,
# assignment certainty, community distinctness) that fail in different ways,
# so K is rarely ambiguous once you see them together. Feed it loo_result
# straight, it pulls out everything it needs on its own.
diagnostics <- eDNA_dmm_k_diagnostics(loo_result)
diagnostics$plot

# Predictive fit and identifiability are separate questions:
loo_result$loo_table[, c("K", "elpd", "lp_rhat", "min_agreement", "identified")]

# All fitted models are stored: no need to refit
# So when you decide K=4 is the correct one, you don't have to run eDNA_dmm() again for that. 
fit <- loo_result$fits[["K4"]]

# From here, you can do all your downstream visualizations and analyses! 
```

> **For a complete walkthrough**, including step-by-step simulation, data formatting, K selection, all visualization options, parameter recovery, and troubleshooting, see the **[full tutorial vignette](vignettes/tutorial.Rmd)**. It is designed to be read start to finish and assumes no prior familiarity with Bayesian mixture models.

---

## Input data format

The model needs **two things**, and only the first is required:

1. **A species/ASV x site table**: the count matrix. This is the only required input.
2. **A covariate table**: *optional*. Supply it only if you want to model how environmental variables drive community membership. `eDNA_dmm()` runs perfectly well without it.

### Count matrix

The primary input to `eDNA_dmm()` is a **sample x taxon** matrix of non-negative integer read counts:

- **Rows** = samples, one row per sample
- **Columns** = taxa or ASVs, taxonomic annotation is not required
- **Values** = raw integer read counts (do not normalize)

```r
data$counts[1:3, 1:5]
#         Sp_1 Sp_2 Sp_3 Sp_4 Sp_7
# STN_001 2054 1794  780  717 1087
# STN_002  620 3256 1492    0    0
# STN_003  192 3306    0    0    0
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

A **sample x covariate** data frame in the same row order as the count matrix. Covariates are Z-score standardized internally by default.

```r
head(data$metadata)
#   sample_id TrueCommunity Depth Distance_shore
# 1   STN_001             1  83.6            137
# 2   STN_002             1  63.7            137
# 3   STN_003             1  66.8            148
```

---

## Functions

### `eDNA_dmm()`: Fit the DMM

The core function. Fits a Dirichlet-Multinomial Mixture model via Stan and returns an `edna_dmm_fit` object.

```r
fit <- eDNA_dmm(
  counts           = my_counts,   # sample x taxon integer count matrix
  covariates       = my_covs,     # sample x covariate data frame, or NULL
  K                = 4,           # number of latent communities to fit
  scale_covariates = TRUE,        # standardise covariates (strongly recommended)
  chains           = 4,           # MCMC chains; their labels are aligned automatically
  cores            = NULL,        # NULL = one core per chain, capped at physical cores
  method           = "STEPHENS",  # label alignment algorithm
  iter             = 4000,        # total iterations per chain, including warmup
  warmup           = 2000,        # warmup iterations to discard
  adapt_delta      = 0.95,        # HMC target acceptance; raise to 0.99 if divergences
  max_treedepth    = 12,          # raise to 14 or 15 if treedepth warnings appear
  seed             = 13,          # random seed
  conc             = 0.5,         # Dirichlet prior concentration; below 1 = sparse
  alpha_shape      = 5,           # Gamma prior shape for overdispersion alpha
  alpha_rate       = 2,           # Gamma prior rate; prior mean = shape / rate = 2.5
  verbose          = TRUE         # print sampling progress
)
```

The returned `edna_dmm_fit` object contains:

| Element | Description |
|---------|-------------|
| `sample_info` | Data frame: membership probability per community, most probable community, and assignment certainty per sample |
| `pi_mean` | Matrix [K x S]: posterior mean community compositions |
| `beta_summary` | Data frame: coefficient summaries with ESS and reliability |
| `beta_draws_aligned` | Array [draws x K x P+1]: aligned, centred coefficient draws |
| `alignment` | Label alignment method, permutations, and the diagnostics below |
| `alpha_mean` | Scalar: posterior mean overdispersion |
| `stan_fit` | Raw `rstan::stanfit` object for advanced diagnostics |

#### Community labels are aligned for you

A mixture likelihood does not change if the community labels are permuted, so
chains that agree completely about the structure of the data can still
disagree about which community is called 1, 2 or 3. Estimates and convergence
diagnostics computed in that state describe the labelling rather than the fit.

`eDNA_dmm()` aligns every posterior draw onto a common labelling before
computing anything, using the algorithm of Stephens (2000) from the
**label.switching** package. Compositions, membership probabilities and
coefficients are all rebuilt from the aligned draws, so `summary()` reports
convergence for the model rather than for its naming.

Alignment fixes disagreement about names. It cannot fix disagreement about the
grouping itself, so two further numbers are reported that relabelling cannot
affect:

```r
fit$alignment$lp_rhat         # convergence of the log posterior density,
                              # invariant to labelling by construction
fit$alignment$min_agreement   # share of samples the two least similar chains
                              # place in the same community
```

Both close to one means the chains found the same solution. A clear departure
in either means they found different ones, which usually indicates a `K` the
data do not support.

#### Covariate coefficients have no reference community

Adding the same constant to every community's linear predictor leaves the
softmax unchanged, so one degree of freedom per covariate has to be pinned
down. This package pins it by making the coefficients sum to zero across
communities, so each is that community's deviation from the average community
and all `K` are estimated.

To read them against one particular community instead, pass `reference`. No
refitting is needed, because the choice is a change of coordinates applied to
the finished draws:

```r
eDNA_dmm_beta_intervals(fit)                  # deviation from the average
eDNA_dmm_beta_intervals(fit, reference = 3)   # contrasts against community 3
```

---

### `dmm_align_labels()`: Re-align labels yourself

`eDNA_dmm()` calls this automatically; you should not normally need it. It
exists for the cases where you do: re-aligning with a different `method` to
sanity-check the default, or aligning a fit that was assembled by hand from a
raw `stanfit` outside `eDNA_dmm()`.

```r
aligned <- dmm_align_labels(
  fit,
  method  = "STEPHENS",   # "STEPHENS" (default), "ECR-pivot", or "ECR-iterative"
  verbose = TRUE
)

aligned$alignment$pct_permuted   # % of draws that were relabelled
aligned$alignment$min_agreement  # smallest pairwise chain agreement, post-alignment
```

`"STEPHENS"` minimises KL divergence against the mean membership matrix using
full membership probabilities, and is the most reliable of the three.
`"ECR-pivot"` anchors to the highest-density draw's allocation. `"ECR-iterative"`
is faster but can settle into a local optimum where draws are aligned within
groups but mismatched between them, symptoms that look exactly like residual
non-convergence; prefer `"STEPHENS"` unless you have a specific reason not to.

---


### `eDNA_dmm_structure()`: Structure bar plot

Produces a STRUCTURE-style plot: one vertical bar per sample, divided into colored segments by posterior community membership probability.

```r
p <- eDNA_dmm_structure(
  fit,
  metadata         = my_metadata,     # data frame with additional sample variables
  sample_id_col    = "sample_id",     # column in metadata matching sample IDs in fit
  facet_var        = "year",          # column facets
  facet_row_var    = "depth_bin",     # row facets: see below
  sort_var         = "Depth",         # sort samples within each panel
  community_colors = NULL,            # named hex vector, e.g. c("Community 1" = "#E63946")
                                      # or NULL for automatic HCL palette
  bar_width        = 0.9,             # bar width (0-1); 1 = no gaps
  x_text           = FALSE,           # show sample ID labels on x-axis?
  base_size        = 11,
  ylab             = "Membership probability",
  title            = NULL,
  subtitle         = NULL,
  show_legend      = TRUE,
  legend_position  = "bottom"
)
```

Returns a `ggplot2` object, or a **cowplot grid object** when `facet_row_var` is supplied. Save with `ggsave()` or extend with additional ggplot2 layers.

#### Layout, labels and legend

`eDNA_dmm_structure()` shares its entire layout parameter surface with `plot_true_compositions()`, by design, so the two plots can be stacked into one figure and line up exactly. All of these behave identically in both functions:

| Group | Parameters |
|---|---|
| Two-way faceting | `facet_var`, `facet_row_var` |
| Row labels | `row_label_fn`, `row_label_position`, `row_label_size`, `panel_labels` |
| Legend layout | `legend_nrow`, `legend_ncol`, `legend_key_size`, `legend_text_size`, `legend_rel_height`, `legend_rel_width` |
| Axis / strip sizing | `ylab`, `ylab_size`, `ylab_rel_width`, `strip_text_size`, `axis_text_size` |
| Reference lines | `vline_var`, `vline_value`, `vline_color`, `vline_linetype`, `vline_linewidth` |

See [`plot_true_compositions()`](#plot_true_compositions-raw-species-composition) above for the full description of each group. Two differences:

- `ylab` defaults to `"Membership probability"` here (vs `"Proportion"`).
- `axis_text_size` (y-axis tick label size) exists here only. Default `NULL` = ggplot2's own scaling from `base_size`.

#### Stacking the observed and fitted plots

The common case, raw composition on top, community assignment below, sharing panel structure:

```r
data <- get_example_data()
fit  <- eDNA_dmm(counts = data$counts, covariates = data$covariates, K = 4)

# A second grouping variable, to get rows of panels as well as columns
meta <- data$metadata
meta$depth_bin <- ifelse(meta$Depth > 50, 80, 10)

common <- list(
  metadata      = meta,
  facet_var     = "TrueCommunity",   # columns
  facet_row_var = "depth_bin",       # rows
  sort_var      = "Depth",
  row_label_fn  = function(x) paste0(x, " m"),
  base_size     = 11
)

p_obs <- do.call(plot_true_compositions,
                 c(list(data$counts, taxa_include = dmm_taxon_order(fit, 20)), common))
p_fit <- do.call(eDNA_dmm_structure,
                 c(list(fit, panel_labels = FALSE), common))

cowplot::plot_grid(p_obs, p_fit, ncol = 1, rel_heights = c(1, 1))
```

Set `panel_labels = FALSE` on the inner plots when the outer figure supplies its own (a), (b) letters, otherwise you get two competing sets.

---


### `eDNA_dmm_nmds()`: NMDS ordination

Runs NMDS on community dissimilarities and plots samples colored by their MAP community assignment. Point **size** reflects assignment certainty: larger points are more confidently assigned to a single community.

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
  size_range       = c(1.5, 5),  # point size range: c(min, max) mapped to
                                 # 50% certainty (smallest) -> 100% certainty (largest)
  alpha            = 0.85,       # point transparency (0 = invisible, 1 = opaque)
  base_size        = 13,
  title            = NULL,
  subtitle         = NULL,
  legend_position  = "right"
)

result$plot   # ggplot2 object
result$nmds   # vegan::metaMDS object (access stress value, species scores, etc.)
```

---

### `eDNA_dmm_beta()`: Covariate effects

Overlays the prior and posterior distributions for each softmax regression coefficient. A posterior pulled away from the prior is evidence that the covariate genuinely predicts community membership.

```r
result <- eDNA_dmm_beta(
  fit,
  layout             = "joint",    # "joint": communities overlaid per covariate panel
                                   # "separate": one row per community, one column per covariate
  covariates_to_plot = NULL,       # character vector of covariate names to include, or NULL for all
  show_intercept     = FALSE,      # include the intercept term?
  reference          = "none",     # "none" (default): centred, deviation from the average community;
                                   # an integer instead contrasts against that community
  beta_prior_sd      = 1.0,        # prior SD: must match the Stan model (default: Normal(0,1))
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

### `eDNA_dmm_beta_intervals()`: Covariate effects as a coefficient plot

A companion to `eDNA_dmm_beta()`. Densities show how far the data moved a
coefficient from its prior, but the tail of a density is hard to read against
zero. This draws the same coefficients as a coefficient plot instead: one row
per community x covariate, a point estimate, two nested credible intervals,
and a dashed line at zero, so "does this interval cross zero" is a glance, not
a squint. This is what the Quick Start example above calls.

```r
result <- eDNA_dmm_beta_intervals(
  fit,
  covariates_to_plot = NULL,       # character vector of covariate names, or NULL for all
  show_intercept     = FALSE,
  intervals          = c(0.5, 0.9), # inner (thick) and outer (thin) interval widths
  reference          = "none",      # "none" (default): centred, deviation from the average community;
                                    # an integer instead contrasts against that community
  point_est          = "median",    # "median" (default) or "mean"
  community_colors   = NULL,
  color_by_community = TRUE,
  facet_scales       = "free_x",    # or "fixed" to compare magnitudes across covariates
  point_size         = 2.8,
  linewidth_inner    = 1.6,
  linewidth_outer    = 0.6,
  base_size          = 13,
  title              = NULL,
  subtitle           = NULL
)

result$plot    # ggplot2 object
result$table   # data frame: community, covariate, estimate, interval bounds, excludes_zero
```

A **filled** point means the outer interval excludes zero (a real effect,
given the chosen credibility level); a **hollow** point means it crosses zero.

---

### `dmm_beta_draws()`: Raw posterior draws of the beta coefficients

What `eDNA_dmm_beta()` and `eDNA_dmm_beta_intervals()` call internally to get
the `reference` argument's re-centering right. Use it directly when you want
the coefficient draws themselves, e.g. to compute a custom summary, rather
than a plot.

```r
d <- dmm_beta_draws(
  fit,
  reference = NULL,   # NULL (default): the fitted centred draws, every community
                      # has a real, generally nonzero coefficient; an integer
                      # instead re-expresses everything as a contrast against
                      # that community, which becomes exactly zero
  n_draws   = NULL    # cap the number of draws returned, or NULL for all
)

dim(d)                        # draws x K x (P+1); third dim named "intercept", covariate names
colMeans(d[, , "Depth"])      # posterior mean Depth effect per community
```

Changing `reference` needs **no refitting**: a softmax is invariant to adding
a constant to every community within a covariate, so this is a normalisation
of the finished draws, not a re-estimate, applied on the draws (not the
summaries) so intervals carry their uncertainty across correctly.

---

### `eDNA_loo()`: K selection via LOO cross-validation

Fits models across a range of K values and compares them using Leave-One-Out cross-validation. Returns an elbow plot and a comparison table to guide K selection.

```r
loo_result <- eDNA_loo(
  counts           = my_counts,
  covariates       = my_covs,
  K_range          = 2:5,        # integer vector of K values to evaluate
  scale_covariates = TRUE,
  chains           = 4,          # per K; their labels are aligned automatically, like eDNA_dmm()
  method           = "STEPHENS", # label alignment algorithm
  iter             = 4000,
  warmup           = 2000,
  adapt_delta      = 0.95,
  seed             = 13,
  conc             = 0.5,
  alpha_shape      = 5,
  alpha_rate       = 2,
  save_dir         = NULL,       # directory to cache each fit as fit_K<k>.rds, resumable
  keep_fits        = TRUE,       # FALSE releases each fit after scoring: flat memory for long sweeps
  verbose          = TRUE
)

loo_result$plot        # LOO-ELPD elbow plot (higher = better; look for the elbow)
loo_result$loo_table   # data frame: one row per K, elpd/se plus convergence diagnostics
loo_result$loo_by_chain # data frame: one row per chain per K; feeds eDNA_dmm_k_diagnostics()
loo_result$loo_compare # loo::loo_compare() output
loo_result$fits        # named list of edna_dmm_fit objects, one per K (e.g. loo_result$fits[["K4"]])
```

**Predictive fit and identifiability are separate questions.** ELPD almost
always keeps improving as K grows, a more flexible mixture fits anything
better, so it rarely stops on its own. `loo_table$identified` (built from
`lp_rhat` and `min_agreement`) says whether the chains at that K actually
agreed on a single grouping of the samples; an ELPD gain at a K that isn't
identified isn't a model you can report. Read both columns together:

```r
loo_result$loo_table[, c("K", "elpd", "lp_rhat", "min_agreement", "identified")]
max(loo_result$loo_table$K[loo_result$loo_table$identified])   # largest identified K
```

---

### `eDNA_dmm_k_diagnostics()`: The full picture for choosing K

`eDNA_loo()`'s elbow plot is one view of "how many communities." This function
is the rest of the view: six panels, side by side, that fail in different ways
on purpose, so K is rarely ambiguous once you see all of them together. No
single panel picks K; the point is that they usually agree.

```r
# The easy way: hand it eDNA_loo()'s return value directly. fits, elpd_by_run
# and convergence are all pulled out of it automatically.
diagnostics <- eDNA_dmm_k_diagnostics(loo_result)

diagnostics$plot     # the assembled 6-panel figure
diagnostics$panels   # panels (a) to (f) individually, for custom layouts
diagnostics$table    # one row per K: elpd_mean/sd, pct_pareto_bad, min_agreement, lp_rhat,
                      # min_distance, mean/median/q10_certainty, mean_excess, pct_confident
```

The full argument list, for anything more custom (e.g. fits cached on disk
rather than held in memory, or a manual `K_values` subset):

```r
diagnostics <- eDNA_dmm_k_diagnostics(
  fits        = loo_result,             # or a plain named list of edna_dmm_fit,
                                        # or a function(k) for on-disk fits
  K_values    = 2:5,                    # required unless `fits` names them (K2, K3, ...)
  level       = "advanced",             # "advanced" (default, 6 panels) or "simple" (2 panels)
  elpd_by_run = NULL,                   # per-chain ELPD; overrides what `fits` would supply
  adjacent    = NULL,                   # paired K-to-K+1 gain + SE for panel (b); auto-computed
                                        # from `fits` when it's a real list, not a function
  pareto_by_run = NULL,                 # panel (c) inputs; auto-built from elpd_by_run + convergence
  convergence = NULL,                   # lp_rhat / min_agreement; overrides what `fits` would supply

  pareto_threshold    = 0.7,  # Pareto k above which a sample's LOO contribution is unreliable
  rhat_threshold       = 1.1, # lp__ Rhat above which a K is flagged unconverged
  agreement_threshold  = 0.9, # reference line for chain agreement
  certainty_threshold  = 0.8, # membership probability counted as "confidently assigned"
  distance             = "aitchison", # or "tv"; distance between community compositions
  base_size            = 12
)
```

**The panels, left to right, top to bottom:**

| Panel | Question | Deteriorates when |
|---|---|---|
| (a) Predictive fit | Does adding a community improve held-out prediction (ELPD)? | Rarely on its own, this is the ambiguous one |
| (b) Marginal gain | Is the *step* from K-1 to K worth it, paired to cancel sample noise? | The gain falls within ~2 SE of zero |
| (c) Reliability | Can (a) and (b) be trusted? Share of samples with Pareto k above threshold | Too many samples are individually influential at that K |
| (d) Convergence | Did the chains at that K find the *same* grouping? | Chain agreement drops, or `lp__` Rhat exceeds `rhat_threshold` |
| (e) Assignment certainty | Are samples confidently assigned, or spread thin across communities? | Membership probabilities flatten toward 1/K |
| (f) Community distinctness | Is the newest community actually different from the others? | The smallest between-community distance collapses toward 0 |

Panels (a) and (b) ask whether the model predicts better; (c) through (f) ask
whether that improvement is trustworthy and interpretable. It is common for
(a) to keep rising after (d)-(f) have already turned, that gap is exactly the
signal that a higher K is buying predictive accuracy at the cost of a
grouping the data do not actually support. See `?eDNA_dmm_k_diagnostics` for
the full reasoning behind each panel.

---

### `plot_true_compositions()`: Raw species composition

Visualizes observed species frequencies per sample as stacked bars, the same layout as `eDNA_dmm_structure()`, allowing direct before/after comparison. Most useful before fitting to inspect the raw community signal, and with simulated data where true community labels are known.

```r
p <- plot_true_compositions(
  counts,
  metadata        = my_metadata,
  sample_id_col   = "sample_id",
  facet_var       = "year",          # column facets
  facet_row_var   = "depth_bin",     # row facets: see "Two-way layouts" below
  sort_var        = "Depth",         # sort samples within each panel
  top_n           = 20,              # top N taxa individually; rest -> "Other"
  taxa_include    = NULL,            # or a taxon vector from dmm_taxon_order()
  bar_width       = 0.9,
  base_size       = 11,
  ylab            = "Proportion",
  title           = NULL,
  subtitle        = NULL,
  show_legend     = FALSE,
  legend_position = "bottom"
)
```

Returns a `ggplot2` object, or a **cowplot grid object** when `facet_row_var` is supplied.

#### Two-way layouts: `facet_row_var`

`facet_var` makes panel *columns*; `facet_row_var` makes panel *rows*. With `facet_row_var` set, one sub-plot is built per row level and they are stacked with cowplot. That architecture exists for a reason: it is the only one that keeps panel widths proportional to sample count (via `space = "free_x"`), so a stratum with 30 samples is drawn three times wider than one with 10 instead of being stretched to match.

This is the recommended layout for publication figures with two-way structure (e.g. depth stratum x year).

#### Row labels: `row_label_fn`, `row_label_position`

| Parameter | Purpose |
|---|---|
| `row_label_fn` | Function applied to each `facet_row_var` level to build its display label. `function(x) paste0(x, " m")` for depth, `function(x) paste0("Year: ", x)` for year. Default `as.character`. |
| `row_label_position` | `"title"` (bold label above each row's panel: **new default**) or `"ylab"` (as the y-axis label: the old, hardcoded behavior). |
| `row_label_size` | Font size of that label. Default `NULL` = `base_size`. |
| `panel_labels` | Label each row panel with a lowercase letter (a, b, c...). Default `TRUE`. Set `FALSE` when this plot is itself a panel inside a larger composite figure, where the outer figure supplies the letters and an inner set would collide. |

`row_label_fn` is what makes the label text generalize across datasets; it is tied to whatever `facet_row_var` means for your data, rather than being hardcoded.

#### Legend layout

Taxon legends get large fast; these control how the legend is laid out and how much room it is given.

| Parameter | Purpose |
|---|---|
| `legend_nrow`, `legend_ncol` | Rows/columns passed to `guide_legend()`. Default `NULL` = ggplot decides. |
| `legend_key_size` | Size (cm) of each color swatch. Default `0.4`. |
| `legend_text_size` | Font size of legend labels. Default `NULL` = inherits `base_size`. |
| `legend_rel_height` | Height of a top/bottom legend as a fraction of the panel stack. Default `0.1`. |
| `legend_rel_width` | Width of a left/right legend as a fraction of the panel stack. Default `0.2`. |
| `n_hues` | Number of hue families in the taxon palette. Default `7`. |

> ⚠️ **`legend_rel_height` is a hard allocation, not a minimum.** The legend gets exactly that
> fraction of the stack, and any rows that do not fit are **clipped without warning**. If legend
> entries go missing, raise `legend_rel_height` rather than assuming the taxa were dropped.

Taxa are sorted and laid out hue block by hue block, so setting `legend_ncol = n_hues` puts roughly one color family per column, which is what makes a 20-taxon legend scannable.

#### Axis and strip sizing

| Parameter | Purpose |
|---|---|
| `ylab` | Shared y-axis label drawn once beside the panel stack. Used only when `row_label_position = "title"`. Default `"Proportion"`. |
| `ylab_size` | Font size of that label. Default `NULL` = `base_size`. |
| `ylab_rel_width` | Width of the label strip as a fraction of the stack. Default `0.03`. |
| `strip_text_size` | Font size of the column-facet strip labels. Default `NULL` = `base_size`. |

#### Reference lines: `vline_*`

Draws a vertical reference line at a threshold in sample-ordering space, e.g. marking where depth crosses 200 m once samples are sorted by depth.

```r
data <- get_example_data()

p <- plot_true_compositions(
  data$counts, metadata = data$metadata, sort_var = "Depth",
  vline_var       = "Depth",     # numeric column to locate the threshold in
  vline_value     = 50,          # the threshold, in that column's units
  vline_color     = "black",
  vline_linetype  = "dashed",
  vline_linewidth = 0.7
)
```

#### Sharing a palette with the fitted plot: `taxa_include`

```r
data <- get_example_data()
fit  <- eDNA_dmm(data$counts, data$covariates, K = 4)

taxa <- dmm_taxon_order(fit, top_n = 20)
plot_true_compositions(data$counts, taxa_include = taxa)
```

`taxa_include` names exactly which taxa to show individually; everything else is pooled into `"Other"`, and `top_n` is ignored. Names absent from `counts` are dropped with a warning. See [`dmm_taxon_order()`](#dmm_taxon_order-shared-taxon-ranking-for-matching-palettes) for why this matters.

---

---
### `eDNA_dmm_compositions()`: Posterior community compositions

Visualizes the posterior mean taxonomic composition of each latent community as stacked bars, the model's estimate of what each community "looks like" in species space. Colors match those used in `eDNA_dmm_structure()` and `plot_true_compositions()` for direct comparison.

```r
p <- eDNA_dmm_compositions(
  fit,
  top_n           = 20,      # show top N taxa; rest collapsed to "Other"
  base_size       = 13,
  title           = NULL,
  subtitle        = NULL,
  legend_position = "right", # "right", "bottom", "left", "top", or "none"
  bar_width       = 0.7      # bar width (0-1)
)
```

Returns a `ggplot2` object. The x-axis labels show community numbers (1, 2, 3...). Pair with `eDNA_dmm_structure()` to connect community identities to sample assignments.

---

### `dmm_taxon_order()`: Shared taxon ranking (for matching palettes)

Returns taxa ranked by total posterior weight across communities (`colSums(pi_mean)`). A taxon that dominates a single community ranks highly even if it is rare overall, the ordering that matters when the question is what *distinguishes* communities.

```r
taxa <- dmm_taxon_order(
  fit,
  top_n = 20        # NULL (default) returns every taxon, ranked
)
```

**Use this whenever two composition figures sit side by side.** The taxon palette is a deterministic function of the sorted taxon set: same set in, same colors out. If `eDNA_dmm_compositions()` picks its own top 20 and `plot_true_compositions()` independently picks its own, the two sets differ and shared taxa get *different colors in each panel*. Driving both from one ranking fixes that:

```r
data <- get_example_data()
fit  <- eDNA_dmm(data$counts, data$covariates, K = 4)

taxa <- dmm_taxon_order(fit, top_n = 20)

p_obs   <- plot_true_compositions(data$counts, taxa_include = taxa)   # observed
p_model <- eDNA_dmm_compositions(fit, top_n = 20)                # fitted

cowplot::plot_grid(p_obs, p_model, ncol = 1)
```

Note this is deliberately **not** the same ranking as mean observed frequency across samples, which rewards being widespread rather than being diagnostic.

---


### `get_example_data()`: Built-in example dataset

A small, deliberately plain example: a site-by-taxon count table plus the two numeric covariates that separate the communities. 20 sites, 32 taxa, 4 true communities.

```r
data <- get_example_data()

data$counts       # integer matrix, 20 sites x 32 taxa (rows = sites)
data$covariates   # data frame, 20 x 2: Depth, Distance_shore: numeric only
data$metadata     # data frame, 20 x 4: sample_id, TrueCommunity, Depth, Distance_shore
```

That's the whole object, three elements. It goes straight into the model with nothing to subset:

```r
fit <- eDNA_dmm(counts = data$counts, covariates = data$covariates, K = 4)
eDNA_dmm_structure(fit, metadata = data$metadata, facet_var = "TrueCommunity")
```

**Why `covariates` and `metadata` are separate.** `eDNA_dmm()` requires every covariate column to be numeric; it will try to fit whatever you hand it. So an identifier like `sample_id`, or a ground-truth label like `TrueCommunity`, cannot live in `covariates`. The plotting functions want exactly those labelling columns. Keeping the two apart means both calls work as written, with no subsetting.

> **Why 32 taxa and not 40?** The simulation draws 40 species, but taxa with zero reads across every sample are dropped from `counts`; they carry no information and the model rejects all-zero columns.

The simulation's internal tables (contributor lists, per-organism shedding, raw long-format reads) are **not** shipped. They are an implementation detail of `simulate_eDNA_survey()`, not an example of what eDNA data look like. Call [`simulate_eDNA_survey()`](#simulation-pipeline) directly if you want them.

---



### `eDNA_clear_stan_cache()`: Reset the compiled model cache

The Stan model is compiled once on first use and cached per machine under
`tools::R_user_dir("eDNAstructure", "cache")`. Clear it to force a fresh compile, 
useful after upgrading `rstan` or `StanHeaders`, which can leave the cached object
mismatched with the new toolchain.

```r
eDNA_clear_stan_cache()
```

The next call to `eDNA_dmm()` recompiles from source and repopulates the cache.

---


### Simulation pipeline

`eDNAstructure` includes a mechanistic simulation pipeline for generating eDNA datasets with known community structure. Use it for method validation, power analysis, or teaching. The full tutorial (`vignettes/tutorial.Rmd`) walks through the pipeline in detail.

```r
# Full pipeline in one call
sim <- simulate_eDNA_survey(
  n_communities             = 4,           # number of communities K
  n_species                 = 40,          # number of species S
  samples_per_community     = 5,           # sampling stations per community
  community_covariate_means = NULL,        # K x P matrix of covariate means per community
                                           # NULL = default 2-covariate depth x shore design
  covariate_sds             = NULL,        # K x P SDs (NULL = 5 for all)
  mean_read_depth           = 10000,       # mean reads per sample
  bio_reps                  = 1,           # biological replicates per station
  seq_reps                  = 1,           # sequencing technical replicates per bio rep
  spillover                 = 0.15,        # fraction of species frequency leaking between communities
  shedding_error            = 0.3,         # lognormal SD on per-organism eDNA shedding
  decay_rate                = 0.1,         # exponential distance decay of eDNA signal
  seed                      = 42
)
# sim$counts: ready for eDNA_dmm()
# sim$covariates: ready for eDNA_dmm()
```

Or run each step individually for full control:

```r
community_mat   <- generate_community_compositions(
  n_communities     = 4,
  n_species         = 40,
  n_dominant_range  = c(2, 3),     # 2-3 high-frequency dominant species per community
  n_unique_low_freq = 4,           # species present only in one community
  community_groups  = c(1,2,1,2),  # group structure for shared species
  n_group_shared    = 4,           # species shared within each group
  spillover         = 0.15,        # cross-community frequency leakage
  seed              = 1
)

contrib_obj <- generate_contributors(
  community_compositions = community_mat,
  samples_per_community  = 5,
  n_contributors_range   = c(20, 100),   # organisms per sample
  n_distribution         = "Negative Binomial",  # overdispersed organism counts
  size_range             = c(1, 10),     # body size range (affects eDNA shedding)
  distance_range         = c(0, 100),    # distance from sampler (affects eDNA decay)
  seed                   = 1
)

eDNA_obj <- generate_eDNA(
  contributors_list = contrib_obj$contributors_list,
  shedding_rate     = 1000,   # baseline molecules shed per unit size
  beta              = 0.75,   # allometric exponent (metabolic scaling)
  decay_rate        = 0.1,    # exponential distance decay
  shedding_error    = 0.3,    # lognormal noise on shedding (SD on log scale)
  bottle_volume     = 0.1,    # fraction of local eDNA pool captured per bottle
  bio_reps          = 1
)

metab_df <- simulate_metabarcoding(
  eDNA_obj,
  mean_read_depth = 10000,   # mean total reads per sample
  read_depth_sd   = 0.2,     # lognormal SD for read depth variation
  error_sd        = 0.05,    # Gaussian noise added to species frequencies
  rep             = 1        # sequencing technical replicates per sample
)

sample_metadata <- generate_sample_covariates(
  contributors_list    = contrib_obj$contributors_list,
  community_covariates = matrix(           # community-specific covariate means
    c(70, 140, 40, 140, 70, 100, 40, 100),
    nrow = 4, byrow = TRUE,
    dimnames = list(NULL, c("Depth", "Distance_shore"))
  ),
  covariate_sds = matrix(6:8, nrow = 4, ncol = 2, byrow = TRUE)  # Depth SD=6, Distance_shore SD=8
  # Keep the gap between community means a moderate multiple (~5x) of the SD:
  # enough to separate communities cleanly without making covariates
  # (quasi-)perfectly separating, which starves the softmax regression of
  # gradient and produces poor chain mixing. Also keep every mean well clear
  # of 0 so sampling noise can't produce a physically impossible negative value.
)
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

1. **Compositions**: π_k ~ Dirichlet(conc * **1**_S) for k = 1...K
2. **Membership**: P(*z*_i = k) = softmax(β_0k + β_1k * x_1i + ... + β_Pk * x_Pi), with β centred so the K coefficients sum to zero for each covariate (no reference community; see [Covariate coefficients have no reference community](#edna_dmm-fit-the-dmm))
3. **Counts**: **x**_i | *z*_i = k ~ DirichletMultinomial(N_i, α * π_k)

The global overdispersion α absorbs both technical (PCR, sequencing) and ecological compositional variance. Marginalizing over *z*_i makes inference exact.

---

## Frequently asked questions

**The chains disagree about which community is which. Is that a problem?**
No, and it is handled for you. A mixture likelihood is unchanged by permuting
the component labels, so chains routinely number the same communities
differently. `eDNA_dmm()` aligns every draw onto a common labelling before
summarising anything, so the estimates and diagnostics you see already
account for it.

**How do I know the chains found the same communities, not just the same names?**
Check `fit$alignment$min_agreement`, the share of samples the two least
similar chains place in the same community, and `fit$alignment$lp_rhat`, which
is invariant to labelling. Both near one means one solution found repeatedly.
If either departs clearly, the chains found genuinely different groupings and
relabelling cannot reconcile them; try other values of `K`. This cuts both
ways, it's not only a "K too large" symptom: K too small forces genuinely
distinct communities to merge, which is itself ambiguous about which samples
belong together and fails the same way. `eDNA_loo()` and
`eDNA_dmm_k_diagnostics()` compare a range of K at once rather than guessing
a direction.

**I have divergent transitions. What do I do?**
Increase `adapt_delta` toward `0.99`. If they persist, try lower K or verify your count matrix has no all-zero rows or columns.

**Can I use raw ASVs instead of taxonomy-collapsed counts?**
Yes. The model treats each column as a compositional unit and does not use taxonomy. ASVs give finer resolution; taxa collapse dimensionality and often converge faster.

**How do I include year as a covariate?**
Pass it as a numeric column. But if you have only a few discrete years, the linearity assumption may be too strong, consider fitting without year and testing it post-hoc via multinomial regression on the posterior assignments.

**The first run takes forever, is something wrong?**
No. Stan compiles the model to C++ on the first call after installation (1-2 minutes). All subsequent calls skip compilation. This is normal behavior for any rstan-based package.

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
