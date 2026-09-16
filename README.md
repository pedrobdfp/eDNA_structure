# eDNAstructure

**Dirichlet-Multinomial Mixture Models for eDNA Metabarcoding Community Structure**

`eDNAstructure` is an R package for fitting Bayesian Dirichlet-Multinomial Mixture (DMM) models to environmental DNA (eDNA) read count data from metabarcoding surveys. Given a sample × taxon count matrix and optional environmental covariates, the model identifies latent ecological communities, estimates their taxonomic compositions, and quantifies how environmental gradients drive community membership — all within a fully Bayesian framework with principled uncertainty quantification.

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
remotes::install_github("pedrobdfp/eDNA_structure", upgrade = "never")
```

During installation, a large amount of black text will appear — this is the Stan model compiling to C++. It only happens once. Every subsequent call to `eDNA_dmm()` goes straight to sampling with no compilation output.

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
library(dplyr)    # for pipe and data manipulation
library(ggplot2)  # for plot customization

# Option A — use the built-in example dataset
data <- get_example_data()

# Option B — simulate your own dataset with known ground truth
data <- simulate_eDNA_survey(
  n_communities         = 4,
  n_species             = 40,
  samples_per_community = 5,
  seed                  = 2026
)

# Inspect raw species composition before fitting
plot_true_compositions(
  data$counts,
  metadata  = data$covariates,
  facet_var = "TrueCommunity"
)

# Fit the model with a given number of communities (K)
fit <- eDNA_dmm(
  counts     = data$counts,
  covariates = data$covariates[, c("Depth", "Distance_shore")],
  K          = 4
)

print(fit)
summary(fit)

# Have PCR / bottle replicates? Don't sum them — pass each replicate as its own
# row and group them with station_id. See "Replicated samples" under Functions.
# fit <- eDNA_dmm(rep_counts, station_id = station_id,
#                 covariates = station_covs, K = 4)

# Or select K using LOO cross-validation
loo_result <- eDNA_loo(data$counts, data$covariates[, c("Depth", "Distance_shore")],
                       K_range = 2:5)
loo_result$plot


# The loo_result object stores all fitted models — no need to refit
# Extract the K=4 model directly:
fit <- loo_result$fits[["K4"]]

# Or using the K value as a number:
K_best <- 4
fit <- loo_result$fits[[paste0("K", K_best)]]

# Confirm what you have:
print(fit)

##You can also plot the results!

# Structure plot — one bar per sample, colored by community membership probability
eDNA_dmm_structure(fit, metadata = data$covariates,
                   facet_var = "TrueCommunity", sort_var = "Depth")

# NMDS ordination colored by community assignment
eDNA_dmm_nmds(fit)$plot

# Prior vs posterior distributions for covariate effects
eDNA_dmm_beta(fit)$plot
```

> **For a complete walkthrough** — including step-by-step simulation, data formatting, K selection, all visualization options, parameter recovery, and troubleshooting — see the **[full tutorial vignette](vignettes/tutorial.Rmd)**. It is designed to be read start to finish and assumes no prior familiarity with Bayesian mixture models.

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

A **sample × covariate** data frame in the same row order as the count matrix. Covariates are Z-score standardized internally by default.

```r
head(data$covariates)
#   sample_id  TrueCommunity  Depth  Distance_shore
# 1   STN_001              1     82             198
# 2   STN_002              1     79             204
# 3   STN_003              2     11             197
```

### Replicates (optional)

If a station was sequenced more than once (PCR or bottle replicates), give each replicate its own row and use `station_id` to say which rows belong together. Do **not** sum them — see [Replicated samples](#replicated-samples--station_id) below for why.

```r
rep_counts[1:4, 1:4]
#                 Sp_1  Sp_2  Sp_3  Sp_4
# STN_001_B1       402   315   118    70
# STN_001_B2       421   298   131    77
# STN_002_B1       380   270   101    58
# STN_002_B2       395   281    95    64

station_id <- c("STN_001", "STN_001", "STN_002", "STN_002")
```

Covariates stay **station-level**: supply either one row per station (in order of first appearance in `station_id`) or one row per row of `counts`, in which case they are collapsed automatically. Replication may be ragged — stations can have different numbers of replicates, and stations with a single replicate are fully supported.

Leave `station_id` unset if each row is an independent sample. That is the standard model and the right choice for unreplicated data.

---

## Functions

### `eDNA_dmm()` — Fit the DMM

The core function. Fits a Dirichlet-Multinomial Mixture model via Stan and returns an `edna_dmm_fit` object.

```r
fit <- eDNA_dmm(
  counts           = my_counts,   # sample × taxon integer count matrix
  covariates       = my_covs,     # sample × covariate data frame, or NULL
  K                = 4,           # number of latent communities to fit
  station_id       = NULL,        # group replicate rows into stations; NULL = each row
                                  # is an independent sample. See "Replicated samples" below
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
  phi_shape        = 2,           # Gamma prior shape for replicate dispersion phi
  phi_rate         = 0.02,        # Gamma prior rate  (prior mean = 100)
                                  # phi_* are used only when station_id is supplied
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
| `theta_mean` | Matrix [N × S]: posterior mean station compositions (replicate model only, else `NULL`) |
| `phi_mean` | Scalar: posterior mean replicate dispersion (replicate model only, else `NULL`) |
| `stan_fit` | Raw `rstan::stanfit` object for advanced diagnostics |

> **On single chains:** Mixture models suffer from label switching across chains — "Community 1" in chain A may map to "Community 2" in chain B, making multi-chain Rhat diagnostics meaningless. A single long chain sidesteps this. Use within-chain ESS (reported by `summary()`) as your convergence criterion.

---

### Replicated samples — `station_id`

If your stations have PCR or bottle replicates, pass each replicate as its own row of `counts` and use `station_id` to group them. This switches `eDNA_dmm()` to a hierarchical model that estimates each station's own composition:

```
theta_i ~ Dirichlet(alpha * pi_k)                      # the station's true composition
y_ir    ~ DirichletMultinomial(N_ir, phi * theta_i)    # each replicate
```

There are then two dispersion parameters answering two different questions:

| Parameter | Meaning |
|-----------|---------|
| `alpha` | How tightly **stations** cluster around their community composition (same as in the standard model) |
| `phi` | How tightly **replicates** cluster around their own station — replicate reproducibility |

```r
fit <- eDNA_dmm(
  counts     = rep_counts,        # one row per replicate
  station_id = station_id,        # which rows belong to the same station
  covariates = station_covs,      # station-level covariates
  K          = 4,
  phi_shape  = 2,                 # Gamma prior on phi (default mean 100)
  phi_rate   = 0.02
)

fit$theta_mean   # [N stations × S taxa] posterior mean station compositions
fit$phi_mean     # posterior mean replicate reproducibility
```

**Why not just sum the replicates?** Summing would be lossless if replicates were plain multinomial draws from the station's composition — their sum is a sufficient statistic, and integrating the station composition out gives you back exactly the standard model. Replicates are informative *precisely because* they are overdispersed relative to multinomial (PCR jackpotting, uneven template, bottle effects). Replicate concordance versus scatter is the signal, and summing discards it.

Practical notes:

- **Ragged replication is fine.** Stations may have different numbers of replicates. A station with only one replicate is still fitted normally — its `theta_i` is informed by that replicate plus shrinkage toward `alpha * pi_k`, using a `phi` learned from the replicated stations.
- **Unreplicated data is unaffected.** If no station has more than one replicate, `alpha` and `phi` are not separately identified and the replicate model would buy you nothing, so `eDNA_dmm()` uses the standard model instead and tells you. Leaving `station_id` unset is always the right choice when each row is an independent sample.
- **`sample_info` is still one row per station**, so `eDNA_dmm_structure()`, `eDNA_dmm_nmds()` and `eDNA_dmm_beta()` all work unchanged.
- **K selection stays on the standard model.** `eDNA_loo()` deliberately uses the summed model: with a per-station `theta_i`, holding out a station leaves its own parameter unidentified, so station-level LOO is not well defined. The replicate model is slower too, which matters when sweeping many K values.

---

### `eDNA_dmm_structure()` — Structure bar plot

Produces a STRUCTURE-style plot: one vertical bar per sample, divided into colored segments by posterior community membership probability.

```r
p <- eDNA_dmm_structure(
  fit,
  metadata         = my_metadata,     # data frame with additional sample variables
  sample_id_col    = "sample_id",     # column in metadata matching sample IDs in fit
  facet_var        = "year",          # column facets
  facet_row_var    = "depth_bin",     # row facets — see below
  sort_var         = "Depth",         # sort samples within each panel
  community_colors = NULL,            # named hex vector, e.g. c("Community 1" = "#E63946")
                                      # or NULL for automatic HCL palette
  bar_width        = 0.9,             # bar width (0–1); 1 = no gaps
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

`eDNA_dmm_structure()` shares its entire layout parameter surface with `plot_true_compositions()` — by design, so the two plots can be stacked into one figure and line up exactly. All of these behave identically in both functions:

| Group | Parameters |
|---|---|
| Two-way faceting | `facet_var`, `facet_row_var` |
| Row labels | `row_label_fn`, `row_label_position`, `row_label_size`, `panel_labels` |
| Legend layout | `legend_nrow`, `legend_ncol`, `legend_key_size`, `legend_text_size`, `legend_rel_height`, `legend_rel_width` |
| Axis / strip sizing | `ylab`, `ylab_size`, `ylab_rel_width`, `strip_text_size`, `axis_text_size` |
| Reference lines | `vline_var`, `vline_value`, `vline_color`, `vline_linetype`, `vline_linewidth` |

See [`plot_true_compositions()`](#plot_true_compositions--raw-species-composition) above for the full description of each group. Two differences:

- `ylab` defaults to `"Membership probability"` here (vs `"Proportion"`).
- `axis_text_size` (y-axis tick label size) exists here only. Default `NULL` = ggplot2's own scaling from `base_size`.

#### Stacking the observed and fitted plots

The common case — raw composition on top, community assignment below, sharing panel structure:

```r
common <- list(
  metadata      = meta,
  facet_var     = "year",
  facet_row_var = "depth_bin",
  sort_var      = "Depth",
  row_label_fn  = function(x) paste0(x, " m"),
  base_size     = 11
)

p_obs <- do.call(plot_true_compositions,
                 c(list(counts, taxa_include = dmm_taxon_order(fit, 20)), common))
p_fit <- do.call(eDNA_dmm_structure,
                 c(list(fit, panel_labels = FALSE), common))

cowplot::plot_grid(p_obs, p_fit, ncol = 1, rel_heights = c(1, 1))
```

Set `panel_labels = FALSE` on the inner plots when the outer figure supplies its own (a), (b) letters — otherwise you get two competing sets.

---


### `eDNA_dmm_nmds()` — NMDS ordination

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
                                 # 50% certainty (smallest) → 100% certainty (largest)
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

Visualizes observed species frequencies per sample as stacked bars — the same layout as `eDNA_dmm_structure()`, allowing direct before/after comparison. Most useful before fitting to inspect the raw community signal, and with simulated data where true community labels are known.

```r
p <- plot_true_compositions(
  counts,
  metadata        = my_metadata,
  sample_id_col   = "sample_id",
  facet_var       = "year",          # column facets
  facet_row_var   = "depth_bin",     # row facets — see "Two-way layouts" below
  sort_var        = "Depth",         # sort samples within each panel
  top_n           = 20,              # top N taxa individually; rest → "Other"
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

#### Two-way layouts — `facet_row_var`

`facet_var` makes panel *columns*; `facet_row_var` makes panel *rows*. With `facet_row_var` set, one sub-plot is built per row level and they are stacked with cowplot. That architecture exists for a reason: it is the only one that keeps panel widths proportional to sample count (via `space = "free_x"`), so a stratum with 30 samples is drawn three times wider than one with 10 instead of being stretched to match.

This is the recommended layout for publication figures with two-way structure (e.g. depth stratum × year).

#### Row labels — `row_label_fn`, `row_label_position`

| Parameter | Purpose |
|---|---|
| `row_label_fn` | Function applied to each `facet_row_var` level to build its display label. `function(x) paste0(x, " m")` for depth, `function(x) paste0("Year: ", x)` for year. Default `as.character`. |
| `row_label_position` | `"title"` (bold label above each row's panel — **new default**) or `"ylab"` (as the y-axis label — the old, hardcoded behavior). |
| `row_label_size` | Font size of that label. Default `NULL` = `base_size`. |
| `panel_labels` | Label each row panel with a lowercase letter (a, b, c…). Default `TRUE`. Set `FALSE` when this plot is itself a panel inside a larger composite figure, where the outer figure supplies the letters and an inner set would collide. |

`row_label_fn` is what makes the label text generalize across datasets — it is tied to whatever `facet_row_var` means for your data, rather than being hardcoded.

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

Taxa are sorted and laid out hue block by hue block, so setting `legend_ncol = n_hues` puts roughly one color family per column — which is what makes a 20-taxon legend scannable.

#### Axis and strip sizing

| Parameter | Purpose |
|---|---|
| `ylab` | Shared y-axis label drawn once beside the panel stack. Used only when `row_label_position = "title"`. Default `"Proportion"`. |
| `ylab_size` | Font size of that label. Default `NULL` = `base_size`. |
| `ylab_rel_width` | Width of the label strip as a fraction of the stack. Default `0.03`. |
| `strip_text_size` | Font size of the column-facet strip labels. Default `NULL` = `base_size`. |

#### Reference lines — `vline_*`

Draws a vertical reference line at a threshold in sample-ordering space — e.g. marking where depth crosses 200 m once samples are sorted by depth.

```r
p <- plot_true_compositions(
  counts, metadata = meta, sort_var = "Depth",
  vline_var       = "Depth",     # numeric column to locate the threshold in
  vline_value     = 200,         # the threshold
  vline_color     = "black",
  vline_linetype  = "dashed",
  vline_linewidth = 0.7
)
```

#### Sharing a palette with the fitted plot — `taxa_include`

```r
taxa <- dmm_taxon_order(fit, top_n = 20)
plot_true_compositions(counts, taxa_include = taxa)
```

`taxa_include` names exactly which taxa to show individually; everything else is pooled into `"Other"`, and `top_n` is ignored. Names absent from `counts` are dropped with a warning. See [`dmm_taxon_order()`](#dmm_taxon_order--shared-taxon-ranking-for-matching-palettes) for why this matters.

---

---
### `eDNA_dmm_compositions()` — Posterior community compositions

Visualizes the posterior mean taxonomic composition of each latent community as stacked bars — the model's estimate of what each community "looks like" in species space. Colors match those used in `eDNA_dmm_structure()` and `plot_true_compositions()` for direct comparison.

```r
p <- eDNA_dmm_compositions(
  fit,
  top_n           = 20,      # show top N taxa; rest collapsed to "Other"
  base_size       = 13,
  title           = NULL,
  subtitle        = NULL,
  legend_position = "right", # "right", "bottom", "left", "top", or "none"
  bar_width       = 0.7      # bar width (0–1)
)
```

Returns a `ggplot2` object. The x-axis labels show community numbers (1, 2, 3…). Pair with `eDNA_dmm_structure()` to connect community identities to sample assignments.

---

### `dmm_taxon_order()` — Shared taxon ranking (for matching palettes)

Returns taxa ranked by total posterior weight across communities (`colSums(pi_mean)`). A taxon that dominates a single community ranks highly even if it is rare overall — the ordering that matters when the question is what *distinguishes* communities.

```r
taxa <- dmm_taxon_order(
  fit,
  top_n = 20        # NULL (default) returns every taxon, ranked
)
```

**Use this whenever two composition figures sit side by side.** The taxon palette is a deterministic function of the sorted taxon set: same set in, same colors out. If `eDNA_dmm_compositions()` picks its own top 20 and `plot_true_compositions()` independently picks its own, the two sets differ and shared taxa get *different colors in each panel*. Driving both from one ranking fixes that:

```r
taxa <- dmm_taxon_order(fit, top_n = 20)

p_obs   <- plot_true_compositions(counts, taxa_include = taxa)   # observed
p_model <- eDNA_dmm_compositions(fit, top_n = 20)                # fitted

cowplot::plot_grid(p_obs, p_model, ncol = 1)
```

Note this is deliberately **not** the same ranking as mean observed frequency across samples, which rewards being widespread rather than being diagnostic.

---


### `get_example_data()` — Built-in example dataset

Returns the built-in simulated dataset: 20 samples across 4 communities separated by depth and distance from shore, generated by `simulate_eDNA_survey()` with known ground truth, so fitted parameters can be compared to the true values.

```r
data <- get_example_data()
# data$counts                 — integer matrix, 20 samples × 33 taxa
# data$covariates             — data frame: sample_id, SampleID, TrueCommunity,
#                               Depth, Distance_shore
# data$community_compositions — true composition matrix, 4 × 40
# data$metab_df               — raw simulated metabarcoding reads
# data$sample_metadata        — full simulation metadata
# data$contributors           — per-sample organism lists behind the reads
```

> **Why 33 taxa and not 40?** The simulation draws 40 species, but taxa that end up
> with zero reads across every sample are dropped from `counts` (they carry no
> information and the model rejects all-zero columns). `community_compositions`
> keeps the full 40-column ground truth, so the two deliberately differ in width.

---

### `eDNA_clear_stan_cache()` — Reset the compiled model cache

The Stan model is compiled once on first use and cached per machine under
`tools::R_user_dir("eDNAstructure", "cache")`. Clear it to force a fresh compile —
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
  community_covariate_means = NULL,        # K × P matrix of covariate means per community
                                           # NULL = default 2-covariate depth × shore design
  covariate_sds             = NULL,        # K × P SDs (NULL = 5 for all)
  mean_read_depth           = 10000,       # mean reads per sample
  bio_reps                  = 1,           # biological replicates per station
  seq_reps                  = 1,           # sequencing technical replicates per bio rep
  spillover                 = 0.15,        # fraction of species frequency leaking between communities
  shedding_error            = 0.3,         # lognormal SD on per-organism eDNA shedding
  decay_rate                = 0.1,         # exponential distance decay of eDNA signal
  seed                      = 42
)
# sim$counts      — ready for eDNA_dmm()
# sim$covariates  — ready for eDNA_dmm()
```

Or run each step individually for full control:

```r
community_mat   <- generate_community_compositions(
  n_communities     = 4,
  n_species         = 40,
  n_dominant_range  = c(2, 3),     # 2–3 high-frequency dominant species per community
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
  rep             = 1        # sequencing technical replicates per bottle
)

sample_metadata <- generate_sample_covariates(
  contributors_list    = contrib_obj$contributors_list,
  community_covariates = matrix(           # community-specific covariate means
    c(80, 200, 10, 200, 80, 20, 10, 20),
    nrow = 4, byrow = TRUE,
    dimnames = list(NULL, c("Depth", "Distance_shore"))
  )
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

1. **Compositions**: π_k ~ Dirichlet(conc · **1**_S) for k = 1…K
2. **Membership**: P(*z*_i = k) = softmax(β_0k + β_1k · x_1i + … + β_Pk · x_Pi), community K = reference
3. **Counts**: **x**_i | *z*_i = k ~ DirichletMultinomial(N_i, α · π_k)

The global overdispersion α absorbs both technical (PCR, sequencing) and ecological compositional variance. Marginalizing over *z*_i makes inference exact.

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
No. Stan compiles the model to C++ on the first call after installation (1–2 minutes). All subsequent calls skip compilation. This is normal behavior for any rstan-based package.

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