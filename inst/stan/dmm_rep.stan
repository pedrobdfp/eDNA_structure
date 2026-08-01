// =============================================================================
// Replicate-aware Dirichlet-Multinomial Mixture Model for eDNA Metabarcoding
// eDNAstructure package — Stan model
// =============================================================================
//
// WHY THIS MODEL EXISTS
//
//   The standard model (dmm.stan) takes one row per station. When a station has
//   PCR / bottle replicates, the usual workaround is to sum their reads. That is
//   only lossless if replicates are plain multinomial draws from the station's
//   true composition — in that case their sum is a sufficient statistic, and
//   integrating the station composition out returns exactly dmm.stan.
//
//   Replicates carry extra information precisely because they are OVERdispersed
//   relative to multinomial (PCR jackpotting, uneven template, bottle effects).
//   Replicate concordance versus scatter is the signal, and summing destroys it.
//   This model keeps replicates as separate observations and adds a genuine
//   hierarchical layer for the station's own composition.
//
// GENERATIVE MODEL:
//
//   Community compositions:
//     pi_k    ~ Dirichlet(conc * 1_S)                  for k = 1..K
//
//   Community membership (covariate-driven softmax regression, station level):
//     Gamma_i = softmax(beta_0k + beta_1k*cov1_i + ...)
//     Community K is the reference category (linear predictor = 0)
//
//   Station composition (z marginalized out):
//     theta_i ~ Dirichlet(alpha * pi_k)                summed over k
//
//   Replicate counts:
//     y_ir    ~ DirichletMultinomial(N_ir, phi * theta_i)
//
//   Two dispersions, answering two different questions:
//     alpha — how tightly STATIONS cluster around their community composition
//             (same meaning as alpha in dmm.stan)
//     phi   — how tightly REPLICATES cluster around their own station
//             (new; large phi = highly reproducible replicates)
//
//   Because replicate counts are conditionally independent of k given theta_i,
//   the posterior community responsibilities depend only on the station layer.
//
// RAGGED REPLICATION:
//   Replicates enter as R rows tagged with a station index, so stations may have
//   different numbers of replicates. Stations with a single replicate are still
//   fully supported — their theta_i is informed by that replicate plus shrinkage
//   toward alpha * pi_k, using a phi learned from the replicated stations.
//   (If NO station has more than one replicate, alpha and phi are not separately
//   identified; the R wrapper detects that case and uses dmm.stan instead.)
//
// PARAMETERS TRACKED:
//   pi[K, S]              — posterior mean community compositions
//   theta[N, S]           — posterior mean STATION compositions
//   beta[K-1, P+1]        — softmax regression coefficients
//   alpha                 — station-level dispersion
//   phi                   — replicate-level dispersion
//   community_probs[N, K] — posterior station-level community membership probs
//   z_hat[N]              — MAP community assignment per station
//   log_lik[N]            — per-station log-likelihood contribution
//                           NOT suitable for LOO / K-selection: theta_i is a
//                           per-station parameter, so holding out a station
//                           leaves its own parameter unidentified. Use dmm.stan
//                           (via eDNA_loo()) for K selection.
//
// =============================================================================

data {
  int<lower=1> N;               // Number of stations (ecological units)
  int<lower=1> R;               // Number of replicate rows (total observations)
  int<lower=1> S;               // Number of taxa (or haplotypes)
  int<lower=2> K;               // Number of latent communities
  int<lower=0> P;               // Number of covariates (0 = intercept-only)

  array[R, S] int<lower=0> X;              // Replicate-level counts [R x S]
  array[R] int<lower=1, upper=N> station;  // Station index of each replicate row
  matrix[N, P] covariates;                 // STATION-level covariates [N x P]

  real<lower=0> conc;           // Dirichlet concentration for pi priors

  real<lower=0> alpha_shape;    // Gamma prior on station-level dispersion
  real<lower=0> alpha_rate;     // Mean = shape/rate. Default: 5, 2 => 2.5

  real<lower=0> phi_shape;      // Gamma prior on replicate-level dispersion
  real<lower=0> phi_rate;       // Default: 2, 0.02 => mean 100 (replicates are
                                // typically much tighter than stations)
}

transformed data {
  // Per-replicate read totals
  array[R] int<lower=0> Nr;
  for (r in 1:R) Nr[r] = sum(X[r]);

  // Precompute log multinomial coefficients
  vector[R] log_mc;
  for (r in 1:R) {
    real lmc = lgamma(Nr[r] + 1.0);
    for (s in 1:S) lmc -= lgamma(X[r, s] + 1.0);
    log_mc[r] = lmc;
  }

  // Counts as real vectors so lgamma() can be applied vectorised in the model
  // block rather than looping over taxa on every leapfrog step.
  array[R] vector[S] X_real;
  for (r in 1:R)
    for (s in 1:S)
      X_real[r, s] = X[r, s];

  // Design matrix: intercept column prepended to station-level covariates
  matrix[N, P + 1] X_design;
  if (P > 0) {
    X_design = append_col(rep_vector(1.0, N), covariates);
  } else {
    X_design = rep_matrix(1.0, N, 1);
  }
}

parameters {
  // K community compositions (K simplices over S taxa)
  array[K] simplex[S] pi;

  // N station compositions — the "true composition" of each sample, estimated
  // jointly from that station's replicates and the community it belongs to.
  array[N] simplex[S] theta;

  // Softmax regression coefficients [K-1 x (P+1)]; community K is the reference
  matrix[K - 1, P + 1] beta;

  real<lower=0> alpha;          // station-around-community dispersion
  real<lower=0> phi;            // replicate-around-station dispersion
}

transformed parameters {
  // Log-scale mixing weights [N x K] (community K = reference = 0)
  matrix[N, K] log_mixing_weights;
  {
    if (K > 1) {
      matrix[N, K - 1] eta = X_design * beta';
      log_mixing_weights = append_col(eta, rep_vector(0.0, N));
    } else {
      log_mixing_weights = rep_matrix(0.0, N, 1);
    }
  }
}

model {
  // Declarations first, so this compiles on older Stan versions too.
  array[K] vector[S] a_pi;        // alpha * pi_k
  vector[K] dir_norm;             // Dirichlet normalising constant per community
  array[N] vector[S] a_theta;     // phi * theta_i
  array[N] vector[S] lg_a_theta;  // lgamma(phi * theta_i)

  // ── Priors ────────────────────────────────────────────────────────────────
  for (k in 1:K)
    pi[k] ~ dirichlet(rep_vector(conc, S));

  to_vector(beta) ~ normal(0, 1.0);

  alpha ~ gamma(alpha_shape, alpha_rate);
  phi   ~ gamma(phi_shape,   phi_rate);

  // ── Precompute community-level quantities (independent of station) ────────
  for (k in 1:K) {
    a_pi[k]     = alpha * pi[k];
    dir_norm[k] = lgamma(alpha) - sum(lgamma(a_pi[k]));
  }

  // ── Station layer: theta_i ~ Dirichlet(alpha * pi_k), k marginalized ──────
  // Written as normalising constant + dot product rather than K*N calls to
  // dirichlet_lpdf(), which is substantially cheaper.
  for (i in 1:N) {
    vector[K] log_weights = log_softmax(log_mixing_weights[i]');
    vector[S] log_theta_i = log(theta[i]);
    vector[K] lp;
    for (k in 1:K)
      lp[k] = log_weights[k] + dir_norm[k] + dot_product(a_pi[k] - 1.0, log_theta_i);
    target += log_sum_exp(lp);
  }

  // ── Replicate layer: y_ir ~ DirichletMultinomial(N_ir, phi * theta_i) ─────
  // Hoist the per-station concentration out of the replicate loop; with several
  // replicates per station this avoids recomputing lgamma() for each of them.
  for (i in 1:N) {
    a_theta[i]    = phi * theta[i];
    lg_a_theta[i] = lgamma(a_theta[i]);
  }
  for (r in 1:R) {
    int i = station[r];
    target += log_mc[r] + lgamma(phi) - lgamma(Nr[r] + phi)
              + sum(lgamma(X_real[r] + a_theta[i]) - lg_a_theta[i]);
  }
}

generated quantities {
  // Posterior community membership probabilities per station
  matrix[N, K] community_probs;
  // MAP community assignment per station
  array[N] int<lower=1, upper=K> z_hat;
  // Per-station log-likelihood contribution.
  // NOTE: not valid for LOO / K-selection — see header.
  vector[N] log_lik;

  {
    array[K] vector[S] a_pi;
    vector[K] dir_norm;
    vector[N] rep_ll = rep_vector(0.0, N);

    for (k in 1:K) {
      a_pi[k]     = alpha * pi[k];
      dir_norm[k] = lgamma(alpha) - sum(lgamma(a_pi[k]));
    }

    // Replicate-layer contribution, accumulated onto each station
    for (r in 1:R) {
      int i = station[r];
      vector[S] a_theta_i = phi * theta[i];
      rep_ll[i] += log_mc[r] + lgamma(phi) - lgamma(Nr[r] + phi)
                   + sum(lgamma(X_real[r] + a_theta_i) - lgamma(a_theta_i));
    }

    for (i in 1:N) {
      vector[K] log_weights = log_softmax(log_mixing_weights[i]');
      vector[S] log_theta_i = log(theta[i]);
      vector[K] lp;
      vector[K] post;
      for (k in 1:K)
        lp[k] = log_weights[k] + dir_norm[k] + dot_product(a_pi[k] - 1.0, log_theta_i);

      post               = softmax(lp);
      community_probs[i] = post';
      z_hat[i]           = categorical_rng(post);
      // Station-layer marginal + this station's replicate contribution
      log_lik[i]         = log_sum_exp(lp) + rep_ll[i];
    }
  }
}
