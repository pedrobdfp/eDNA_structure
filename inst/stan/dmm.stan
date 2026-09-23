// =============================================================================
// Dirichlet-Multinomial Mixture Model for eDNA Metabarcoding
// =============================================================================
// eDNAstructure package. Stan model.
//
// GENERATIVE MODEL:
//
//   Community compositions:
//     pi_k  ~ Dirichlet(conc * 1_S)              for k = 1..K
//
//   Community membership (covariate-driven softmax regression):
//     Gamma_ik = softmax(x_i' beta_k)
//
//   Observed counts (z marginalized out):
//     x_i ~ DirichletMultinomial(N_i, alpha * pi_k)   summed over k
//
//   Overdispersion:
//     alpha ~ gamma(alpha_shape, alpha_rate)
//       alpha >> 1 : near-multinomial (low overdispersion)
//       alpha ~  1 : high overdispersion (typical for eDNA)
//
// IDENTIFYING THE SOFTMAX:
//   Adding the same constant to every community's linear predictor leaves the
//   softmax unchanged, so for each covariate the K coefficients hold one
//   redundant degree of freedom that some constraint must remove.
//
//   This model removes it by CENTRING: for each covariate the K coefficients
//   sum to zero, so each is a deviation from the average community. The usual
//   alternative fixes one community at zero and reads every other coefficient
//   as a contrast against it.
//
//   Both are the same model in different coordinates. Centring keeps the prior
//   exchangeable across communities and leaves the posterior less correlated,
//   which samples better; fixing a reference gives one community a different
//   prior from the rest and makes every coefficient carry that community's
//   sampling noise. To read coefficients against a particular community
//   instead, subtract that community's row from the draws afterwards.
//
//   Consequently beta is K x (P+1) here, not (K-1) x (P+1).
//
// PARAMETERS TRACKED:
//   pi[K, S]              : community compositions over taxa
//   beta[K, P+1]          : centred softmax regression coefficients
//   alpha                 : global overdispersion scalar
//   community_probs[N, K] : posterior sample-level community membership probs
//   z_hat[N]              : sampled community assignment per sample
//   log_lik[N]            : per-sample log-likelihood (for LOO-CV)
//
// =============================================================================

data {
  int<lower=1> N;               // Number of samples
  int<lower=1> S;               // Number of taxa (or haplotypes)
  int<lower=2> K;               // Number of latent communities
  int<lower=0> P;               // Number of covariates (0 = intercept-only)

  array[N, S] int<lower=0> X;   // Count matrix [N x S], raw reads
  matrix[N, P] covariates;      // Covariate matrix [N x P], should be Z-scored

  real<lower=0> conc;           // Dirichlet concentration for pi priors
                                // conc < 1 = sparse communities (recommended for eDNA)
                                // conc = 1 = flat/uninformative
                                // conc > 1 = even compositions

  real<lower=0> alpha_shape;    // Gamma prior shape for overdispersion alpha
  real<lower=0> alpha_rate;     // Gamma prior rate  for overdispersion alpha
                                // Mean = shape/rate. Default: shape=5, rate=2 => mean=2.5
}

transformed data {
  // Per-sample read totals (N_i in the model)
  array[N] int<lower=0> Ni;
  for (i in 1:N) Ni[i] = sum(X[i]);              // total reads in sample i

  // Precompute log multinomial coefficients for efficiency. Constant in the
  // parameters, but kept so that log_lik is a genuine log density.
  vector[N] log_mc;
  for (i in 1:N) {
    real lmc = lgamma(Ni[i] + 1.0);              // log(N_i !)
    for (s in 1:S) lmc -= lgamma(X[i, s] + 1.0); // minus sum_s log(x_is !)
    log_mc[i] = lmc;
  }

  // Design matrix: intercept column prepended to covariates
  // When P=0 (intercept-only model), covariates is N x 0 matrix
  matrix[N, P + 1] X_design;
  if (P > 0) {
    X_design = append_col(rep_vector(1.0, N), covariates);  // [1 | covariates]
  } else {
    X_design = rep_matrix(1.0, N, 1);            // intercept only
  }
}

parameters {
  // K community compositions (K simplices over S taxa)
  array[K] simplex[S] pi;

  // Uncentred softmax coefficients [K x (P+1)]: all K communities, not K-1.
  // Each column's mean is redundant for the likelihood and is pinned by the
  // prior in the model block rather than by a constraint.
  matrix[K, P + 1] beta_raw;

  // Global overdispersion: single scalar shared across all communities
  real<lower=0> alpha;
}

transformed parameters {
  // Centred coefficients: within each covariate the K values sum to zero, so
  // each is that community's deviation from the average community.
  matrix[K, P + 1] beta;
  for (j in 1:(P + 1))
    beta[, j] = beta_raw[, j] - mean(beta_raw[, j]);  // subtract the column mean

  // Log-scale mixing weights [N x K], one linear predictor per sample and
  // community. No reference column is appended: all K are estimated.
  matrix[N, K] log_mixing_weights = X_design * beta';
}

model {
  // ── Priors ────────────────────────────────────────────────────────────────
  for (k in 1:K)
    pi[k] ~ dirichlet(rep_vector(conc, S));    // symmetric prior on each composition

  // Placed on the UNCENTRED coefficients, which gives the redundant column
  // means a proper prior and leaves the centred beta identified.
  to_vector(beta_raw) ~ normal(0, 1.0);

  alpha ~ gamma(alpha_shape, alpha_rate);      // overdispersion

  // ── Marginalized likelihood ───────────────────────────────────────────────
  for (i in 1:N) {
    vector[K] log_weights = log_softmax(log_mixing_weights[i]');  // log P(z_i = k)
    vector[K] lp;                                                 // log joint per community
    for (k in 1:K) {
      vector[S] alpha_k = alpha * pi[k];                          // DM concentration vector
      real dm = log_mc[i] + lgamma(alpha) - lgamma(Ni[i] + alpha);// DM normalising terms
      for (s in 1:S)
        dm += lgamma(X[i, s] + alpha_k[s]) - lgamma(alpha_k[s]);  // per-taxon DM terms
      lp[k] = log_weights[k] + dm;                                // membership x likelihood
    }
    target += log_sum_exp(lp);                                    // sum over communities
  }
}

generated quantities {
  // Posterior community membership probabilities per sample
  matrix[N, K] community_probs;
  // Sampled community assignment per sample
  array[N] int<lower=1, upper=K> z_hat;
  // Per-sample log-likelihood for LOO-CV
  vector[N] log_lik;

  for (i in 1:N) {
    vector[K] log_weights = log_softmax(log_mixing_weights[i]');  // log P(z_i = k)
    vector[K] lp;                                                 // log joint per community
    for (k in 1:K) {
      vector[S] alpha_k = alpha * pi[k];                          // DM concentration vector
      real dm = log_mc[i] + lgamma(alpha) - lgamma(Ni[i] + alpha);// DM normalising terms
      for (s in 1:S)
        dm += lgamma(X[i, s] + alpha_k[s]) - lgamma(alpha_k[s]);  // per-taxon DM terms
      lp[k] = log_weights[k] + dm;                                // membership x likelihood
    }
    log_lik[i]         = log_sum_exp(lp);       // marginal log density of sample i
    vector[K] post     = softmax(lp);           // P(z_i = k | data), normalised
    community_probs[i] = post';                 // stored as a row of the N x K matrix
    z_hat[i]           = categorical_rng(post); // one draw of the label, not the mode
  }
}
