// Non-spatial counterpart to bym2_weighted.stan: same Gaussian likelihood,
// design matrix, and priors, with the county effect either dropped (pooled,
// the Bayesian analogue of OLS) or iid (BYM2 with rho fixed at 0).
data {
  int<lower=1> N;                      // observations (county x wave)
  int<lower=0> K;                      // predictors (centred, no intercept)
  matrix[N, K] X;
  vector[N] y;
  int<lower=1> J;                      // counties
  array[N] int<lower=1, upper=J> county;
  int<lower=0, upper=1> has_county;    // 1 = iid county effect, 0 = pooled
  real y_median;
}
parameters {
  real Intercept;
  vector[K] b;
  array[has_county] real<lower=0> sd_county;  // SD of the iid county effect
  vector[has_county ? J : 0] z;              // standardised county effects
  real<lower=0> sigma;
}
transformed parameters {
  vector[J] r = rep_vector(0, J);
  if (has_county)
    r = sd_county[1] * z;
}
model {
  z ~ std_normal();

  Intercept ~ student_t(3, y_median, 2.5);
  b ~ normal(0, 1);
  sd_county ~ student_t(3, 0, 2.5);
  sigma ~ student_t(3, 0, 2.5);

  y ~ normal(Intercept + X * b + r[county], sigma);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | Intercept + X[n] * b + r[county[n]], sigma);
}
