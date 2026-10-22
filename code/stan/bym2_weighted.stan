// BYM2 with a weighted ICAR component (Riebler et al. 2016 parameterisation).
// Gaussian outcome, repeated observations per area. Priors mirror the brms
// defaults used in 06_Models_Spatial.R so the two scripts are comparable.
data {
  int<lower=1> N;                      // observations (county x wave)
  int<lower=0> K;                      // predictors (centred, no intercept)
  matrix[N, K] X;
  vector[N] y;
  int<lower=1> J;                      // counties
  array[N] int<lower=1, upper=J> county;
  int<lower=1> N_edges;
  array[N_edges] int<lower=1, upper=J> node1;
  array[N_edges] int<lower=1, upper=J> node2;
  vector<lower=0>[N_edges] weight;     // symmetric edge weights
  real<lower=0> scaling_factor;        // geometric mean of generalised-inverse diag(Q)
  real y_median;
}
parameters {
  real Intercept;
  vector[K] b;
  real<lower=0> sdcar;                 // SD of the combined county effect
  real<lower=0, upper=1> rhocar;       // share of variance that is spatial
  vector[J] phi;                       // weighted ICAR component
  vector[J] theta;                     // iid component
  real<lower=0> sigma;
}
transformed parameters {
  vector[J] rcar = sdcar * (sqrt(rhocar / scaling_factor) * phi + sqrt(1 - rhocar) * theta);
}
model {
  // weighted ICAR: -0.5 * sum_ij w_ij (phi_i - phi_j)^2, soft sum-to-zero
  target += -0.5 * dot_product(weight, square(phi[node1] - phi[node2]));
  sum(phi) ~ normal(0, 0.001 * J);
  theta ~ std_normal();

  Intercept ~ student_t(3, y_median, 2.5);
  b ~ normal(0, 1);
  sdcar ~ student_t(3, 0, 2.5);
  sigma ~ student_t(3, 0, 2.5);

  y ~ normal(Intercept + X * b + rcar[county], sigma);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | Intercept + X[n] * b + rcar[county[n]], sigma);
}
