// BYM2 with a weighted ICAR component (Riebler et al. 2016 parameterisation).
// Gaussian outcome, county x wave panel. The ICAR graph changes by wave (one
// commuting vintage per wave), so the spatial component phi is wave-specific;
// the iid component theta is shared across waves, which keeps it separate
// from the residual sigma. Priors mirror the brms defaults used in
// 06_Models_Spatial.R so the two scripts are comparable.
//
// County-wave effects are stacked by wave: element (t - 1) * J + j.
data {
  int<lower=1> N;                      // observations (county x wave)
  int<lower=0> K;                      // predictors (centred, no intercept)
  matrix[N, K] X;
  vector[N] y;
  int<lower=1> J;                      // counties
  int<lower=1> T;                      // waves
  array[N] int<lower=1, upper=J * T> cw;   // stacked county-wave index
  int<lower=1> N_edges;                // edges summed over all wave graphs
  array[N_edges] int<lower=1, upper=J * T> node1;  // stacked indices
  array[N_edges] int<lower=1, upper=J * T> node2;
  vector<lower=0>[N_edges] weight;     // symmetric edge weights
  vector<lower=0>[T] scaling_factor;   // per wave: geometric mean of generalised-inverse diag(Q_t)
  real y_median;
}
parameters {
  real Intercept;
  vector[K] b;
  real<lower=0> sdcar;                 // SD of the combined county effect
  real<lower=0, upper=1> rhocar;       // share of variance that is spatial
  vector[J * T] phi;                   // weighted ICAR component, one field per wave
  vector[J] theta;                     // iid component, per county
  real<lower=0> sigma;
}
transformed parameters {
  vector[J * T] rcar;
  for (t in 1:T)
    rcar[((t - 1) * J + 1):(t * J)] =
      sdcar * (sqrt(rhocar / scaling_factor[t]) * phi[((t - 1) * J + 1):(t * J)]
               + sqrt(1 - rhocar) * theta);
}
model {
  // weighted ICAR: -0.5 * sum_ij w_ij (phi_i - phi_j)^2, soft sum-to-zero per wave
  target += -0.5 * dot_product(weight, square(phi[node1] - phi[node2]));
  for (t in 1:T)
    sum(phi[((t - 1) * J + 1):(t * J)]) ~ normal(0, 0.001 * J);
  theta ~ std_normal();

  Intercept ~ student_t(3, y_median, 2.5);
  b ~ normal(0, 1);
  sdcar ~ student_t(3, 0, 2.5);
  sigma ~ student_t(3, 0, 2.5);

  y ~ normal(Intercept + X * b + rcar[cw], sigma);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | Intercept + X[n] * b + rcar[cw[n]], sigma);
}
