//
// This Stan program defines a finite mixture model
// code copied from:
// https://mc-stan.org/docs/stan-users-guide/summing-out-the-responsibility-parameter.html


data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  array[N] real y;         // observations
}

parameters {
  simplex[K] theta;          // mixing proportions
  ordered[K] scale;          // locations of mixture components
  vector<lower=0>[K] shape;  // scales of mixture components
}

model {
  vector[K] log_theta = log(theta);  // cache log calculation
  scale ~ normal(5, 3);
  shape ~ lognormal(0, 2);
  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += normal_lpdf(y[n] | scale[k], shape[k]);
    }
    target += log_sum_exp(lps);
  }
}
