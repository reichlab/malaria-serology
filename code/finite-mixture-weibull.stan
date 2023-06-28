//
// This Stan program defines a finite mixture model
// code copied from:
// https://mc-stan.org/docs/stan-users-guide/summing-out-the-responsibility-parameter.html


data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  array[N] real y;         // observations
}

// k = alpha = shape
// lambda = sigma = scale
parameters {
  simplex[K] theta;          // mixing proportions
  ordered[K] scale;         // locations of mixture components
  vector<lower=0>[K] shape;  // scales of mixture components
}

model {
  vector[K] log_theta = log(theta);  // cache log calculation
  scale ~ lognormal(0, 2);
  shape ~ lognormal(0, 2);
  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += weibull_lpdf(y[n] | shape[k], scale[k]);
    }
    target += log_sum_exp(lps);
  }
}
