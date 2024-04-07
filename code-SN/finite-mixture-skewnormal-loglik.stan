//
// This Stan program defines a finite mixture model
// code copied from:
// https://mc-stan.org/docs/stan-users-guide/summing-out-the-responsibility-parameter.html


data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  array[N] real y;         // observations
}

// xi = location (real)
// omega = scale (positive, real)
// alpha = skewness (real)

parameters {
  simplex[K] theta;          // mixing proportions
  ordered[K] location;    // locations of mixture components   
  vector<lower=0.00001>[K] scale;  // scales of mixture components
  vector[K] skewness;        // skewnesses of mixture components
}

// add transformed_parameters 

transformed parameters{
  vector[N] log_lik;
  vector[N] log_prob_noinfection; 
  vector[K] log_theta = log(theta);  // cache log calculation
     
  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += skew_normal_lpdf(y[n] | location[k], scale[k], skewness[k]); // skew-normal log pdf
    }
    log_lik[n] = log_sum_exp(lps); //store target for use in LOO IC
    log_prob_noinfection[n] = log_theta[1] +  skew_normal_lpdf(y[n] | location[1], scale[1], skewness[1]) - log_lik[n];
  }
}


model {
  location ~ normal(0, 2);
  scale ~ lognormal(0, 2);
  skewness ~ normal(0, 2);
  
  target += sum(log_lik);
}
