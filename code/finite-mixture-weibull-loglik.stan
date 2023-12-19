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

// add transformed_parameters 

transformed parameters{
  vector[N] log_lik;
  vector[N] log_prob_noinfection; 
    vector[K] log_theta = log(theta);  // cache log calculation
     
   for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += weibull_lpdf(y[n] | shape[k], scale[k]);
    }
    log_lik[n] = log_sum_exp(lps); //store target for use in LOO IC
    log_prob_noinfection[n] = log_theta[1] +  weibull_lpdf(y[n] | shape[1], scale[1]) - log_lik[n];
    
  }
}

model {
  scale ~ lognormal(0, 2);
  shape ~ lognormal(0, 2);
  
  target += sum(log_lik);
 
}
