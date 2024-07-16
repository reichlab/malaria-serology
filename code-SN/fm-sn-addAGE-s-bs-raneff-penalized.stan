/*
  Finite mixture model for Malaria serology data
  With changing probability of infection over time, s(t)
  Use splines for s(t) : logit(s(t)) = w0 + w * B
  where w0 is the intercept and w is the weight for the basis functions in B
  Note: w is scaled by tau, as a transformation for HMC
  Note: w_i's are penalized to reduce overfitting
  code reference: https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
*/


data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  vector<lower=0>[N]  y;   // observations
  array[N] int a;          // ages
  int num_basis;           // number of basis functions
  int num_times;           // number of time points
  matrix[num_basis, num_times] B;  // basis matrix
  }


// xi = location (real)
// omega = scale (positive, real)
// alpha = skewness (real)

parameters {
  ordered[K] location;    // locations of mixture components   
  vector<lower=0.00001>[K] scale;  // scales of mixture components
  vector[K] skewness;        // skewnesses of mixture components
  real<lower=0> mu; // mean of skew-normal distribution location for ever infected individuals (i.e., random effect)
  real<lower=0, upper=1> w0; // intercept for spline
  row_vector[num_basis] w_raw; // raw weights for basis functions in B
  real<lower=0> tau; 
}


// add transformed_parameters 
transformed parameters{
  row_vector[num_basis] w; // weights for basis functions in B
  vector<lower=0, upper=1>[100] s_hat; // probability of infection
  w[1] = w_raw[1];
  for (i in 2:num_basis)
        w[i] = w[i-1] + w_raw[i] * tau;
        s_hat = inv_logit(w0 + to_vector(w * B));
  
  vector<lower=0, upper=1>[N] theta; // probability of no infection by age a
  vector[N] logtheta; //log of theta
  for (n in 1:N) {
    //theta[n] = prod(1-s_hat[(100 - a[n]):100]); // direct calculation for theta using product
    logtheta[n] = sum(log1m(s_hat[(100 - a[n]):100]));
    theta[n] = exp(logtheta[n]);
  }
  
  vector[N] log_lik;
  vector[N] log_prob_noinfection;
  vector[K] lps;
     
  for (n in 1:N) {
    lps[1] = log(theta[n]);
    lps[2] = log1m(theta[n]);
    for (k in 1:K) {
      lps[k] += skew_normal_lpdf(y[n] | location[k], scale[k], skewness[k]); // skew-normal log pdf
    }
    log_lik[n] = log_sum_exp(lps); //store target for use in LOO IC
    log_prob_noinfection[n] = log(theta[n]) +  skew_normal_lpdf(y[n] | location[1], scale[1], skewness[1]) - log_lik[n];
  }
}


model {
  // priors
  location[1] ~ lognormal(1.5, 1);
  location[2] ~ normal(mu, 1);
  mu ~ lognormal(2.5, 1);
  scale ~ lognormal(1, 1);
  skewness ~ normal(0, 2);
 
  // spline priors
  w0 ~ beta(3, 30) ;
  tau ~ normal( 0 , 1 );
  
  target += sum(log_lik);
}

// add generated quantities
generated quantities {
  vector[N] y_tilde; // predicted data
  real<lower=0, upper=1> z; // status of infection
  for (n in 1:N) {
    z = bernoulli_rng(1-theta[n]); // status of infection
    if (z == 0) {
      y_tilde[n] = skew_normal_rng(location[1], scale[1], skewness[1]);
    } else {
      y_tilde[n] = skew_normal_rng(location[2], scale[2], skewness[2]);
    }
  }
}
 
 
