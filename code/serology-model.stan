//
// This Stan program defines a model that estimates two distributions
// of antibody levels for individuals ever and never infected with
// malaria. True infection status is not known, but is strongly associated
// with age, which is also observed. In addition to the two distributions,
// the model also estimates whether individuals were ever/never infected
// and what the temporal structure is of malaria risk over time.
//

// The input data 
data {
  int<lower=0> N;    // number of observations
  vector[N] y;       // observed levels of a particular marker
  array[N] int age;  // observed ages
  int<lower=0> maxage; // something that should be larger than the max observed age 
  matrix[maxage, 3] B; // design matrix with basis functions
}

transformed data{
  // int<lower=1> n_basis = 3; // number of basis functions, hard-coded
}

// The parameters accepted by the model.
parameters {
  // parameters for never infected lognormal antibody distribution
  real mu0;
  real<lower=0,upper=5> sigma0;
  // parameters for never infected lognormal antibody distribution
  real mu1;
  real<lower=0,upper=5> sigma1;
  
  // parameters for the smooth spline
  vector[3] beta; 
}

transformed parameters{
  // precompute all needed values of the S function
  vector<lower=0,upper=1>[maxage] Svals = inv_logit(B * beta);
  vector<lower=0,upper=1>[N] pi;
  for(n in 1:N) {
    // product involves computing the annual probabilities for the last age[n] values of the Svals vector
    pi[n] = 1 - prod(1-Svals[(maxage - age[n]+1):maxage]);
    
  }
}

// The model 
model {
  // priors
  mu0 ~ normal(log(3), 0.4); // log-scale mean of never infecteds
  mu1 ~ normal(log(7), 0.4); // log-scale mean of ever infecteds
  sigma0 ~ lognormal(0, .1);
  sigma1 ~ lognormal(0, .1);
  beta[1] ~ normal(1, 0.5);  // inv_logit(1) = 0.73
  beta[2] ~ normal(0, 0.5);  // inv_logit(0) = 0.5
  beta[3] ~ normal(-3, 0.5);  // inv_logit(-3) = 0.05

  for(n in 1:N) {
    target += log_mix(pi[n],
                      lognormal_lpdf(y[n] | mu1, sigma1),  // ever infected
                      lognormal_lpdf(y[n] | mu0, sigma0)); // never infected
  };
}

