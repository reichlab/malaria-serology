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
  int<lower=0> T;    // the year in which data was collected
  int year0;         // the year which should be treated as "timezero" in the logistic function 

}

transformed data{
  int<lower=0> maxage = max(age);
  array[maxage] int tvec = linspaced_int_array(maxage, T-year0-(maxage-1), T-year0);
}

// The parameters accepted by the model.
parameters {
  // parameters for never infected lognormal antibody distribution
  real mu0;
  real<lower=0> sigma0;
  // parameters for never infected lognormal antibody distribution
  real mu1;
  real<lower=0> sigma1;
  // parameters for the logistic curve for force of infection
  real logit_alpha_l;
  real logit_alpha_u;
  real<lower=0> beta;
  real<lower=0> neg_gamma;
}

transformed parameters{
  // precompute all needed values of the S function
  vector<lower=0,upper=1>[maxage] Svals = inv_logit(logit_alpha_l) + (inv_logit(logit_alpha_u)-inv_logit(logit_alpha_l))/(1 + beta*exp(neg_gamma*to_vector(tvec)));
  vector<lower=0,upper=1>[N] pi;
  for(n in 1:N) {
    // product involves computing the annual probabilities for the last age[n] values of the Svals vector
    pi[n] = 1 - prod(1-Svals[(maxage - age[n]+1):maxage]);
  }
}

// The model 
model {
  // priors
  mu0 ~ normal(log(3), 0.5);
  mu1 ~ normal(log(7), 0.5);
  sigma0 ~ lognormal(0, .1);
  sigma1 ~ lognormal(0, .1);
  logit_alpha_l ~ normal(-5, 2); // lower logistic asymptote on logit scale
  logit_alpha_u ~ normal(0, 2); // upper logistic asymptote on logit scale
  beta ~ lognormal(0, 2);
  neg_gamma ~ lognormal(0,1);

  for(n in 1:N) {
    target += log_mix(pi[n],
                      lognormal_lpdf(y[n] | mu0, sigma0),
                      lognormal_lpdf(y[n] | mu1, sigma1));
  };
}

