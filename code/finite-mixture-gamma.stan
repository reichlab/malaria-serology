//
// This Stan program defines a finite mixture model
// code copied from:
// https://mc-stan.org/docs/stan-users-guide/summing-out-the-responsibility-parameter.html


data {
  int<lower=1> K;          // number of mixture components (columns in model matrix)
  int<lower=1> N;          // number of data points
  real y[N]; //the response
  matrix[N,K] X; //the model matrix
}

parameters {
 vector[K] betas; //the regression parameters
  real phi; //the variance parameter
}

transformed parameters {
  vector[N] mu; //the expected values (linear predictor)
  vector[N] alpha; //shape parameter for the gamma distribution
  vector[N] beta; //rate parameter for the gamma distribution
  
  mu = exp(X*betas); //using the log link 
  alpha = mu .* mu / phi; 
  beta = mu / phi;
}



model {  
  betas[1] ~ cauchy(0,10); //prior for the intercept following Gelman 2008

  for(i in 2:K)
   betas[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008
  
  y ~ gamma(alpha,beta);
}
generated quantities {
 vector[N] y_rep;
 for(n in 1:N){
  y_rep[n] = gamma_rng(alpha[n],beta[n]); //posterior draws to get posterior predictive checks
 }
}
