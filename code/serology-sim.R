library(cmdstanr)
library(posterior)

mdl <- cmdstan_model(stan_file = "code/serology-model.stan")

nsim <- 10

N <- 1000       ## num of observations
age <- rnbinom(N, mu=35, size=5)
maxage <- max(age)
T <- 2016
year0 <- 1990
tvec <- (T-year0-(maxage-1)):(T-year0)
yearvec <- year0+tvec

for(i in 1:nsim){
  message(paste('iteration', i, Sys.time()))
  
  ## simulate data directly in R, couldn't get prior simulation to work natively in stan
  
  ## priors
  mu0 <- rnorm(1, log(3), .5)
  mu1 <- rnorm(1, log(7), .5)
  sigma0 <- rlnorm(1, 0, .1)
  sigma1 <- rlnorm(1, 0, .1)
  logit_alpha_l <- rnorm(1, -5, 3)
  logit_alpha_u<- rnorm(1, 0, 3)
  beta <- rlnorm(1, 0, 2)
  neg_gamma <- rlnorm(1, 0,1)
  
  Svals <- boot::inv.logit(logit_alpha_l) + 
    (boot::inv.logit(logit_alpha_u)-boot::inv.logit(logit_alpha_l))/(1 + beta*exp(neg_gamma*tvec))
  oneminusSvals = 1-Svals
  curve(logisticfunc(x, a_l = boot::inv.logit(logit_alpha_l), a_u = boot::inv.logit(logit_alpha_u), gamma = -neg_gamma, v=1, beta = beta), from=min(tvec), to=max(tvec))
  
  for(j in 1:N) {
    pi[j] = 1 - prod(oneminusSvals[(maxage - age[j]+1):maxage])
  }
  z <- rbinom(N, 1, pi)
  y <- (1-z)*rlnorm(N, mu0, sigma0) + z*rlnorm(N, mu1, sigma1)
  
  fit <- mdl$sample(data = list(N = N, ## num observations
                                y = y,   
                                age = age,  
                                T = T,
                                year0 = year0),
                    chains=4,
                    parallel_chains=4)
  
  pi_samples <- as_draws_df(fit$draws("pi")) |> 
    apply(MARGIN = 2, FUN=function(x) quantile(x, probs=c(0.1, 0.9)))
  for(k in 1:K){
    coverage_pi[i,k] <- probs[k]>pi_samples[1,k] & probs[k]<pi_samples[2,k]
  }
  
  beta_samples <- as_draws_df(fit$draws("beta")) |> 
    apply(MARGIN = 2, FUN=function(x) quantile(x, probs=c(0.1, 0.9)))
  for(k in 1:K){
    coverage_beta[i,k] <- beta[k]>beta_samples[1,k] & beta[k]<beta_samples[2,k]
  }
}
