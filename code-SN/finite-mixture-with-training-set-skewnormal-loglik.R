library(cmdstanr)
library(posterior)
library(splines)
library(tidyverse)
library(sn) #package for skew-normal dist.
library(loo)
options(mc.cores = 4)

training_data <- read.csv("data/training_serology_sample.csv")

mdl <- cmdstan_model(stan_file = "code-SN/finite-mixture-skewnormal-loglik.stan")
N <- nrow(training_data)     ## num of observations

## k = alpha = shape
## lambda = sigma = scale

## for extracting weibull params, K=2
fit_k2 <- mdl$sample(data = list(N = N, ## num observations
                                 y = log(training_data$PvAMA1),   
                                 K = 2),
                     chains=4,
                     parallel_chains=4) #, adapt_delta = 0.99

#shinystan::launch_shinystan(fit_k2)
#fit_k2$cmdstan_diagnose()

draws_k2 <- as_draws_matrix(fit_k2$draws(c("theta", "location", "scale", "skewness"))) |> 
  data.frame()
mean_params_k2 <- colMeans(draws_k2)  


hist(log(training_data$PvAMA1), breaks = 100, freq=FALSE)
dist1 <- function(x) dsn(x, xi=mean_params_k2["location.1."], omega = mean_params_k2["scale.1."], alpha = mean_params_k2["skewness.1."])*mean_params_k2["theta.1."]
dist2 <- function(x) dsn(x, xi=mean_params_k2["location.2."], omega = mean_params_k2["scale.2."], alpha = mean_params_k2["skewness.2."])*mean_params_k2["theta.2."]
curve(dist1(x) + dist2(x), add=TRUE)
curve(dist1(x), add=TRUE, col="blue")
curve(dist2(x), add=TRUE, col="red")



hist(log(training_data$PvAMA1), breaks = 100, freq=FALSE, main = "Fitting Skew Normal with k=2")
for(i in 1:500){
        dist_G0 <- function(x) dsn(x, xi = draws_k2[i, "location.1."], omega = draws_k2[i, "scale.1."], alpha = draws_k2[i,"skewness.1."])*draws_k2[i,"theta.1."]
        dist_G1 <- function(x) dsn(x, xi = draws_k2[i, "location.2."], omega = draws_k2[i, "scale.2."], alpha = draws_k2[i,"skewness.2."])*draws_k2[i,"theta.2."]
        curve(dist_G0(x) + dist_G1(x), add=TRUE, col=alpha("black", 0.05))
        curve(dist_G0(x), add=TRUE, col=alpha("blue", 0.05))
        curve(dist_G1(x), add=TRUE, col=alpha("red", 0.05))
}


## for extracting weibull params, K=3
fit_k3 <- mdl$sample(data = list(N = N, ## num observations
                                 y = log(training_data$PvAMA1), 
                                 K = 3),
                     chains=4,
                     parallel_chains=4)

#shinystan::launch_shinystan(fit_k3)

draws_k3 <- as_draws_matrix(fit_k3$draws(c("theta", "location", "scale", "skewness"))) |> 
        data.frame()
mean_params_k3 <- colMeans(draws_k3)  


hist(log(training_data$PvAMA1), breaks = 100, freq=FALSE)
dist1 <- function(x) dsn(x, xi=mean_params_k3["location.1."], omega = mean_params_k3["scale.1."], alpha = mean_params_k3["skewness.1."])*mean_params_k3["theta.1."]
dist2 <- function(x) dsn(x, xi=mean_params_k3["location.2."], omega = mean_params_k3["scale.2."], alpha = mean_params_k3["skewness.2."])*mean_params_k3["theta.2."]
dist3 <- function(x) dsn(x, xi=mean_params_k3["location.3."], omega = mean_params_k3["scale.3."], alpha = mean_params_k3["skewness.3."])*mean_params_k3["theta.3."]
curve(dist1(x) + dist2(x) + dist3(x), add=TRUE)
curve(dist1(x), add=TRUE, col="blue")
curve(dist2(x), add=TRUE, col="red")
curve(dist3(x), add=TRUE, col="green")


hist(log(training_data$PvAMA1), breaks = 100, freq=FALSE, main = "Fitting Skew Normal with k=3")
for(i in 1:1000){
        dist_G0 <- function(x) dsn(x, xi = draws_k3[i, "location.1."], omega = draws_k3[i, "scale.1."], alpha = draws_k3[i,"skewness.1."])*draws_k3[i,"theta.1."]
        dist_G1 <- function(x) dsn(x, xi = draws_k3[i, "location.2."], omega = draws_k3[i, "scale.2."], alpha = draws_k3[i,"skewness.2."])*draws_k3[i,"theta.2."]
        dist_G2 <- function(x) dsn(x, xi = draws_k3[i, "location.3."], omega = draws_k3[i, "scale.3."], alpha = draws_k3[i,"skewness.3."])*draws_k3[i,"theta.3."]
        curve(dist_G0(x) + dist_G1(x) + dist_G2(x), add=TRUE, col=alpha("black", 0.05))
        curve(dist_G0(x), add=TRUE, col=alpha("blue", 0.05))
        curve(dist_G1(x), add=TRUE, col=alpha("red", 0.05))
        curve(dist_G2(x), add=TRUE, col=alpha("green", 0.05))
}


# for extracting log likelihood 
draws_k2_loglik <- as_draws_matrix(fit_k2$draws(c("log_lik")))

loo_k2 <- loo(draws_k2_loglik, save_psis = TRUE)
print(loo_k2)

# for extracting log likelihood 
draws_k3_loglik <- as_draws_matrix(fit_k3$draws(c("log_lik")))

loo_k3 <- loo(draws_k3_loglik, save_psis = TRUE)
print(loo_k3)

loo_compare(loo_k2, loo_k3)
