library(cmdstanr)
library(posterior)
library(splines)
library(tidyverse)

training_data <- read.csv("datasets/training_serology_sample.csv")

mdl <- cmdstan_model(stan_file = "code/finite-mixture-weibull.stan")
N <- nrow(training_data)     ## num of observations

## k = alpha = shape
## lambda = sigma = scale

## for extracting weibull params, K=2
fit_k2 <- mdl$sample(data = list(N = N, ## num observations
                                 y = log(training_data$PvAMA1),   
                                 K = 2),
                     chains=4,
                     parallel_chains=4)

#shinystan::launch_shinystan(fit)

draws_k2 <- as_draws_matrix(fit_k2$draws(c("scale", "shape", "theta"))) |> 
  data.frame()
mean_params_k2 <- colMeans(draws_k2)  

hist(log(training_data$PvAMA1), breaks = 100, freq=FALSE)
dist1 <- function(x) dweibull(x, shape=mean_params_k2["shape.1."], scale=mean_params_k2["scale.1."])*mean_params_k2["theta.1."]
dist2 <- function(x) dweibull(x, shape=mean_params_k2["shape.2."], scale=mean_params_k2["scale.2."])*mean_params_k2["theta.2."]
curve(dist1(x) + dist2(x), add=TRUE)
curve(dist1(x), add=TRUE, col="blue")
curve(dist2(x), add=TRUE, col="red")

hist(log(training_data$PvAMA1), breaks = 100, freq=FALSE)
for(i in 1:500){
  dist1 <- function(x) dweibull(x, shape=draws_k2[i,"shape.1."], scale=draws_k2[i,"scale.1."])*draws_k2[i,"theta.1."]
  dist2 <- function(x) dweibull(x, shape=draws_k2[i,"shape.2."], scale=draws_k2[i,"scale.2."])*draws_k2[i,"theta.2."]
  curve(dist1(x) + dist2(x), add=TRUE, col=alpha("black", 0.05))
  curve(dist1(x), add=TRUE, col=alpha("blue", 0.05))
  curve(dist2(x), add=TRUE, col=alpha("red", 0.05))
}

## for extracting weibull params, K=3
fit_k3 <- mdl$sample(data = list(N = N, ## num observations
                                 y = log(training_data$PvAMA1),   
                                 K = 3),
                     chains=4,
                     parallel_chains=4)

#shinystan::launch_shinystan(fit)

draws_k3 <- as_draws_matrix(fit$draws(c("scale", "shape", "theta"))) |> 
  data.frame()
mean_params_k3 <- colMeans(draws_k3)  

hist(log(training_data$PvAMA1), breaks = 100, freq=FALSE)
dist1 <- function(x) dweibull(x, shape=mean_params_k3["shape.1."], scale=mean_params_k3["scale.1."])*mean_params_k3["theta.1."]
dist2 <- function(x) dweibull(x, shape=mean_params_k3["shape.2."], scale=mean_params_k3["scale.2."])*mean_params_k3["theta.2."]
dist3 <- function(x) dweibull(x, shape=mean_params_k3["shape.3."], scale=mean_params_k3["scale.3."])*mean_params_k3["theta.3."]
curve(dist1(x) + dist2(x) + dist3(x), add=TRUE)
curve(dist1(x), add=TRUE, col="blue")
curve(dist2(x), add=TRUE, col="red")
curve(dist3(x), add=TRUE, col="green")

hist(log(training_data$PvAMA1), breaks = 100, freq=FALSE)
for(i in 1:500){
  dist1 <- function(x) dweibull(x, shape=draws_k3[i,"shape.1."], scale=draws_k3[i,"scale.1."])*draws_k3[i,"theta.1."]
  dist2 <- function(x) dweibull(x, shape=draws_k3[i,"shape.2."], scale=draws_k3[i,"scale.2."])*draws_k3[i,"theta.2."]
  dist3 <- function(x) dweibull(x, shape=draws_k3[i,"shape.3."], scale=draws_k3[i,"scale.3."])*draws_k3[i,"theta.3."]
  curve(dist1(x) + dist2(x) + dist3(x), add=TRUE, col=alpha("black", 0.05))
  curve(dist1(x), add=TRUE, col=alpha("blue", 0.05))
  curve(dist2(x), add=TRUE, col=alpha("red", 0.05))
  curve(dist3(x), add=TRUE, col=alpha("green", 0.05))
}

