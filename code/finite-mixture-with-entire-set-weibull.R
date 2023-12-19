library(cmdstanr)
library(posterior)
library(splines)
library(tidyverse)
library(loo)

#entire_data <- read.csv("datasets/training_serology_sample.csv")
entire_data <- read.csv("/Users/ecramer/Desktop/Lover Lab /N-Lao-Serology/datasets/Epi_seropos_merge.csv") %>%
  select(!contains("_pos")) %>% filter(AGE > 1)

mdl <- cmdstan_model(stan_file = "code/finite-mixture-weibull-loglik.stan")
N <- nrow(entire_data)     ## num of observations

## k = alpha = shape
## lambda = sigma = scale

## for extracting weibull params, K=2
fit_k2 <- mdl$sample(data = list(N = N, ## num observations
                                 y = log(entire_data$PfMSP119),   
                                 K = 2),
                     chains=4,
                     parallel_chains=4)

#shinystan::launch_shinystan(fit)

draws_k2 <- as_draws_matrix(fit_k2$draws(c("scale", "shape", "theta"))) |> 
  data.frame()
mean_params_k2 <- colMeans(draws_k2)  

hist(log(entire_data$PfMSP119), breaks = 100, freq=FALSE, main = "", xlab = "") 
dist1 <- function(x) dweibull(x, shape=mean_params_k2["shape.1."], scale=mean_params_k2["scale.1."])*mean_params_k2["theta.1."]
dist2 <- function(x) dweibull(x, shape=mean_params_k2["shape.2."], scale=mean_params_k2["scale.2."])*mean_params_k2["theta.2."]
curve(dist1(x) + dist2(x), add=TRUE)
curve(dist1(x), add=TRUE, col="blue")
curve(dist2(x), add=TRUE, col="red")

hist(log(entire_data$PfMSP119), breaks = 100, freq=FALSE, main = "Histogram of log(PfMSP119)\n k = 2",
     xlab = "log(PfMSP119)")
for(i in 1:500){
  dist1 <- function(x) dweibull(x, shape=draws_k2[i,"shape.1."], scale=draws_k2[i,"scale.1."])*draws_k2[i,"theta.1."]
  dist2 <- function(x) dweibull(x, shape=draws_k2[i,"shape.2."], scale=draws_k2[i,"scale.2."])*draws_k2[i,"theta.2."]
  curve(dist1(x) + dist2(x), add=TRUE, col=alpha("black", 0.05))
  curve(dist1(x), add=TRUE, col=alpha("blue", 0.05))
  curve(dist2(x), add=TRUE, col=alpha("red", 0.05))
}

## for extracting weibull params, K=3
fit_k3 <- mdl$sample(data = list(N = N, ## num observations
                                 y = log(entire_data$PfMSP119),   
                                 K = 3),
                     chains=4,
                     parallel_chains=4)

#shinystan::launch_shinystan(fit)

draws_k3 <- as_draws_matrix(fit_k3$draws(c("scale", "shape", "theta"))) |> 
  data.frame()

mean_params_k3 <- round(100*colMeans(draws_k3),2)
quantile_params_k3 <- 100*round(quantile(draws_k3$theta.1., probs = c(0.1, 0.9)),4)


hist(log(entire_data$PfMSP119), breaks = 100, freq=FALSE, main = "", xlab = "")
dist1 <- function(x) dweibull(x, shape=mean_params_k3["shape.1."], scale=mean_params_k3["scale.1."])*mean_params_k3["theta.1."]
dist2 <- function(x) dweibull(x, shape=mean_params_k3["shape.2."], scale=mean_params_k3["scale.2."])*mean_params_k3["theta.2."]
dist3 <- function(x) dweibull(x, shape=mean_params_k3["shape.3."], scale=mean_params_k3["scale.3."])*mean_params_k3["theta.3."]
curve(dist1(x) + dist2(x) + dist3(x), add=TRUE)
curve(dist1(x), add=TRUE, col="blue")
curve(dist2(x), add=TRUE, col="red")
curve(dist3(x), add=TRUE, col="green")

hist(log(entire_data$PfMSP119), breaks = 100, freq=FALSE, main = "Histogram of log(PfMSP119)\n k = 3", xlab = "log(PfMSP119)")
for(i in 1:500){
  dist1 <- function(x) dweibull(x, shape=draws_k3[i,"shape.1."], scale=draws_k3[i,"scale.1."])*draws_k3[i,"theta.1."]
  dist2 <- function(x) dweibull(x, shape=draws_k3[i,"shape.2."], scale=draws_k3[i,"scale.2."])*draws_k3[i,"theta.2."]
  dist3 <- function(x) dweibull(x, shape=draws_k3[i,"shape.3."], scale=draws_k3[i,"scale.3."])*draws_k3[i,"theta.3."]
  curve(dist1(x) + dist2(x) + dist3(x), add=TRUE, col=alpha("black", 0.05))
  curve(dist1(x), add=TRUE, col=alpha("blue", 0.05))
  curve(dist2(x), add=TRUE, col=alpha("red", 0.05))
  curve(dist3(x), add=TRUE, col=alpha("green", 0.05))
}


# for extracting log likelihood 
draws_k2_loglik <- as_draws_matrix(fit_k2$draws(c("log_lik")))
loo_k2 <- loo(draws_k2_loglik, save_psis = TRUE)
print(loo_k2)

# for extracting log likelihood 
draws_k3_loglik <- as_draws_matrix(fit_k3$draws(c("log_lik")))
loo_k3 <- loo(draws_k3_loglik, save_psis = TRUE)
print(loo_k3)


# for extracting draws 

# for extracting log likelihood from best fit (all markers)
probs_per_person <- as_draws_matrix(fit_k3$draws(c("log_prob_noinfection")))
avg_prob_no_infect <- colMeans(exp(probs_per_person))
hist(avg_prob_no_infect)


# # for extracting log likelihood from best fit (PfMSP119)
# probs_per_person <- as_draws_matrix(fit_k2$draws(c("log_prob_noinfection")))
# avg_prob_no_infect <- colMeans(exp(probs_per_person))
# hist(avg_prob_no_infect)

# 
# prob_noinfect_add_PvAMA1 <- cbind(entire_data, prob_noinfect_PvAMA1 = avg_prob_no_infect)
# 

# prob_noinfect_add_PvMSP119 <- cbind(prob_noinfect_add_PvAMA1, prob_noinfect_PvMSP119 = avg_prob_no_infect)


# prob_noinfect_add_PfAMA1 <- cbind(prob_noinfect_add_PvMSP119, prob_noinfect_PfAMA1 = avg_prob_no_infect)


prob_noinfect_all4_markers <- cbind(prob_noinfect_all4_markers, prob_noinfect_PfMSP119 = avg_prob_no_infect)

write.csv(prob_noinfect_all4_markers, "datasets/entire_data_marker_calls.csv")





