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

extract_loo_function_k2 <- function(x) {
## for extracting weibull params, K=2
fit_k2 <- mdl$sample(data = list(N = N, ## num observations
                                 y = log(x + 1),   
                                 K = 2),
                     chains=4,
                     parallel_chains=4)

draws_k2 <- as_draws_matrix(fit_k2$draws(c("scale", "shape", "theta"))) |> 
  data.frame()
mean_params_k2 <- colMeans(draws_k2)  

# for extracting log likelihood 
draws_k2_loglik <- as_draws_matrix(fit_k2$draws(c("log_lik")))
loo_k2 <- loo(draws_k2_loglik, save_psis = TRUE)
print(loo_k2)
}


## for extracting weibull params, K=3
extract_loo_function_k3 <- function(x)  {
fit_k3 <- mdl$sample(data = list(N = N, ## num observations
                                 y = log(x + 1),   
                                 K = 3),
                     chains=4,
                     parallel_chains=4)

# for extracting log likelihood 
draws_k3_loglik <- as_draws_matrix(fit_k3$draws(c("log_lik")))
loo_k3 <- loo(draws_k3_loglik, save_psis = TRUE)
print(loo_k3)
}

PvAMA1_k2_loo <- extract_loo_function_k2(entire_data$PvAMA1)
PvMSP119_k2_loo <- extract_loo_function_k2(entire_data$PvMSP119)
PfAMA1_k2_loo <- extract_loo_function_k2(entire_data$PfAMA1)
PfMSP119_k2_loo <- extract_loo_function_k2(entire_data$PfMSP119)

PvAMA1_k3_loo <- extract_loo_function_k3(entire_data$PvAMA1)
PvMSP119_k3_loo <- extract_loo_function_k3(entire_data$PvMSP119)
PfAMA1_k3_loo <- extract_loo_function_k3(entire_data$PfAMA1)
PfMSP119_k3_loo <- extract_loo_function_k3(entire_data$PfMSP119)






