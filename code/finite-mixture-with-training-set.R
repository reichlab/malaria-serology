library(cmdstanr)
library(posterior)
library(splines)
library(tidyverse)

training_data <- read.csv("datasets/training_serology_sample.csv")

mdl <- cmdstan_model(stan_file = "code/finite-mixture-gamma.stan")
N <- nrow(training_data)     ## num of observations
K = 2 #number of categories (clusters)

X<-model.matrix(~x1*x2,dat)


fit <- mdl$sample(data = list(N = N, ## num observations
                              y = log(training_data$PvAMA1),   
                              K = K),
                  chains=4,
                  parallel_chains=4)

shinystan::launch_shinystan(fit)


hist(log(training_data$PvAMA1), breaks = 100, freq = FALSE)
curve(dnorm(x, 3.3, 0.5)*.3, add=TRUE)
curve(dnorm(x, 4.9, 2)*.7, add=TRUE)


# model with age 



mdl2 <- cmdstan_model(stan_file = "code/serology-model.stan")

n_basis <- 3
maxage = max(training_data$AGE) + 1
time_vec <- seq.int(maxage)
B <- bs(time_vec, df = n_basis)


fit2 <- mdl2$sample(data = list(N = N, ## num observations
                              y = log(training_data$PvAMA1),   
                              age = training_data$AGE,
                              maxage = maxage,
                              B = B),
                   
                  chains=4,
                  parallel_chains=4)

shinystan::launch_shinystan(fit2)



## summarizing data and pi estimates
pi_samples <- as_draws_matrix(fit2$draws("pi")) |>
  apply(MARGIN = 2, FUN=function(x) quantile(x, probs=c(0.1, 0.5, 0.9))) |> t()
colnames(pi_samples) <- c("q10", "q50", "q90")
dat <- tibble(y = log(training_data$PvAMA1),   
              age = training_data$AGE) |>
  bind_cols(pi_samples)

## summarizing posteriors of the group means
mu0_samples <- as_draws_matrix(fit2$draws("mu0")) |>
  apply(MARGIN = 2, FUN=function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))|> t()

colnames(mu0_samples) <- c("q10", "q50", "q90")
mu1_samples <- as_draws_matrix(fit2$draws("mu1")) |>
  apply(MARGIN = 2, FUN=function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))|>
  t()
colnames(mu1_samples) <- c("q10", "q50", "q90")
mu_ests <- tibble(z = 0:1, mu = c(mu0, mu1)) |>
  bind_cols(rbind(mu0_samples, mu1_samples))


ggplot(dat, aes(x=age)) +
  geom_point(aes(y=q50), col="red", alpha=.1) +
  geom_segment(aes(y=q10, yend=q90, xend=age), alpha=.1, col="red")

## plotting svals from the FOI model
Svals_samples <- as_draws_matrix(fit2$draws("Svals")) |> 
  apply(MARGIN = 2, FUN=function(x) quantile(x, probs=c(0.1, 0.5, 0.9))) |> 
  t() 
colnames(Svals_samples) <- c("q10", "q50", "q90")
Svals_qs <- tibble(yr= 1:nrow(Svals_samples)) |> 
  bind_cols(Svals_samples)

ggplot(Svals_qs, aes(x=yr)) +
  geom_point(aes(y=q50)) +
  geom_segment(aes(y=q10, yend=q90, xend=yr)) +
  ylim(0,1)