library(cmdstanr)
library(posterior)
library(splines)
library(tidyverse)
library(sn) #package for skew-normal dist.
library(bayesplot) # for ppc_dens_overlay
library(loo)
options(mc.cores = 4)

sim_data <- readRDS("./simulation-study/simulated-datasets/simulated_data-1.rds")

N <- nrow(sim_data)
time_points <- 1:100
num_times <- length(time_points)
num_knots <- 5
# knot_list <- quantile(time_points, probs = seq(0, 1, length.out = num_knots)) # equally spaced knots
knot_list <- quantile(time_points, probs = c(0, 0.4, 0.6, 0.8, 1)) 

# Create B-spline basis matrix
B <- t(bs(time_points, 
          knots = knot_list[-c(1, num_knots)], 
          degree = 3, 
          intercept = TRUE))
num_basis <- nrow(B)

# Fit the model with transformed parameters of B-spline with penalized priors
mdl_age_unfixedprob_bs <- cmdstan_model(stan_file = "./code-SN/fm-sn-addAGE-s-bs-raneff-penalized.stan")

fit_age_unfixedprob_bs <- mdl_age_unfixedprob_bs$sample(
  data = list(N = N,
              y = log(sim_data$PvAMA1),
              K = 2,
              a = sim_data$age,
              num_basis = num_basis,
              num_times = num_times,
              B = B),
  chains=4,
  parallel_chains=4,
  iter_warmup = 1000,
  iter_sampling = 3000)

#shinystan::launch_shinystan(fit_age_unfixedprob_bs)
#fit_age_unfixedprob_bs$cmdstan_diagnose()


# Posterior predictive checks 
y <- log(sim_data$PvAMA1)
# extract the samples from the posterior predictive distribution and compare these to the observed data
yrep <- as_draws_matrix(fit_age_unfixedprob_bs$draws(c("y_tilde")))
rnd_idx <- sample(nrow(yrep), 100)
yrep_sample <- yrep[rnd_idx, ]

ppc_dens_overlay(y, yrep_sample) +
  scale_x_continuous(limits = c(-2, 13))


# Posterior draws of $s(t)$
draws_s_hat <- as_draws_matrix(fit_age_unfixedprob_bs$draws("s_hat")) |> 
  data.frame()
rnd_idx_s <- sample(nrow(draws_s_hat), 100) 
draws_s_hat_sample <- draws_s_hat[rnd_idx_s,]

plot(t <- 1:100, 0.00055 + (0.055-0.00055)/(1+2*exp(-0.1*(70-t))), 
     ylim = c(0, 0.15),
     col="red", type='l', lwd=3,
     main = "Posterior draws of s, changing probability of infection over time",
     xlab = "Time", ylab = "Annual Probability of Infection")
for (i in 1:nrow(draws_s_hat_sample)) {
  s_hat_sample <- t(draws_s_hat_sample[i,]) %>% as.data.frame()
  colnames(s_hat_sample) <- "s_hat"
  ss <- s_hat_sample %>% mutate(time = 1:100)
  lines(ss$time, ss$s_hat, lwd=2, add=TRUE, col=alpha("lightblue", 0.4))
}