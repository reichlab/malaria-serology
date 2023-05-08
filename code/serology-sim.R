library(cmdstanr)
library(posterior)
library(splines)
library(tidyverse)

mdl <- cmdstan_model(stan_file = "code/serology-model.stan")

N <- 500       ## num of observations
age <- rnbinom(N, mu=30, size=5)
maxage <- max(age)+1
n_basis <- 3
time_vec <- seq.int(maxage)
B <- bs(time_vec, df = n_basis)

message(paste('iteration', i, Sys.time()))


## simulate data directly in R, couldn't get prior simulation to work natively in stan

## priors
mu0 <- rnorm(1, log(3), 0.2)
mu1 <- rnorm(1, log(7), 0.2)
# y1 <- exp(rnorm(1000, log(3),.2))
# y2 <- exp(rnorm(1000, log(7),.2))
# boxplot(y1, y2)
sigma0 <- rlnorm(1, 0, .1)
sigma1 <- rlnorm(1, 0, .1)
beta <- rnorm(3, c(1, 0, -3), 0.5)

Svals <- boot::inv.logit(B %*% beta)
## for checking structure of pis
# for(j in 1:N) {
#   ## product involves computing the annual probabilities for the last age[j] values of the Svals vector
#   pi[j] = 1 - prod(1-Svals[(maxage - age[j]+1):maxage])
# }

pi <- rep(NA, N)
for(j in 1:N) {
  pi[j] = 1 - prod(1-Svals[(maxage - age[j]+1):maxage])
}
z <- rbinom(N, 1, pi)
y <- (1-z)*rlnorm(N, mu0, sigma0) + z*rlnorm(N, mu1, sigma1)
plot(age, pi)

fit <- mdl$sample(data = list(N = N, ## num observations
                              y = y,   
                              age = age,  
                              maxage = maxage,
                              B = B),
                  chains=4,
                  parallel_chains=4)

pi_samples <- as_draws_matrix(fit$draws("pi")) |> 
  apply(MARGIN = 2, FUN=function(x) quantile(x, probs=c(0.1, 0.5, 0.9))) |> 
  t()
colnames(pi_samples) <- c("q10", "q50", "q90")
dat <- tibble(y, z, age, pi) |> 
  bind_cols(pi_samples)
ggplot(dat, aes(x=age)) +
  geom_point(aes(y=pi)) +
  geom_point(aes(y=q50), col="red", alpha=.1) +
  geom_segment(aes(y=q10, yend=q90, xend=age), alpha=.1, col="red")

mu0_samples <- as_draws_matrix(fit$draws("mu0")) |> 
  apply(MARGIN = 2, FUN=function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))|> 
  t()
colnames(mu0_samples) <- c("q10", "q50", "q90")

mu1_samples <- as_draws_matrix(fit$draws("mu1")) |> 
  apply(MARGIN = 2, FUN=function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))|> 
  t()
colnames(mu1_samples) <- c("q10", "q50", "q90")

mu_ests <- tibble(z = 0:1, mu = c(mu0, mu1)) |> 
  bind_cols(rbind(mu0_samples, mu1_samples))

ggplot(dat, aes(x=factor(z))) + 
  geom_boxplot(aes(y=y)) + 
  scale_y_log10() +
  geom_point(data=mu_ests, aes(y=exp(q50)), color="blue", shape=3, size=2) +
  geom_segment(data=mu_ests, aes(y=exp(q10), yend=exp(q90), xend=factor(z)), col="blue", size=2) +
  geom_point(data=mu_ests, aes(y=exp(mu)), color="red", shape=4) 


