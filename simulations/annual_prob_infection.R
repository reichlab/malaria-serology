# Annual probability of infection

# Parameters
# lower: the infection rate at the end of the time period
# upper: the infection rate at the beginning of the time period
# beta: a parameter that shifts the curve's transition earlier or later in time
# gamma: a parameter that governs the steepness of the transition
# t0: sigmoid midpoint
annual_prob_infection <- function(lower, upper, beta, gamma, t0) {
  s <- function(t) lower + (upper - lower) / (1 + beta * exp(-gamma * (-t + t0)))
  return(s)
}

# cumulative probability of being never infected by age
cum_prob_never_infected <- function(age, max_time, infection_fun) {
  return(prod(1 - infection_fun(max_time - 0:age)))
}


##-- example
##-- s(t) = alpha_l + (alpha_u - alpha_l) / (1 + beta * exp(-gamma * (-t + t0)))

# alpha_l <- mean_params_k2["alpha_l"]
# alpha_u <- mean_params_k2["alpha_u"]
# beta <- mean_params_k2["beta"]
# gamma <- mean_params_k2["gamma"]
# t0 <- mean_params_k2["t0"]
# s <- annual_prob_infection(lower=alpha_l, upper=alpha_u, beta=beta, gamma=gamma, t0=t0)
# t <- 1:100
# plot(t, s(t), type = "l", ylab = "probability of infection")
# cum_prob_never_infected(80, 100, s) # should be 0.08-0.1
