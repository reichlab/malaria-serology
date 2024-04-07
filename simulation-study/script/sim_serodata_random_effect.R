# Simulate Lao malaria serology dataset
# considering random effects on the infected individuals
# log(PvAMA1) ~ either Normal or Skew Normal mixture distribution

# Parameters ----------------------------------------------------------------
# n: number of individuals
# times: number of time points
# birth_times: vector of birth times of n individuals
# infection_rate_fun: a smooth function of infection rates over time
# dist: distribution for the mixture for antibody levels
# ...: additional arguments to the distribution: first list for infected, second list for not infected
# for example, if dist = "norm", then ... = list(mean = 6, sd = 1, mean = 3, sd = 0.3)
# if dist = "sn", then ... = list(xi = 6, omega = 1, alpha = 3, xi = 3, omega = 0.3, alpha = 0.3)
# the first argument in the list indicates the mean of the normal distribution
# from which a mean is drawn for the infected individuals

sim_serodata_random_effect <- function(n, times, birth_times,
                                       prob_infection_fun,
                                       dist = "norm", ...) {
  arg <- list(...)
  len_list <- length(arg)
  
  # set.seed(123)
  dat <- data.frame(birth_times = birth_times) %>%
    mutate(
      age = max(times) - birth_times
    ) %>%
    rowwise() %>%
    mutate(
      prob_never_infected = prod(1 - prob_infection_fun(max(times) - 0:age))
    ) %>%
    mutate(
      infection = rbinom(1, 1, 1 - prob_never_infected)
    ) %>%
    mutate(
      infection = factor(infection,
                         levels = c(0, 1),
                         labels = c("Not infected", "Infected")
      )
    ) %>%
    mutate(
      infected_mean = ifelse(
        infection == "Infected",
        rnorm(1, arg[[1]], 1),
        arg[[(len_list / 2) + 1]]
      )
    ) %>%
    ungroup()
  
  if (dist == "norm") {
    dat <- dat %>%
      rowwise() %>%
      mutate(
        log_PvAMA1 = ifelse(
          infection == "Infected",
          do.call(paste0("r", dist), c(list(n = 1, mean = infected_mean), arg[2:(len_list / 2)])),
          do.call(paste0("r", dist), c(list(n = 1), arg[(len_list / 2) + 1:(len_list / 2)]))
        )
      ) %>%
      ungroup()
  } else {
    dat <- dat %>%
      rowwise() %>%
      mutate(
        log_PvAMA1 = ifelse(
          infection == "Infected",
          do.call(paste0("r", dist), c(list(n = 1, xi = infected_mean), arg[2:(len_list / 2)])),
          do.call(paste0("r", dist), c(list(n = 1), arg[(len_list / 2) + 1:(len_list / 2)]))
        )
      ) %>%
      ungroup()
  }
  out <- dat %>% mutate(PvAMA1 = exp(log_PvAMA1))
  
  return(out)
}



