# Simulate Lao malaria serology dataset
# Assumptions:
# 1) two distinct distributions of mixture for antibody levels
# 2) no Random effects
# log(PvAMA1) ~ either Normal or Skew Normal mixture distribution

# Parameters ----------------------------------------------------------------
# n: number of individuals
# times: number of time points
# birth_times: vector of birth times of n individuals
# infection_rate_fun: a smooth function of infection rates over time
# include_RE: logical, whether to include random effects (default FALSE)
# dist: distribution for the mixture for log-transformed antibody levels
# ...: additional arguments to the distribution:
# first set of arguments in the list is for infected people,
# second set of arguments in the list is for uninfected people
# For example,
# if dist = "norm", then ... = list(mean1 = 6, sd1 = 1, mean0 = 3, sd0 = 0.3)
# `mean1` indicates the mean of the normal distribution from which
# a mean is drawn for the infected individuals
# similarly, if dist = "sn", then
# ... = list(xi1 = 6, omega1 = 1, alpha1 = 3, xi0 = 3, omega0 = 0.3, alpha0 = 0.3)


sim_serodata <- function(n, times, birth_times,
                         prob_infection_fun,
                         include_randeff = FALSE, wanining = FALSE, 
                         dist = "norm", ...) {
  arg <- list(...)
  len_list <- length(arg)
  
  # base data frame with infection status
  base_dat <- data.frame(birth_times = birth_times) |>
    mutate(
      age = max(times) - birth_times
    ) |>
    rowwise() |>
    mutate(
      # calculate probability of never being infected by the age
      prob_never_infected = prod(1 - prob_infection_fun(max(times) - 0:age)),
      # infection status based on the probability of ever being infected
      # 1 indicates infected, 0 indicates not infected
      infection = rbinom(1, 1, 1 - prob_never_infected)
    ) |>
    # label infection status
    mutate(
      infection = factor(infection,
        levels = c(0, 1),
        labels = c("Not infected", "Infected")
      )
    ) |>
    ungroup()

  # determine infected_mean based on infection status and include_RE
  if (!include_randeff) {
    # without random effects
    dat <-  base_dat |>
      mutate(
        infected_mean = ifelse(infection == "Infected", arg[[1]],
          arg[[(len_list / 2) + 1]]
        )
      )
  } else if (include_randeff && !wanining) {
    # with random effects but no wanining
    dat <-  base_dat |>
      mutate(
        infected_mean = ifelse(infection == "Infected", rnorm(1, arg[[1]], 1),
          arg[[(len_list / 2) + 1]]
        )
      )
  } else if (include_randeff && wanining){
    # with random effects and wanining
    dat <- base_dat
  } else {
    stop("This is the setting not considered")
  }

  # log_PvAMA1 based on the choice of 'dist'
  if (dist == "norm") {
    df <- dat |>
      rowwise() |>
      mutate(
        log_PvAMA1 = ifelse(
          infection == "Infected",
          do.call(paste0("r", dist), list(n = 1, mean = infected_mean, sd = arg$sd1)),
          do.call(paste0("r", dist), list(n = 1, mean = infected_mean, sd = arg$sd0))
        )
      ) |>
      ungroup()
  } else {
    df <- dat |>
      rowwise() |>
      mutate(
        log_PvAMA1 = ifelse(
          infection == "Infected",
          do.call(
            paste0("r", dist),
            list(n = 1, xi = infected_mean, omega = arg$omega1, alpha = arg$alpha1)
          ),
          do.call(
            paste0("r", dist),
            list(n = 1, xi = infected_mean, omega = arg$omega0, alpha = arg$alpha0)
          )
        )
      ) |>
      ungroup()
  }
  out <- df |> mutate(PvAMA1 = exp(log_PvAMA1))

  return(out)
}
