N <- 10000
age <- sample(1:80, size=N, replace=TRUE)
age_grp <- cut(age, breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80))

par(mfrow = c(4,1), mar=c(2,4,2,1))

## annual risk of infection (constant for now), with plot
p_infection <- 0.03
cum_risk_any_infection <- function(age, p_infection) 1-(1-p_infection)^age
curve(cum_risk_any_infection(x, p_infection = p_infection), 
      from=0, to=80, 
      ylim=c(0,1), ylab="prob ever infected")

## distributions of titers for infected vs. naives
meanlog_naive <- log(4)      ## 4 is median of never-infecteds
meanlog_infected <- log(6.5) ## 7 is median of never-infecteds
sdlog_naive <- 0.2
sdlog_infected <- 0.1

## plot the implied distributions
curve(dlnorm(x, meanlog=meanlog_naive, sdlog=sdlog_naive), from=0, to=15, ylab="density")
curve(dlnorm(x, meanlog=meanlog_infected, sdlog=sdlog_infected), from=0, to=15, ylab="density")

## compute the infections
p_any_infection <- 1-(1-p_infection)^age
any_infection <- rbinom(n=N, size=1, prob=p_any_infection)

## compute distribution of antibody levels
antibody_levels <- (1-any_infection) * rlnorm(N, meanlog = meanlog_naive, sdlog=sdlog_naive) +
  any_infection * rlnorm(N, meanlog = meanlog_infected, sdlog=sdlog_infected)

plot(age_grp, antibody_levels)
