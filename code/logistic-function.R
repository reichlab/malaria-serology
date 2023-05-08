
## logistic function

logisticfunc <- function(x, a_l, a_u, gamma, v, beta) a_l + (a_u-a_l)/(1+beta*exp(-gamma*x))^(1/v)
curve(logisticfunc(x, a_l=.1, gamma=-.1, a_u=.9, v=1, beta=5), 
      from=-70, to=70, ylim=c(0,1))


## priors for alphas

## lower asymptote
boot::logit(c(.001, .01, .3)) 
boot::inv.logit(summary(rnorm(1000, -5, 3)))


## upper asymptote
boot::logit(c(.6, .8, .9)) 
boot::inv.logit(summary(rnorm(1000, 0, 3)))

## beta, shift parameter
summary(rlnorm(10000, meanlog = 0, sdlog=2))
curve(tmp(x, a_l=.1, gamma=-.1, a_u=.9, v=1, beta=.02), 
      from=-70, to=70, ylim=c(0,1))

## gamma
qlnorm(c(0.01, .1, .25, .5, .75, .9, .99), meanlog = 0, sdlog=1)
curve(tmp(x, a_l=.1, gamma=-.01, a_u=.9, v=1, beta=1), 
      from=-70, to=70, ylim=c(0,1))
curve(tmp(x, a_l=.1, gamma=-10, a_u=.9, v=1, beta=1), 
      from=-70, to=70, ylim=c(0,1))


## year indices
T=2018
max_age = 100
year0 = 1990

## centers the s() function with zero at 99
t <- (T-year0-(max_age-1)):(T-year0)# -71:28
yr <- (T-(max_age-1)):T
idx <- 1:max_age
