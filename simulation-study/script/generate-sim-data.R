# generate simulated data

library(tidyverse)
library(sn)

source("simulation-study/script/annual_prob_infection.R")
source("simulation-study/script/sim_serodata_random_effect.R")

set.seed(123)

s <- annual_prob_infection(0.00055, 0.055, 2, 0.1, 70)
df <- read.csv("data/training_serology_sample.csv")
sample_size <- 1000
age_distribution <- table(df$AGE) / nrow(df)
times <- seq(1, 100, by = 1) # simulate yearly time steps across a 100 year period

for (i in 1:10){
sampled_ages <- sample(
  names(age_distribution),
  size = sample_size,
  replace = TRUE,
  prob = age_distribution
) %>%
  as.numeric()

birth_times <- max(times) - sampled_ages
d_sn <- sim_serodata_random_effect(sample_size, times, birth_times,
                                   prob_infection_fun = s,
                                   dist = "sn", 6, 1, 1, 2, 1.75, 4
)
save(d_sn, file = paste0("simulation-study/simulated-datasets/simulated_data-", i, ".rds"))
}

