#load libraries
library(tidyverse)


#read in data
training_serology <- read.csv("/Users/ecramer/Desktop/Lover Lab /N-Lao-Serology/datasets/Epi_seropos_merge.csv") %>%
  select(!contains("_pos"))



#set training dataset of 1000 rows
set.seed(490)

training_serology  <- training_serology %>% 
  slice_sample(n = 1000)


write.csv(training_serology, "datasets/training_serology_sample.csv", row.names = FALSE)
