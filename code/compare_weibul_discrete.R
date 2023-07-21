#load libraries
library(tidyverse)


#load table for graph 
table_dat <- read.csv("datasets/weibul_discrete_seroneg.csv")



ggplot(data = table_dat, aes(x = marker, y = seroneg_weibul, color = marker, group = marker)) + 
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95)) +
  geom_point(aes(x = marker, y = seroneg_discrete), shape = 17, size = 3) +
  ylab("% seronegative") +
  theme_minimal()