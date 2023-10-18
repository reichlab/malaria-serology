#load libraries
library(tidyverse)


#load table for graph 
table_dat <- read.csv("datasets/seroneg-probabilities.csv")




boxplot <- ggplot(data = table_dat, aes(x = marker, y = seroneg_weibul, color = marker, group = marker)) + 
  geom_point(aes(shape = "Bayesian"), size = 3) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95)) +
  geom_point(aes(x = marker, y = seroneg_discrete, shape = "Discrete"), size = 3) +
  ylab("% seronegative") +
  theme_bw() +
  guides(color = 'none') + 
  scale_shape_discrete(name = "Model Classification") +
  xlab("Seromarker") + ylab("% Seronegative")


jpeg(file = "figures/seroneg_boxplot.jpg", width=8, height=6, units="in", res=200)
boxplot
dev.off()