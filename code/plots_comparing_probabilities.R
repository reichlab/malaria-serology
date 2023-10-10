# Analysis with probabilities for each person from the training set 


library(tidyverse)
library(betareg)

calls_per_person <- read.csv("datasets/training_data_marker_calls.csv")




#graph of Pv data 
ggplot(data = calls_per_person, aes(x = AGE, y = prob_noinfect_PvAMA1)) +
  geom_point() + 
  stat_smooth(method ="loess", se = FALSE, color = "black") + 
  geom_point(aes(y = prob_noinfect_PvMSP119), color = "blue") +
  stat_smooth(aes(y = prob_noinfect_PvMSP119), method ="loess", se = FALSE, color = "blue") 


#graph of Pf data 
ggplot(data = calls_per_person, aes(x = AGE, y = prob_noinfect_PfAMA1)) +
  geom_point() + 
  stat_smooth(method ="loess", se = FALSE, color = "black") + 
  geom_point(aes(
    y = prob_noinfect_PfMSP119), color = "blue") +
  stat_smooth(aes(y = prob_noinfect_PfMSP119), method ="loess", se = FALSE, color = "blue") 




ggplot(data = calls_per_person, aes(x = prob_noinfect_PvAMA1, y = prob_noinfect_PvMSP119, color = AGE)) + 
  geom_point()



ggplot(data = calls_per_person, aes(x = prob_noinfect_PfAMA1, y = prob_noinfect_PfMSP119, color = AGE)) + 
  geom_point()





#### Beta Regression 

library(betareg) #install.packages("betareg)



raw_demo_data <- read.csv("Datasets/Parademo_individual_master_v15a.csv") %>% 
  select(id_individual, id_HH, district_id, pcr_result, pcr_species, BARCODE = barcode, sex = q4, occupation = q7a, village_id, hohh_id, agebands, fever, ITN_coverage_bin) %>%
  mutate(ITN_coverage_bin = fct_relevel(ITN_coverage_bin, ref = "No"),
         sex = fct_relevel(sex, ref = "Female"),
         district_id = fct_relevel(district_id, ref = "KH")) %>%
  mutate("District" = fct_recode(factor(district_id),
                                 "Paktha "= "PA", 
                                 "Nambak" = "NB",
                                 "Khua" = "KH", 
                                 "Muang Et" = "ET")) %>%
  mutate(occupation_agg = fct_recode(factor(occupation),
                          `Unknown / Unemployed` = "",
                          `Unknown / Unemployed` = "Not working",
                          `Small-scale farmer` = "Small-scale farmer",
                          `Farm / hunt` = "Plantation work (for cash)",
                          `Farm / hunt` = "Collecting plants/hunting in forest",
                          `Student / Teacher` = "Student",
                          `Student / Teacher` = "Teacher",
                          `Manual labor` = "Manual labor",
                          `Small business owner` = "Small business owner"),
         occupation_agg = fct_relevel(occupation_agg, ref = "Small-scale farmer"))
  


merge_data <- calls_per_person %>%
  left_join(raw_demo_data)


m1_all <- betareg(prob_noinfect_PvAMA1 ~ AGE + district_id + SEX + ITN_coverage_bin + occupation_agg, data = merge_data)
summary(m1_all)

m2_rITN <- betareg(prob_noinfect_PvAMA1 ~ AGE + district_id + SEX + occupation_agg, data = merge_data)
summary(m2_rITN)

m3_rOccupation <- betareg(prob_noinfect_PvAMA1 ~ AGE + district_id + SEX, data = merge_data)
summary(m3_rOccupation)

m4_rSex <- betareg(prob_noinfect_PvAMA1 ~ AGE + district_id, data = merge_data)
summary(m4_rSex)


AIC(m1_all)
AIC(m2_rITN)
AIC(m3_rOccupation)
AIC(m4_rSex)


BIC(m1_all)
BIC(m2_rITN)
BIC(m3_rOccupation)
BIC(m4_rSex)

plot(m1_all)
plot(m2_rITN)
plot(m3_rOccupation)




ggplot(merge_data, aes(x = AGE, y = prob_noinfect_PvAMA1 )) +
  geom_point(size = 4, aes(fill = district_id), shape = 21) +
  scale_fill_grey() + 
  geom_line(aes(y = predict(m1_all, merge_data), 
                colour = "m1_all", linetype = "m1_all")) +
  geom_line(aes(y = predict(m2_rITN, merge_data), 
                colour = "m2_rITN", linetype = "m2_rITN")) +
  geom_line(aes(y = predict(m3_rOccupation, merge_data), 
                colour = "m3_rOccupation", linetype = "m3_rOccupation")) +
  geom_line(aes(y = predict(m4_rSex, merge_data),
                colour = "m4_rSex", linetype = "m4_rSex")) +
  theme_bw()









merge_data <- merge_data %>%
  mutate(prob_noinfect_PfMSP119 = prob_noinfect_PfMSP119 + 0.00001)

m1_all <- betareg(prob_noinfect_PfMSP119 ~ AGE + district_id + SEX + ITN_coverage_bin + occupation_agg, data = merge_data)
summary(m1_all)

m2_rITN <- betareg(prob_noinfect_PfMSP119 ~ AGE + district_id + SEX + occupation_agg, data = merge_data)
summary(m2_rITN)

m3_rOccupation <- betareg(prob_noinfect_Pfama1_edit ~ AGE + district_id + SEX, data = merge_data_pfama1)
summary(m3_rOccupation)

m4_rSex <- betareg(prob_noinfect_PfMSP119 ~ AGE + district_id, data = merge_data)
summary(m4_rSex)


merge_data_pfmsp119 <- merge_data %>% 
  mutate(prob_noinfect_PfMSP119_edit = prob_noinfect_PfMSP119 + 0.000001)
m0_age <- betareg(prob_noinfect_PfMSP119_edit ~ AGE, data = merge_data_pfmsp119)
summary(m0_age)


merge_data_pfama1 <- merge_data %>% 
  mutate(prob_noinfect_Pfama1_edit = prob_noinfect_PfAMA1 + 0.00001)
m0_age <- betareg(prob_noinfect_Pfama1_edit ~ AGE, data = merge_data_pfama1)
summary(m0_age)



AIC(m1_all)
AIC(m2_rITN)
AIC(m3_rOccupation)
AIC(m4_rSex)


BIC(m1_all)
BIC(m2_rITN)
BIC(m3_rOccupation)
BIC(m4_rSex)

plot(m1_all)
plot(m2_rITN)
plot(m3_rOccupation)


ggplot(merge_data, aes(x = AGE, y = prob_noinfect_PfMSP119)) +
  geom_point(size = 4, aes(fill = district_id), shape = 21) +
  scale_fill_grey() + 
  # geom_line(aes(y = predict(m1_all, merge_data), 
  #               colour = "m1_all", linetype = "m1_all")) +
  geom_line(aes(y = predict(m2_rITN, merge_data), 
                colour = "m2_rITN", linetype = "m2_rITN")) +
  geom_line(aes(y = predict(m3_rOccupation, merge_data), 
                colour = "m3_rOccupation", linetype = "m3_rOccupation")) +
  geom_line(aes(y = predict(m4_rSex, merge_data),
                colour = "m4_rSex", linetype = "m4_rSex")) +
  theme_bw() 



ggplot(merge_data, aes(x = AGE, y = prob_noinfect_PfAMA1)) +
  geom_point(size = 4, shape = 21) +
  scale_fill_grey() + 
  ylab("Probability of no infection") + 
  xlab("Age") +
  geom_line(aes(y = predict(m0_age, merge_data), 
                color = "Beta-Regression model")) +
  theme_bw() 




#distribution plot
ggplot(data = calls_per_person, aes(x = log(PfMSP119))) + 
  geom_histogram(color = "white", fill = "darkgrey") + 
  ylab("Count of individuals") +
  xlab("log(MFI) for PfMSP1-19") +
  theme_bw()





or <- exp(coef(m3_rOccupation))
ci <- exp(confint(m3_rOccupation))


data <- data.frame(predictor = names(or), odds_ratio = or, lower_ci = ci[,1], upper_ci = ci[,2])

ggplot(data %>% filter(predictor != "(phi)" & predictor != "(Intercept)"), aes(x = predictor, y = odds_ratio, ymin = lower_ci, ymax = upper_ci)) +
  geom_hline(aes(yintercept = 1), size = .25, linetype = "dashed") + 
  geom_pointrange() +
  xlab("Predictor Variable") +
  ylab("Odds Ratio") +
  ggtitle("Odds Ratios for Beta Regression Model") + 
  coord_flip() +
  theme_bw()
  
