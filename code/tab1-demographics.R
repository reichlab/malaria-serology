#Table 1 Demographic Data for Ch 3


# libraries
library(tidyverse)
library(table1)

#calls_per_person <- read.csv("datasets/training_data_marker_calls.csv")

individuals_wSerology <- read.csv("/Users/ecramer/Desktop/Lover Lab /N-Lao-Serology/datasets/Epi_seropos_merge.csv") %>%
  select(BARCODE, PvAMA1, PvMSP119, PfAMA1, PfMSP119, AGE)

raw_demo_data <- read.csv("datasets/Parademo_individual_master_v15a.csv") %>% 
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



merge_data <- individuals_wSerology %>%
  left_join(raw_demo_data) %>% 
  filter(AGE > 1)



table1(~agebands +  sex + occupation_agg + factor(fever) + ITN_coverage_bin + 
         pcr_result + pcr_species | District, data = merge_data)



variable.names(merge_data)
