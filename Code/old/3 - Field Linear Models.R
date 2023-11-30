
library(tidyverse)
library(lubridate)
library(lmerTest)
library(emmeans)

#### Data ####
field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                      show_col_types = FALSE) %>%
  mutate(across(domain:genus, str_replace_na)) 


#### Model each ASV ####
asv_models <- field_data %>%
  nest_by(asv_id, across(domain:genus))



lm(log2_cpm_norm ~ health * year * season * site, 
   data = asv_models$data[[1]]) %>%
  anova
