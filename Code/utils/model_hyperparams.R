
library(tidyverse)
library(tidymodels)
library(tidyposterior)
library(broom)

the_models <- read_csv('../../Results/equivilant_top_models.csv.gz', show_col_types = FALSE) %>%
  select(model, pract_equiv) %>%
  filter(pract_equiv > 0.8) %>%
  mutate(model_file = str_c('../../Results/model_tuning/models/', model, '.rds.gz')) %>%
  rowwise(model) %>%
  summarise(model_fit = list(read_rds(model_file)),
            .groups = 'drop')


the_models$model_fit[[1]]
the_models$model_fit[[2]]
the_models$model_fit[[3]]
the_models$model_fit[[4]]
