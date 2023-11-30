library(tidyverse)
library(tidyposterior)
library(broom)


all_models <- read_rds('../../Results/model_tuning_metrics.rds.gz') %>%
  group_by(wflow_id) %>%
  filter(.config == .config[train_test == 'test'],
         train_test != 'tank') %>%
  group_by(wflow_id, train_test) %>%
  summarise(across(c(accuracy, roc_auc, brier_class, overall), list(mean = mean, se = ~sd(.) / sqrt(n()))),
            .groups = 'drop') %>%
  left_join(read_csv('../../Results/equivilant_top_models.csv.gz', show_col_types = FALSE) %>%
              select(model, pract_equiv),
            by = c('wflow_id' = 'model')) %>%
  mutate(wflow_id = str_remove(wflow_id, 'base_') %>%
           str_replace_all(c('forest' = 'Random Forest',
                             'knn' = 'KNN',
                             'lasso' = 'LASSO',
                             'mlp' = 'MLP',
                             'null' = 'Null',
                             'pls' = 'PLS',
                             'svmLinear' = 'SVM')))



out_table <- all_models %>%
  pivot_longer(cols = -c(wflow_id, train_test, pract_equiv),
               names_to = c('metric', '.value'),
               names_pattern = '(.*)_(mean|se)') %>%
  mutate(across(c(mean, se), ~case_when(metric %in% c('brier_class', 'roc_auc') ~ round(., 3) %>% 
                                          as.character,
                                        TRUE ~ scales::percent(., accuracy = 0.1))),
         se = str_replace_na(se, ''),
         out = str_c(mean, se, sep = ' +- ') %>%
           str_remove(' \\+- $')) %>%
  select(-mean, -se) %>%
  pivot_wider(names_from = train_test,
              values_from = out) %>%
  mutate(out = str_c(train, ' (', test, ')')) %>%
  select(-test, -train) %>%
  pivot_wider(names_from = metric, 
              values_from = 'out') %>%
  select(wflow_id, overall, everything()) %>%
  mutate('Practical Equivilance' = scales::percent(pract_equiv, accuracy = 0.1), .keep = 'unused') %>%
  left_join(filter(all_models, train_test == 'train') %>%
              select(wflow_id, overall_mean),
            by = 'wflow_id') %>%
  arrange(-overall_mean) %>%
  select(-overall_mean) %>%
  rename('Algorithm' = 'wflow_id') %>%
  rename_with(~str_replace(., '_', ' ') %>%
                str_to_title()) %>%
  rename('ROC/AUC' = 'Roc Auc')
write_csv(out_table, '../../Results/Table1_modelMetrics.csv')





model_equivilance_model <- read_rds('../../intermediate_files/model_equivilance.rds.gz')

tidy(model_equivilance_model) %>%
  summary(model, prob = 0.95) %>%
  arrange(-mean) %>%
  mutate(train_test = 'train',
         model = fct_reorder(model, -mean))

