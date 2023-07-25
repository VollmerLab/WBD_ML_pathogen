
library(tidyverse)
library(tidymodels)
library(themis)
library(finetune)
library(tidytext)
library(vip)
library(kernelshap)
library(shapviz)

#### Data ####
field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                      show_col_types = FALSE) %>%
  mutate(across(domain:genus, str_replace_na)) 

#### Split data ####
field_data_wide <- select(field_data, asv_id, sample_id, log2_cpm_norm, health) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm)

taxonomy <- select(field_data, asv_id, domain:genus) %>%
  distinct

metadata <- select(field_data, sample_id, health, year, season, site) %>%
  distinct

#### Visualize Differences btw healthy/diseased and sites etc ####
library(vegan)

nmds <- metaMDS(select(field_data_wide, -health) %>%
                  column_to_rownames('sample_id'))

scores(nmds, 'sites') %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(metadata,
            by = 'sample_id') %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  
  geom_point(aes(colour = health), alpha = 1, size = 2) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black'),
        axis.title = element_text(colour = 'black', size = 14),
        panel.background = element_rect(colour = 'black'))

#### Train/Test Split ####
# not ideal because of sample size...but what can you do
# use 10 samples to build model and assess model on 5 (9 and 6 from 2/3 1/3)
# get parameters using 3 fold cv repeated 10 times
coral_split <- initial_split(field_data_wide, prop = 3/4,
                             strata = health)
coral_train <- training(coral_split)
coral_test <- testing(coral_split)

count(coral_train, health)
count(coral_test, health)

#### Random Forest Preprocessing ####
preprocess_recipie <- recipe(health ~ ., data = coral_train) %>%
  update_role(sample_id, new_role = "ID") %>%
  step_downsample(all_outcomes()) %>%
  step_nzv(all_predictors()) %>%
  step_log(all_predictors(), base = 2, offset = 0.5) %>%
  step_normalize(all_predictors()) %>%
  identity()

#### Build Random Forest Model ####
forest_model <- rand_forest() %>%
  set_engine('ranger', importance = 'permutation', probability = TRUE) %>%
  set_mode('classification') %>%
  set_args(mtry = tune(),
           trees = tune(),
           min_n = tune())

disease_wflow <- workflow() %>% 
  add_model(forest_model) %>% 
  add_recipe(preprocess_recipie)

#### Train Hyperparameters ####
coral_folds <- vfold_cv(coral_train, 
                        v = 5, repeats = 1, 
                        strata = health)

tuning_params <- parameters(list(mtry = mtry(range = c(10, 250)), 
                                 trees = trees(range = c(1e2, 1e4)),
                                 min_n = min_n(range = c(1, 4))))


# tune_aov <- tune_race_anova(disease_wflow,
#                             resamples = coral_folds,
#                             param_info = tuning_params,
#                             metrics = metric_set(accuracy),
#                             grid = 20,
#                             control = control_race(verbose = TRUE,
#                                                    verbose_elim = TRUE,
#                                                    alpha = 0.1,
#                                                    burn_in = 5))
# autoplot(tune_aov)

initial_grid_tune <- tune_grid(disease_wflow,
                               resamples = coral_folds,
                               param_info = tuning_params,
                               metrics = metric_set(accuracy),
                               grid = 10,
                               control = control_grid(verbose = TRUE))

bayes_tune <- tune_bayes(disease_wflow,
                         resamples = coral_folds,
                         param_info = tuning_params,
                         metrics = metric_set(accuracy),
                         iter = 100,
                         initial = initial_grid_tune,
                         objective = exp_improve(),
                         control = control_bayes(verbose = TRUE,
                                                 verbose_iter = TRUE,
                                                 no_improve = 10,
                                                 uncertain = 3))
autoplot(bayes_tune)

#### Finalize Hyper Parameters ####
final_rf <- finalize_workflow(disease_wflow,
                              select_best(bayes_tune, "accuracy"))


#### Assess Model Accuracy on Test Dataset ####
disease_fit <- final_rf %>%
  last_fit(coral_split)

disease_fit %>%
  collect_metrics()

disease_fit %>%
  collect_predictions() %>%
  conf_mat(health, .pred_class)

#### Fit Final Model for Variable Importance ####
fit_rf <- final_rf %>%
  fit(data = field_data_wide)

#### Get Variable Importance ####
fit_rf %>%
  extract_fit_parsnip() %>% 
  vi() %>%
  rename(asv_id = Variable) %>%
  left_join(taxonomy, by = c('asv_id')) %>%
  mutate(var = str_c(family, genus, asv_id, sep = '_')) %>%
  arrange(-Importance) %>%
  slice(1:24) %>%
  mutate(var = fct_reorder(var, Importance)) %>%
  ggplot(aes(x = Importance, y = var)) +
  geom_segment(xend = 0, aes(yend = var)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(y = NULL) +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank())

