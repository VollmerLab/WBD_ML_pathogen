##TODO - output diagnostic plots & conf mats etc.
##TODO - decide mcc/brier/logLike/bal_accuracy for both tuning and determining equivilent models
##TODO - decide repeated vfold CV vs monte carlo CV vs Bootstrap CV
##TODO - have base model without subsampling??


library(doFuture)
library(tidyverse)
library(magrittr)
library(tidymodels)
library(tidyposterior)
library(discrim)
library(bonsai)
library(baguette)
library(plsmod)
library(rules)
library(themis)
library(finetune)

#### Data ####
field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  mutate(across(domain:genus, str_replace_na)) 

field_data %>%
  select(asv_id, domain:genus) %>%
  distinct %>%
  write_csv('../intermediate_files/taxonomy.csv.gz')

#### Split data ####
field_data_wide <- select(field_data, asv_id, sample_id, log2_cpm_norm, health) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm) %>%
  mutate(health = factor(health))

taxonomy <- select(field_data, asv_id, domain:genus) %>%
  distinct

metadata <- select(field_data, sample_id, health, year, season, site) %>%
  distinct

#### Train/Test Split ####
coral_split <- initial_split(field_data_wide, prop = 3/4,
                             strata = health)
coral_train <- training(coral_split) 
coral_test <- testing(coral_split)
write_rds(coral_split, '../intermediate_files/coral_split.rds.gz')

count(coral_train, health)
count(coral_test, health)

coral_folds <- vfold_cv(coral_train,
                        v = 10,
                        repeats = 10,
                        strata = health)


# coral_folds <- bootstraps(coral_train,
#                           times = 10,
#                           strata = health)

# coral_folds <- mc_cv(coral_train,
#                      prop = 3/4,
#                      times = 25,
#                      strata = health)



#### Preprocessing ####
preprocess_recipie <- recipe(health ~ ., data = coral_train) %>%
  update_role(sample_id, new_role = "ID") %>%
  # step_string2factor(all_outcomes()) %>%
  # step_downsample(all_outcomes(), under_ratio = tune()) %>%
  # step_downsample(all_outcomes(), under_ratio = 1) %>%
  step_nzv(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_YeoJohnson(all_predictors()) %>%
  step_normalize(all_predictors()) %>%
  step_zv(all_predictors())

correlation_recipie <- preprocess_recipie %>%
  step_corr(all_predictors(), threshold = tune()) %>%
  step_zv(all_predictors())

pca_recipie <- preprocess_recipie %>%
  step_pca(all_predictors(), num_comp = tune()) %>%
  step_zv(all_predictors())

pls_recipie <- preprocess_recipie %>%
  step_pls(all_predictors(), outcome = 'health', 
           num_comp = tune(), predictor_prop = tune()) %>%
  step_zv(all_predictors())

#### Build Models ####
default_model <- null_model() %>%
  set_engine("parsnip") %>% 
  set_mode('classification')

## "Regression" Models
ridge_model <- logistic_reg() %>%
  set_engine('glmnet') %>%
  set_mode('classification') %>%
  set_args(penalty = tune(),
           mixture = 0)

lasso_model <- logistic_reg() %>%
  set_engine('glmnet') %>%
  set_mode('classification') %>%
  set_args(penalty = tune(),
           mixture = 1)

elasticNet_model <- logistic_reg() %>%
  set_engine('glmnet') %>%
  set_mode('classification') %>%
  set_args(penalty = tune(),
           mixture = tune())

## Rules Models
c5_model <- C5_rules() %>%
  set_mode('classification') %>%
  set_args(trees = tune(),
           min_n = tune())

rule_model <- rule_fit() %>%
  set_engine('xrf',
             counts = FALSE) %>%
  set_mode('classification') %>%
  set_args(mtry = tune(),
           min_n = tune(),
           tree_depth = tune(),
           learn_rate = tune(),
           loss_reduction = tune(),
           sample_size = tune(),
           penalty = tune())

## Tree Models
forest_model <- rand_forest() %>%
  set_engine('ranger', 
             importance = 'impurity', 
             probability = TRUE) %>%
  set_mode('classification') %>%
  set_args(mtry = tune(),
           trees = tune(),
           min_n = tune()) 

boostTree_model <- boost_tree() %>%
  set_engine('xgboost',
             counts = FALSE,
             validation = 0.2) %>%
  set_mode('classification') %>%
  set_args(mtry = tune(),
           trees = tune(),
           min_n = tune(),
           tree_depth = tune(),
           learn_rate = tune(),
           loss_reduction = tune(),
           sample_size = tune(),
           stop_iter = tune())

decision_model <- decision_tree() %>%
  set_engine('rpart') %>%
  set_mode('classification') %>%
  set_args(cost_complexity = tune(),
           tree_depth = tune(),
           min_n = tune())

bart_model <- bart() %>%
  set_mode('classification') %>%
  set_args(trees = tune(),
           prior_terminal_node_coef = tune(),
           prior_terminal_node_expo = tune(),
           prior_outcome_range = tune())


## Discriminant Models
lda_model <- discrim_linear() %>%
  set_engine('MASS') %>%
  set_mode('classification') 

fda_model <- discrim_flexible() %>%
  set_engine('earth',
             nfold = 5) %>%
  set_mode('classification') %>%
  set_args(num_terms = tune(),
           prod_degree = tune(),
           prune_method = tune())

qda_model <- discrim_quad() %>%
  set_engine('sparsediscrim') %>%
  set_mode('classification') %>%
  set_args(regularization_method = tune())

rda_mode <- discrim_regularized() %>%
  set_mode('classification') %>%
  set_args(frac_common_cov = tune(),
           frac_identity = tune())

nb_model <- naive_Bayes() %>%
  set_engine('klaR') %>%
  set_mode('classification') %>%
  set_args(smoothness = tune(),
           Laplace = tune())

## SVM Models 
svmLinear_model <- svm_linear() %>%
  set_engine('kernlab') %>%
  set_mode('classification') %>%
  set_args(cost = tune())

svmPOLY_model <- svm_poly() %>%
  set_engine('kernlab') %>%
  set_mode('classification') %>%
  set_args(cost = tune(),
           degree = tune(),
           scale_factor = tune())

svmRBF_model <- svm_rbf() %>%
  set_engine('kernlab') %>%
  set_mode('classification') %>%
  set_args(cost = tune(),
           rbf_sigma = tune())

## Bagged Models 
bagMARS_model <- bag_mars() %>%
  set_engine('earth',
             nfold = 5) %>%
  set_mode('classification') %>%
  set_args(num_terms = tune(),
           prod_degree = tune(),
           prune_method = tune())

bagMLP_model <- bag_mlp() %>%
  set_mode('classification') %>%
  set_args(hidden_units = tune(),
           penalty = tune(),
           epochs = tune())

bagTree_model <- bag_tree() %>%
  set_mode('classification') %>%
  set_args(cost_complexity = tune(),
           tree_depth = tune(),
           min_n = tune(),
           class_cost = tune())

## Other Models
knn_model <- nearest_neighbor() %>%
  set_mode('classification') %>%
  set_args(neighbors = tune(),
           weight_func = tune(),
           dist_power = tune())

pls_model <- pls() %>%
  set_mode('classification') %>%
  set_args(predictor_prop = tune(),
           num_comp = tune())

mlp_model <- mlp() %>%
  set_engine('brulee') %>%
  set_mode('classification') %>%
  set_args(hidden_units = tune(),
           # penalty = tune(),
           epochs = tune(),
           dropout = tune(),
           learn_rate = tune(),
           activation = tune())

mars_model <- mars() %>%
  set_engine('earth',
             nfold = 5) %>%
  set_mode('classification') %>%
  set_args(num_terms = tune(),
           prod_degree = tune(),
           prune_method = tune())

#### Make Workflow Set ####
disease_wflow <- workflow_set(
  preproc = list(base = preprocess_recipie,
                 corr = correlation_recipie,
                 pca = pca_recipie,
                 pls = pls_recipie),
  # preproc = list(base = preprocess_recipie),
  models = list(null = default_model,
                ridge = ridge_model,
                lasso = lasso_model,
                elasticNet = elasticNet_model,
                c5 = c5_model,
                rule = rule_model,
                forest = forest_model, 
                boostTree = boostTree_model,
                decision = decision_model,
                bart = bart_model,
                lda = lda_model,
                fda = fda_model,
                qda = qda_model,
                rda = rda_mode,
                nb = nb_model,
                svmLinear = svmLinear_model,
                svmPOLY = svmPOLY_model,
                svmRBF = svmRBF_model,
                bagMARS = bagMARS_model,
                bagMLP = bagMLP_model,
                bagTree = bagTree_model,
                knn = knn_model,
                pls = pls_model,
                mlp = mlp_model,
                mars = mars_model),
  cross = TRUE) %>%
  # sample_frac(1) %>%
  # filter(wflow_id == 'base_rf') %>%
  
  #Remove combos that don't make sense
  anti_join(tibble(wflow_id = c("pca_lasso", "corr_lasso", 'pls_lasso',
                                "pca_ridge", "corr_ridge", 'pls_ridge',
                                "pca_elasticNet", "corr_elasticNet", 'pls_elasticNet')),
            by = "wflow_id") %>%
  sample_frac(1) %>% #randomize order
  
  # filter(str_detect(wflow_id, 'bagMARS', negate = TRUE)) %>%
  
  identity()

#### Update some annoying hyperparam settings ####
for(model_id in str_subset(disease_wflow$wflow_id, 'forest$|boostTree$|fda$|qda$|mars|MARS')){ #
  
  if(str_detect(model_id, 'fda')){
    update_pars <- extract_workflow(disease_wflow, model_id) %>%
      extract_parameter_set_dials() %>%
      update(num_terms = num_terms(c(1, ncol(coral_train) - 1)),
             prune_method = prune_method(c('backward', 'none', 'forward', 'seqrep')))
  } else if(str_detect(model_id, 'forest$')) {
    update_pars <- extract_workflow(disease_wflow, model_id) %>%
      extract_parameter_set_dials() %>%
      update(mtry = mtry(c(1, ncol(coral_train) - 1)))
  } else if(str_detect(model_id, 'boostTree$')) {
    update_pars <- extract_workflow(disease_wflow, model_id) %>%
      extract_parameter_set_dials() %>%
      update(mtry = mtry_prop())
  } else if(str_detect(model_id, 'qda$')) {
    update_pars <- extract_workflow(disease_wflow, model_id) %>%
      extract_parameter_set_dials() %>%
      update(regularization_method = regularization_method(c('diagonal',
                                                             'shrink_cov',
                                                             'shrink_mean')))
  } else if(str_detect(model_id, 'mars|MARS')) {
    update_pars <- extract_workflow(disease_wflow, model_id) %>%
      extract_parameter_set_dials() %>%
      update(prune_method = prune_method(c('backward', 'none', 'forward', 'seqrep')))
  } 
  
  disease_wflow <- disease_wflow %>%
    option_add(param_info = update_pars,
               id = model_id)
}

# disease_wflow <- disease_wflow %>%
#   # filter(str_detect(wflow_id, 'base')) %>%
#   filter(str_detect(wflow_id, 'lda'))


#### Train Hyperparameters ####
registerDoFuture()
plan(cluster, workers = parallel::makeCluster(parallel::detectCores(logical = FALSE)))

if(file.exists('../intermediate_files/model_tune_run.rds.gz')){
  disease_model_tune <- read_rds('../intermediate_files/model_tune_run.rds.gz')
} else {
  disease_model_tune <- workflow_map(disease_wflow,
                                     "tune_bayes", 
                                     resamples = coral_folds, 
                                     initial = 25,
                                     iter = 100,
                                     metrics = metric_set(mn_log_loss, brier_class,
                                                          bal_accuracy, accuracy,
                                                          j_index, roc_auc,
                                                          mcc, kap,
                                                          sens, spec),
                                     verbose = TRUE,
                                     control = control_bayes(verbose = FALSE,
                                                             allow_par = TRUE,
                                                             no_improve = 25L,
                                                             uncertain = 10L))
  write_rds(disease_model_tune, '../intermediate_files/model_tune_run.rds.gz')
}

autoplot(disease_model_tune, metric = "mn_log_loss")
autoplot(disease_model_tune, select_best = TRUE, 
         metric = "mn_log_loss", rank_metric = "mn_log_loss")
autoplot(disease_model_tune, id = 'base_lda', metric = "mn_log_loss")

rank_results(disease_model_tune, rank_metric = "mn_log_loss", select_best = TRUE) %>%
  filter(.metric == 'mn_log_loss') 
 
collect_metrics(disease_model_tune, summarize = FALSE)  %>%
  pivot_wider(names_from = '.metric',
              values_from = '.estimate') %>%
  write_csv('../intermediate_files/model_tuning_metrics.csv.gz')


#### Model Model Equivilance ####
get_best_fit_stats <- function(wflows, metric){
  top_fits <- rank_results(wflows, rank_metric = metric, select_best = TRUE) %>%
    filter(.metric == metric)  %>%
    dplyr::select(wflow_id, .config)
  
  collect_metrics(wflows, summarize = FALSE) %>%
    inner_join(top_fits,
               by = join_by(wflow_id, .config)) %>%
    filter(.metric == metric) %>%
    dplyr::select(wflow_id, starts_with('id'), .estimate) %>%
    pivot_wider(names_from = wflow_id,
                values_from = .estimate)
}

model_class_tweak <- function(model, metric, direction = 'minimize'){
  class(model) <- c("perf_mod_workflow_set", "perf_mod")
  model$metric$name <- metric
  model$metric$direction <- direction
  model
}

model_equivilance_model <- get_best_fit_stats(disease_model_tune, 'bal_accuracy') %>%
  
  perf_mod(#formula = ,
           #transform = logit_trans,
           # family = 'binomial',
           hetero_var = FALSE,
           
           chains = 4,
           cores = 4,
           iter = 2000,
           warmup = 1000,
           refresh = 100) %>%
  model_class_tweak(metric = 'bal_accuracy', direction = 'maximize')

autoplot(model_equivilance_model, type = "intervals", prob = 0.9)
autoplot(model_equivilance_model, size = 0.05, type = "ROPE") #1% difference in accuracy

#### Pick Equivilant Models ####
find_equivilent_models <- function(wflowset, size, cutoff = NA){
  if(is.na(cutoff)){
    cutoff <- 1
  }
  min_max <- wflowset$metric$direction
  if(min_max == 'minimize'){
    posterior <- tidy(wflowset) %>%
      summary %>%
      arrange(mean)
  } else {
    posterior <- tidy(wflowset) %>%
      summary %>%
      arrange(-mean)
  }
  
  
  equivilance <- contrast_models(wflowset, 
                                 list_1 = rep(posterior$model[[1]], nrow(posterior)),
                                 list_2 = posterior$model)  %>%
    summary(size = size) %>%
    mutate(contrast = str_remove(contrast, str_c(posterior$model[[1]], ' vs '))) %>%
    dplyr::select(contrast, pract_equiv)
  
  left_join(posterior,
            equivilance,
            by = c('model' = 'contrast')) %>%
    filter(pract_equiv >= 1 - cutoff)
}

top_wflows <- find_equivilent_models(model_equivilance_model, 0.05, 0.3)

top_disease_models <- top_wflows %>%
  rowwise %>%
  mutate(wflowset = list(extract_workflow_set_result(disease_model_tune, model)),
         # race_plot = list(plot_race(wflowset)),
         best_params = list(select_best(wflowset, 'bal_accuracy')),
         final_wflow = list(extract_workflow(disease_model_tune, model) %>%
                              finalize_workflow(best_params)),
         
         final_wflow = list(last_fit(final_wflow, coral_split, 
                                    metrics = metric_set(bal_accuracy, accuracy, 
                                                         brier_class, mn_log_loss,
                                                         f_meas, j_index,
                                                         mcc, kap,
                                                         precision, recall,
                                                         sens, spec))),
         
         collect_metrics(final_wflow) %>%
           dplyr::select(.metric, .estimate) %>%
           pivot_wider(names_from = '.metric',
                       values_from = '.estimate'),
         
         conf_mat = list(collect_predictions(final_wflow) %>%
                           conf_mat(health, .pred_class))) %>%
  dplyr::select(-wflowset)

#### Assess Model Accuracy on Test Dataset ####
top_disease_models %>%
  dplyr::select(model, bal_accuracy:mn_log_loss) %>%
  arrange(-bal_accuracy)

top_disease_models  %>%
  arrange(-bal_accuracy) %>%
  dplyr::select(model, conf_mat) %$%
  set_names(conf_mat, model)

#### Write outputs ####
dplyr::select(top_disease_models, -mean:-upper, 
              -contains('race_plot'), -conf_mat) %>%
  ungroup %>%
  write_rds('../intermediate_files/top_disease_models.rds.gz')
