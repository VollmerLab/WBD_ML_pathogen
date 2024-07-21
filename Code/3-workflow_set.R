if(!interactive()){
  rerun_process <- TRUE
  reduce_modelset <- FALSE
  long_tune <- TRUE
} else {
  rerun_process <- FALSE
  reduce_modelset <- FALSE
  long_tune <- FALSE
}

tune_variant <- 'bayes' #can be bayes or race
maxZ <- 2

# To run on Discovery 
# sbatch --output=/scratch/j.selwyn/Panama_Tank_Field/slurm/workflowset_%j.out \
#   --mem=200G \
#   /work/vollmer/software/jds_scripts/runRscript.slurm \
#   /scratch/j.selwyn/Panama_Tank_Field/Code/3-workflow_set.R

library(doFuture)
library(batchtools)
library(tidyverse)
library(magrittr)
library(patchwork)
library(tidymodels)
library(tidyposterior)
library(discrim)
library(bonsai)
library(baguette)
library(plsmod)
library(rules)
library(themis)
library(finetune)
library(desirability2)



if(Sys.info()['sysname'] != 'Windows'){
  setwd('/scratch/j.selwyn/Panama_Tank_Field/Code')
} else {
  setwd("~/Google Drive/Research/Vollmer Lab PostDoc/Panama_Tank_Field/Code")
}

#### Data ####
field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  mutate(across(domain:genus, str_replace_na)) 

tank_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                      show_col_types = FALSE)

if(!file.exists('../intermediate_files/taxonomy.csv.gz')){
  field_data %>%
    select(asv_id, domain:genus) %>%
    distinct %>%
    write_csv('../intermediate_files/taxonomy.csv.gz')
}

#### Split data ####
field_data_wide <- select(field_data, asv_id, sample_id, 
                          log2_cpm_norm, health) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm) %>%
  mutate(health = factor(health))

taxonomy <- select(field_data, asv_id, domain:genus) %>%
  distinct

metadata <- select(field_data, sample_id, health, year, season, site) %>%
  distinct

tank_data_wide <- tank_data %>%
  filter(geno != 'GE', time == '8_exp',
         anti == 'N', exposure == 'D') %>%
  select(asv_id, sample_id, health, log2_cpm_norm) %>%
  mutate(health = factor(health, levels = c('D', 'H'))) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm)

#### Train/Test Split ####
if(file.exists('../intermediate_files/coral_split.rds.gz') & !rerun_process){
  coral_split <- read_rds('../intermediate_files/coral_split.rds.gz')
} else {
  coral_split <- initial_split(field_data_wide, prop = 3/4,
                               strata = health)
  write_rds(coral_split, '../intermediate_files/coral_split.rds.gz')
}

coral_train <- training(coral_split) 
coral_test <- testing(coral_split)

count(coral_train, health)
count(coral_test, health)

if(file.exists('../intermediate_files/coral_folds.rds.gz') & !rerun_process){
  coral_folds <- read_rds('../intermediate_files/coral_folds.rds.gz')
} else {
  coral_folds <- vfold_cv(coral_train,
                          v = 10,
                          repeats = 10,
                          strata = health)
  
  write_rds(coral_folds, '../intermediate_files/coral_folds.rds.gz')
}


#### Preprocessing ####
preprocess_recipie <- recipe(health ~ ., data = coral_train) %>%
  update_role(sample_id, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_YeoJohnson(all_predictors()) %>%
  step_normalize(all_predictors()) %>%
  step_zv(all_predictors())

#### Build Models ####
default_model <- null_model() %>%
  set_engine("parsnip") %>% 
  set_mode('classification')

## "Regression" Models
lasso_model <- logistic_reg() %>%
  set_engine('glmnet') %>%
  set_mode('classification') %>%
  set_args(penalty = tune(),
           mixture = 1)

## Tree Models
forest_model <- rand_forest() %>%
  set_engine('ranger', 
             importance = 'impurity', 
             probability = TRUE) %>%
  set_mode('classification') %>%
  set_args(mtry = tune(),
           trees = tune(),
           min_n = tune()) 

## SVM Models 
svmLinear_model <- svm_linear() %>%
  set_engine('kernlab') %>%
  set_mode('classification') %>%
  set_args(cost = tune())

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
           epochs = tune(),
           dropout = tune(),
           learn_rate = tune(),
           activation = tune())

#### Make Workflow Set ####
disease_wflow <- workflow_set(
  preproc = list(base = preprocess_recipie),
  models = list(null = default_model,
                lasso = lasso_model,
                forest = forest_model, 
                svmLinear = svmLinear_model,
                knn = knn_model,
                pls = pls_model,
                mlp = mlp_model),
  cross = TRUE) %>%
  sample_frac(1) %>% #randomize order
  identity()

#### Update some annoying hyperparam settings ####
for(model_id in str_subset(disease_wflow$wflow_id, 'forest$')){ #

  update_pars <- extract_workflow(disease_wflow, model_id) %>%
    extract_parameter_set_dials() %>%
    update(mtry = mtry(c(1, floor((2/3) * (ncol(coral_train) - 1)))))

  disease_wflow <- disease_wflow %>%
    option_add(param_info = update_pars,
               id = model_id)
}

#### Train Hyperparameters ####
if(reduce_modelset){
  disease_wflow <- sample_n(disease_wflow, 2)
}

if(tune_variant == 'bayes'){
  workflow_tuning <- function(wflowset){
    workflow_map(wflowset,
                 "tune_bayes", 
                 resamples = coral_folds, 
                 initial = if_else(long_tune, 50, 5),
                 iter = if_else(long_tune, 200, 10),
                 metrics = metric_set(brier_class, mn_log_loss,
                                      bal_accuracy, accuracy,
                                      j_index, roc_auc,
                                      mcc, kap,
                                      sens, spec),
                 verbose = TRUE,
                 control = control_bayes(verbose = if_else(Sys.info()['sysname'] != 'Windows', FALSE, TRUE),
                                         allow_par = if_else(Sys.info()['sysname'] != 'Windows', TRUE, FALSE),
                                         no_improve = if_else(long_tune, 25L, 2L),
                                         uncertain = if_else(long_tune, 10L, 5L)))
  }
} else if(tune_variant == 'race'){
  workflow_tuning <- function(wflowset){
    workflow_map(wflowset,
                 "tune_race_anova", 
                 resamples = coral_folds, 
                 grid = if_else(long_tune, 200, 10),
                 metrics = metric_set(brier_class, mn_log_loss,
                                      bal_accuracy, accuracy,
                                      j_index, roc_auc,
                                      mcc, kap,
                                      sens, spec),
                 verbose = TRUE,
                 control = control_race(verbose = if_else(Sys.info()['sysname'] != 'Windows', FALSE, TRUE),
                                        allow_par = if_else(Sys.info()['sysname'] != 'Windows', TRUE, FALSE),
                                        burn_in = if_else(long_tune, 5, 2),
                                        num_ties = if_else(long_tune, 10, 10),
                                        alpha = 0.05))
  }
}


workflow_post_process <- function(wflows, rank_metric, train_test_split, tankDat, outdir){
  top_fits <- rank_results(wflows, rank_metric = rank_metric, select_best = TRUE) %>%
    dplyr::select(wflow_id, .config) %>%
    distinct
  
  wfset_results <- extract_workflow_set_result(wflows, top_fits$wflow_id) 
  
  train_metrics <- collect_metrics(wfset_results, summarize = FALSE) %>%
    pivot_wider(names_from = '.metric', 
                values_from = '.estimate') %>%
    dplyr::select(-.estimator)
  
  best_params <- select_best(wfset_results, rank_metric)
  
  if(str_detect(top_fits$wflow_id, 'null')){
    final_wflow_fit <- extract_workflow(wflows, top_fits$wflow_id) 
  } else {
    final_wflow_fit <- extract_workflow(wflows, top_fits$wflow_id) %>%
      finalize_workflow(best_params)
  }
  
  write_rds(final_wflow_fit, paste0(outdir, '/models/', top_fits$wflow_id, '.rds.gz'))
  
  final_wflow <- last_fit(final_wflow_fit, train_test_split, 
                          metrics = metric_set(bal_accuracy, accuracy, 
                                               brier_class, mn_log_loss,
                                               f_meas, j_index,
                                               mcc, kap, roc_auc,
                                               precision, recall,
                                               sens, spec))
  
  test_metrics <- collect_metrics(final_wflow) %>%
    dplyr::select(.metric, .estimate) %>%
    pivot_wider(names_from = '.metric',
                values_from = '.estimate') %>%
    bind_cols(best_params, .)
  
  #add tank tankDat fit & metrics
  all_metrics <- metric_set(bal_accuracy, accuracy, 
                            brier_class, mn_log_loss,
                            f_meas, j_index,
                            mcc, kap, roc_auc,
                            precision, recall,
                            sens, spec)
  
  tank_metrics <- bind_cols(dplyr::select(tankDat, sample_id, health), 
                            predict(final_wflow$.workflow[[1]], new_data = tankDat, type = 'class'),
                            predict(final_wflow$.workflow[[1]], new_data = tankDat, type = 'prob')) %>%
    all_metrics(truth = health, estimate = .pred_class, .pred_D) %>%
    dplyr::select(.metric, .estimate) %>%
    pivot_wider(names_from = '.metric',
                values_from = '.estimate') %>%
    bind_cols(best_params, .)
  
  
  out_metrics <- bind_rows(tank = tank_metrics,
                           test = test_metrics,
                           train = train_metrics,
                           .id = 'train_test') %>%
    relocate(.config, starts_with('id'), contains('.iter'), 
             .after = train_test) 
  
  write_csv(out_metrics, paste0(outdir, '/metrics/', top_fits$wflow_id, '_metrics.csv.gz'))
  
  predictions_out <- bind_rows(
    train = bind_cols(dplyr::select(analysis(train_test_split), sample_id, health), 
                      predict(final_wflow$.workflow[[1]], new_data = analysis(train_test_split), type = 'class'),
                      predict(final_wflow$.workflow[[1]], new_data = analysis(train_test_split), type = 'prob')),
    
    test = bind_cols(dplyr::select(assessment(train_test_split), sample_id, health), 
                     predict(final_wflow$.workflow[[1]], new_data = assessment(train_test_split), type = 'class'),
                     predict(final_wflow$.workflow[[1]], new_data = assessment(train_test_split), type = 'prob')),
    
    tank = bind_cols(dplyr::select(tankDat, sample_id, health), 
                     predict(final_wflow$.workflow[[1]], new_data = tankDat, type = 'class'),
                     predict(final_wflow$.workflow[[1]], new_data = tankDat, type = 'prob')),
    .id = 'train_test'
  ) 
  
  write_csv(predictions_out, paste0(outdir, '/metrics/', top_fits$wflow_id, '_prediction.csv.gz'))
  
  nest(out_metrics, params = colnames(dplyr::select(best_params, -.config))) %>%
    mutate(wflow_id = top_fits$wflow_id, .before = 'train_test')
}

slurm_wflow_tune <- function(wflowset){
  if(Sys.info()['sysname'] != 'Windows'){
    registerDoFuture()
    plan(cluster, workers = 20)
  }
  
  trained_wflow <- workflow_tuning(wflowset)
  
  all_error <- map_lgl(trained_wflow$result, ~'try-error' %in% class(.)) %>%
    which
  
  if(length(all_error) > 0){
    none_worked <- logical()
  } else {
    none_worked <- map_lgl(trained_wflow$result, ~all(map_lgl(.x$.metrics, is.null))) %>%
      which
  }
  
  
  if(!(length(none_worked) > 0 | length(all_error) > 0)){
    ## Output plot
    plot_location <- paste0(file_out, '/plots')
    top_config <- rank_results(trained_wflow, rank_metric = 'brier_class', select_best = TRUE) %>%
      pull(.config) %>%
      unique
    
    # plot_race(trained_wflow$result[[1]])
    
    metric_table <- collect_metrics(trained_wflow, summarize = TRUE) 
    has_iterations <- '.iter' %in% colnames(metric_table)
    
    the_plot <- metric_table %>%
      {if(has_iterations) . else mutate(., .iter = 1)} %>%
      
      ggplot(aes(x = .iter, y = mean, ymin = mean - 1.96 * std_err, ymax = mean + 1.96 * std_err,
                 colour = .config == top_config)) +
      geom_pointrange(position = position_dodge(width = 0.5)) +
      facet_wrap(~.metric, scales = 'free_y') +
      scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
      theme_classic() +
      theme(legend.position = 'none')
    ggsave(paste0(plot_location, '/', trained_wflow$wflow_id[[1]], '.png'),
           plot = the_plot,
           height = 15, width = 15)
    
    
    ## Get model metrics
    ## Fit Final Model
    ## Get test metrics
    out <- workflow_post_process(trained_wflow, rank_metric = 'brier_class', 
                                 train_test_split = coral_split,
                                 tankDat = tank_data,
                                 outdir = file_out)
    
  } else {
    out <- tibble(wflow_id = trained_wflow$wflow_id[[1]])
  }
  
  out
}

remove_errors <- function(wfset){
  
  all_error <- map_lgl(wfset$result, ~'try-error' %in% class(.)) %>%
    which
  
  if(length(all_error) > 0){
    message('Models failed with all errors: ', str_c(wfset$wflow_id[all_error], collapse = '; '))
    wfset <- slice(wfset, -1 * all_error)
  } 
  
  
  none_worked <- map_lgl(wfset$result, ~all(map_lgl(.x$.metrics, is.null))) %>%
    which
  
  if(length(none_worked) > 0){
    message('Models failed because none worked: ', str_c(wfset$wflow_id[none_worked], collapse = '; '))
    wfset <- slice(wfset, -1 * none_worked)
  } 
  
  wfset
}

calc_desire <- function(data){
  data %>%
    mutate(acc_d = d_max(accuracy, use_data = TRUE),
           # loss_d = d_min(mn_log_loss, use_data = TRUE),
           brier_d = d_min(brier_class, use_data = TRUE),
           rocAUC_d = d_max(roc_auc, use_data = TRUE),
           overall = d_overall(across(ends_with("_d")))) %>%
    select(-ends_with("_d"))
}

if(file.exists('../Results/model_tuning_metrics.rds.gz') & !rerun_process){
  disease_model_tune <- read_rds('../Results/model_tuning_metrics.rds.gz')
} else {
  if(Sys.info()['sysname'] != 'Windows'){
    reg <- makeRegistry(file.dir = paste0(getwd(), '/../batch_files'), 
                        packages = c('batchtools', 'doFuture', 'tidyverse', 'tidymodels',
                                     'discrim', 'bonsai', 'baguette', 
                                     'plsmod', 'rules', 'themis', 'finetune'))
    
    reg$cluster.functions <- makeClusterFunctionsSlurm(template = "~/slurm_template.tmpl",
                                                       array.jobs = TRUE,
                                                       nodename = "localhost",
                                                       scheduler.latency = 1,
                                                       fs.latency = 65)
    
    successful_jobs <- character()
    z <- 0
    newly_out <- vector('list', maxZ)
    walk(paste0(getwd(), '/../Results/model_tuning/', c('plots', 'metrics', 'models')),
         ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))
    
    
    while(length(successful_jobs) < nrow(disease_wflow) & z < maxZ){
      z <- z + 1
      
      incomplete_wflows <- filter(disease_wflow, !wflow_id %in% successful_jobs)
      all_wflows <- split(incomplete_wflows, 1:nrow(incomplete_wflows))
      
      
      if(z > 1){
        message('Running WorkFlow IDs: ', str_c(incomplete_wflows$wflow_id, collapse = '; '))
      }
      
      
      batchMap(fun = slurm_wflow_tune, wflowset = all_wflows)
      batchExport(list(coral_folds = coral_folds,
                       coral_split = coral_split,
                       tank_data = tank_data_wide,
                       workflow_tuning = workflow_tuning,
                       workflow_post_process = workflow_post_process,
                       file_out = paste0(getwd(), '/../Results/model_tuning'),
                       long_tune = long_tune))
      submitJobs(resources = list(max.concurrent.jobs = 30))
      waitForJobs()
      
      newly_finished <- dplyr::as_tibble(getJobTable()) %>%
        dplyr::filter(is.na(error) & !is.na(done)) %>%
        dplyr::select(job.pars) %>%
        unnest(job.pars) %>%
        unnest(job.pars)
      
      if(nrow(newly_finished) > 0){
        newly_finished <- pull(newly_finished, wflow_id)
        newly_out[[z]] <- purrr::map_dfr(which(incomplete_wflows$wflow_id %in% newly_finished), loadResult)
        successful_jobs <- c(successful_jobs, newly_finished)
      } 
      
      clearRegistry()
    }
    
    disease_model_tune <- bind_rows(newly_out) %>%
      full_join(dplyr::select(disease_wflow, wflow_id),
                by = 'wflow_id') %>%
      calc_desire()
    
    #write csv/rds of all metrics
    write_rds(disease_model_tune, '../Results/model_tuning_metrics.rds.gz')
    
    #plot top model metrics
    model_tune_dat <- disease_model_tune %>%
      
      group_by(wflow_id) %>%
      #get rid of failures
      filter(n() > 1) %>%
      filter(.config == .config[train_test == 'test']) %>%
      filter(!train_test %in% c('test', 'tank')) %>%
      ungroup %>%
      dplyr::select(-params) %>%
      mutate(wflow_id = fct_reorder(wflow_id, overall, .desc = TRUE,
                                    .fun = mean, na.rm = TRUE)) 
    
    model_tune_plot <- disease_model_tune %>%
      select(-params) %>%
      group_by(wflow_id) %>%
      #get rid of failures
      filter(n() > 1) %>%
      filter(.config == .config[train_test == 'test']) %>%
      ungroup %>%
      pivot_longer(cols = c(bal_accuracy:overall),
                   names_to = '.metric') %>%
      filter(!is.na(value)) %>%
      group_by(train_test, wflow_id, .metric) %>%
      summarise(mean = mean(value),
                se = sd(value) / sqrt(n()),
                .groups = 'drop') %>%
      mutate(wflow_id = factor(wflow_id, levels = levels(model_tune_dat$wflow_id))) %>%
      
      ggplot(aes(x = as.integer(wflow_id), y = mean, ymin = mean - se, ymax = mean + se,
                 colour = wflow_id, shape = train_test)) +
      geom_pointrange() +
      facet_wrap(~ .metric, scales = 'free_y') +
      theme_classic()
    
    ggsave('../Results/model_tuning_results.png', plot = model_tune_plot, 
           height = 21, width = 21)
    
  } else {
    incomplete_wflows <- filter(disease_wflow, wflow_id %in% c('base_null', 'base_lasso'))
    all_wflows <- split(incomplete_wflows, 1:nrow(incomplete_wflows))
    
    wflowset <- all_wflows[[1]]
  }
}

#### Model Model Equivilance ####
if(file.exists('../intermediate_files/model_equivilance.rds.gz') & !rerun_process){
  model_equivilance_model <- read_rds('../intermediate_files/model_equivilance.rds.gz')
} else {
  
  top_stats_all_models <- disease_model_tune %>%
    group_by(wflow_id) %>%
    #get rid of failures
    filter(n() > 1) %>%
    filter(.config == .config[train_test == 'test']) %>%
    filter(!train_test %in% c('test', 'tank')) %>%
    ungroup %>%
    dplyr::select(-params) %>%
    
    dplyr::select(wflow_id, starts_with('id'), overall) %>%
    pivot_wider(names_from = 'wflow_id',
                values_from = 'overall') %>%
    identity()
  
  model_class_tweak <- function(model, metric, direction = 'minimize'){
    class(model) <- c("perf_mod_workflow_set", "perf_mod")
    model$metric$name <- metric
    model$metric$direction <- direction
    model
  }
  
  model_equivilance_model <- top_stats_all_models %>%
    perf_mod(hetero_var = FALSE,
             
             chains = 20,
             cores = 20,
             iter = 2000,
             warmup = 1000,
             refresh = 10,
             adapt_delta = 0.99, 
             control = list(max_treedepth = 20)) %>%
    model_class_tweak(metric = 'overall', direction = 'maximize')
  
  write_rds(model_equivilance_model, '../intermediate_files/model_equivilance.rds.gz')
}

perf_sum_out <- tidy(model_equivilance_model) %>%
  summary(model, prob = 0.9) %>%
  arrange(-mean) %>%
  mutate(train_test = 'train',
         model = fct_reorder(model, -mean))

model_accuracies <- filter(disease_model_tune, train_test %in% c('test', 'tank')) %>%
  select(train_test, wflow_id, overall) %>%
  rename(mean = overall,
         model = wflow_id) %>%
  bind_rows(perf_sum_out) %>%
  left_join(perf_sum_out %>%
              mutate(model_rank = as.integer(model)) %>%
              select(model, model_rank),
            by = 'model') %>%
  filter(train_test != 'tank') %>%
  mutate(train_test = str_to_sentence(train_test)) %>%
  ggplot(aes(x = model_rank, y = mean, ymin = lower, ymax = upper,
             colour = model, shape = train_test)) +
  geom_errorbar(width = 0.1, show.legend = FALSE) +
  geom_point() +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format()) +
  guides(colour = 'none') +
  labs(x = 'Model Rank',
       y = 'Model Quality',
       shape = NULL) +
  theme_classic()


model_ropes <- autoplot(model_equivilance_model, size = 0.01, type = "ROPE") + #5% difference in joint thing
  geom_hline(yintercept = 1 - 0.2, linetype = 'dashed') +
  scale_colour_discrete('Algorithm',
                        labels = ~str_remove(., 'base_') %>% 
                          str_replace_all(c('forest' = 'Random\nForest',
                                            'knn' = 'KNN',
                                            'lasso' = 'LASSO',
                                            'mlp' = 'MLP',
                                            'null' = 'Null',
                                            'pls' = 'PLS',
                                            'svmLinear' = 'SVM')))

combined_model_plot <- model_accuracies + model_ropes + plot_layout(nrow = 1) & theme_classic()
ggsave('../Results/post_tuning_model_quality.png', plot = combined_model_plot, height = 7, width = 14)

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

top_wflows <- find_equivilent_models(model_equivilance_model, 0.01)
write_csv(top_wflows, '../Results/equivilant_top_models.csv.gz')


#### Remove Directory ####
if(Sys.info()['sysname'] != 'Windows'){
  unlink(paste0(getwd(), '/../batch_files'), recursive = TRUE)
}
