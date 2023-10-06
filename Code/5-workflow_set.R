##TODO - output diagnostic plots & conf mats etc.
##TODO - decide mcc/brier/logLike/bal_accuracy for both tuning and determining equivilent models

# if(Sys.info()['sysname'] != 'Windows'){
#   .libPaths("~/R/x86_64-pc-linux-gnu-library/4.2")
# }

##TODO - figure out wtf happened with the corr_rda 

if(!interactive()){
  rerun_process <- TRUE
  reduce_modelset <- FALSE
} else {
  rerun_process <- FALSE
  reduce_modelset <- TRUE
}

maxZ <- 2

# To run on Discovery 
# sbatch --output=/scratch/j.selwyn/Panama_Tank_Field/slurm/workflowset_%j.out \
#   --mem=200G \
#   /work/vollmer/software/jds_scripts/runRscript.slurm \
#   /scratch/j.selwyn/Panama_Tank_Field/Code/5-workflow_set.R

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



if(Sys.info()['sysname'] != 'Windows'){
  setwd('/scratch/j.selwyn/Panama_Tank_Field/Code')
}


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
  
  # coral_folds <- bootstraps(coral_train,
  #                           times = 10,
  #                           strata = health)

  # coral_folds <- mc_cv(coral_train,
  #                      prop = 3/4,
  #                      times = 25,
  #                      strata = health)
  
  write_rds(coral_folds, '../intermediate_files/coral_folds.rds.gz')
}


#### Preprocessing ####
preprocess_recipie <- recipe(health ~ ., data = coral_train) %>%
  update_role(sample_id, new_role = "ID") %>%
  # step_pct_missing(pct = tune())
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
  
  #Remove combos that don't make sense
  anti_join(tibble(wflow_id = c("pca_lasso", "corr_lasso", 'pls_lasso',
                                "pca_ridge", "corr_ridge", 'pls_ridge',
                                "pca_elasticNet", "corr_elasticNet", 'pls_elasticNet',
                                'base_mars', 'base_bagMARS', 'pls_pls')),
            by = "wflow_id") %>%
  sample_frac(1) %>% #randomize order
  
  #remove models with excessive errors fitting sample models
  filter(str_detect(wflow_id, 'rule', negate = TRUE)) %>% 
  
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

#### Train Hyperparameters ####
if(reduce_modelset){
  disease_wflow <- disease_wflow %>% #pull(wflow_id) %>% sort
    filter(!wflow_id %in% c('base_knn', 'base_bagMLP', 'base_svmLinear', 'base_svmPOLY', 'base_svmRBF',
                           'base_nb', 'base_lda', 'base_qda', 'base_bagMLP'),
           !str_detect(wflow_id, 'corr|pca'))
}


workflow_tuning <- function(wflowset){
  workflow_map(wflowset,
               "tune_bayes", 
               resamples = coral_folds, 
               initial = 50,
               iter = 200,
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
}


if(file.exists('../intermediate_files/model_tune_run.rds.gz') & !rerun_process){
  disease_model_tune <- read_rds('../intermediate_files/model_tune_run.rds.gz')
} else {
  
  if(Sys.info()['sysname'] != 'Windows'){
    
    workflow_post_process <- function(wflows, rank_metric, train_test_split, outdir){
      top_fits <- rank_results(wflows, rank_metric = rank_metric, select_best = TRUE) %>%
        dplyr::select(wflow_id, .config) %>%
        distinct
      
      wfset_results <- extract_workflow_set_result(wflows, top_fits$wflow_id) 
      
      train_metrics <- collect_metrics(wfset_results, summarize = FALSE) %>%
        pivot_wider(names_from = '.metric', 
                    values_from = '.estimate') %>%
        dplyr::select(-.estimator)
      
      best_params <- select_best(wfset_results, rank_metric)
      
      final_wflow_fit <- extract_workflow(wflows, top_fits$wflow_id) %>%
        finalize_workflow(best_params)
      #This breaks with null models
      
      write_rds(final_wflow_fit, paste0(outdir, '/models/', top_fits$wflow_id, '.rds.gz'))
      
      final_wflow <- last_fit(final_wflow_fit, train_test_split, 
                              metrics = metric_set(bal_accuracy, accuracy, 
                                                   brier_class, mn_log_loss,
                                                   f_meas, j_index,
                                                   mcc, kap,
                                                   precision, recall,
                                                   sens, spec))
      
      test_metrics <- collect_metrics(final_wflow) %>%
        dplyr::select(.metric, .estimate) %>%
        pivot_wider(names_from = '.metric',
                    values_from = '.estimate') %>%
        bind_cols(best_params, .)
      
      
      out_metrics <- bind_rows(test = test_metrics,
                               train = train_metrics,
                               .id = 'train_test') %>%
        relocate(.config, starts_with('id'), contains('.iter'), 
                 .after = train_test) 
      
      write_csv(out_metrics, paste0(outdir, '/metrics/', top_fits$wflow_id, '_metrics.csv.gz'))
      
      collect_predictions(final_wflow) %>%
        write_csv(paste0(outdir, '/metrics/', top_fits$wflow_id, '_prediction.csv.gz'))
      
      nest(out_metrics, params = colnames(dplyr::select(best_params, -.config))) %>%
        mutate(wflow_id = top_fits$wflow_id, .before = 'train_test')
    }
    
    slurm_wflow_tune <- function(wflowset){
      registerDoFuture()
      plan(cluster, workers = 20)
      
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
      top_config <- rank_results(trained_wflow, rank_metric = 'mn_log_loss', select_best = TRUE) %>%
        pull(.config) %>%
        unique
      
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
      out <- workflow_post_process(trained_wflow, rank_metric = 'mn_log_loss', 
                                   train_test_split = coral_split, 
                                   outdir = file_out)
      
      } else {
        out <- tibble(wflow_id = trained_wflow$wflow_id[[1]])
      }
      
      out
    }
    
    reg <- makeRegistry(file.dir = paste0(getwd(), '/../batch_files'), 
                        packages = c('doFuture', 'tidyverse', 'tidymodels',
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
                       workflow_tuning = workflow_tuning,
                       workflow_post_process = workflow_post_process,
                       file_out = paste0(getwd(), '/../Results/model_tuning')))
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
                by = 'wflow_id')
    
    #plot top model metrics
    model_tune_plot <- disease_model_tune %>%
      
      group_by(wflow_id) %>%
      #get rid of failures
      filter(n() > 1) %>%
      filter(.config == .config[train_test == 'test']) %>%
      filter(train_test != 'test') %>%
      ungroup %>%
      dplyr::select(-params) %>%
      mutate(wflow_id = fct_reorder(wflow_id, mn_log_loss, .fun = mean, .na_rm = TRUE)) %>%
      pivot_longer(cols = -c(wflow_id:.iter),
                   names_to = '.metric') %>%
      filter(!is.na(value)) %>%
      group_by(wflow_id, .metric) %>%
      summarise(mean = mean(value),
                se = sd(value) / sqrt(n()),
                .groups = 'drop') %>%
      # arrange(wflow_id) %>% filter(.metric == 'mn_log_loss')
      
      ggplot(aes(x = as.integer(wflow_id), y = mean, ymin = mean - se, ymax = mean + se,
                 colour = wflow_id)) +
      geom_pointrange() +
      facet_wrap(~ .metric, scales = 'free_y') +
      theme_classic()
    ggsave('../Results/model_tuning_results.png', plot = model_tune_plot, height = 21, width = 21)
    
    #write csv/rds of all metrics
    write_rds(disease_model_tune, '../Results/model_tuning_metrics.rds.gz')
    
  } else {
    registerDoFuture()
    plan(cluster, workers = parallel::makeCluster(parallel::detectCores(logical = FALSE)))
    
    disease_model_tune <- workflow_tuning(disease_wflow)
    
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
    
    disease_model_tune <- remove_errors(disease_model_tune)
    
    # all_models_tune <- autoplot(disease_model_tune, metric = "bal_accuracy")
    model_tune_plot <- autoplot(disease_model_tune, select_best = TRUE, 
                                # metric = c('accuracy', 'bal_accuracy', "mn_log_loss"), 
                                rank_metric = "mn_log_loss")
    ggsave('../Results/model_tuning_results.png', plot = model_tune_plot, height = 21, width = 21)
    # autoplot(disease_model_tune, id = 'pls_qda', metric = "bal_accuracy")
    
    collect_metrics(disease_model_tune, summarize = FALSE)  %>%
      pivot_wider(names_from = '.metric',
                  values_from = '.estimate') %>%
      write_csv('../intermediate_files/model_tuning_metrics.csv.gz')
   
    write_rds(disease_model_tune, '../intermediate_files/model_tune_run.rds.gz') 
  }
}

#### Model Model Equivilance ####
if(file.exists('../intermediate_files/model_equivilance.rds.gz') & !rerun_process){
  model_equivilance_model <- read_rds('../intermediate_files/model_equivilance.rds.gz')
} else {
  
  if(Sys.info()['sysname'] != 'Windows'){
    top_stats_all_models <- disease_model_tune %>%
      group_by(wflow_id) %>%
      #get rid of failures
      filter(n() > 1) %>%
      filter(.config == .config[train_test == 'test']) %>%
      filter(train_test != 'test') %>%
      ungroup %>%
      dplyr::select(-params) %>%
      dplyr::select(wflow_id, starts_with('id'), bal_accuracy) %>%
      pivot_wider(names_from = 'wflow_id',
                  values_from = 'bal_accuracy')
    
  } else {
    get_best_fit_stats <- function(wflows, metric, rank_metric){
      top_fits <- rank_results(wflows, rank_metric = rank_metric, select_best = TRUE) %>%
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
    
    top_stats_all_models <- get_best_fit_stats(disease_model_tune, 'bal_accuracy', rank_metric = 'mn_log_loss')
  }

  model_class_tweak <- function(model, metric, direction = 'minimize'){
    class(model) <- c("perf_mod_workflow_set", "perf_mod")
    model$metric$name <- metric
    model$metric$direction <- direction
    model$transform <- logit_trans
    model
  }
  
  model_equivilance_model <- top_stats_all_models %>%
    mutate(across(where(is.numeric), ~if_else(. == 1, . - 1e-6, .)),
           across(where(is.numeric), ~binomial(link = 'logit')$linkfun(.))) %>%
    perf_mod(#formula = ,
      #transform = logit_trans,
      # family = 'binomial',
      hetero_var = FALSE,
      
      chains = 20,
      cores = 20,
      iter = 2000,
      warmup = 1000,
      refresh = 10) %>%
    model_class_tweak(metric = 'bal_accuracy', direction = 'maximize')
  
  write_rds(model_equivilance_model, '../intermediate_files/model_equivilance.rds.gz')
}

model_accuracies <- autoplot(model_equivilance_model, type = "intervals", prob = 0.9) + guides(colour = 'none') 
model_ropes <- autoplot(model_equivilance_model, size = 0.01, type = "ROPE") + #1% difference in accuracy
  geom_hline(yintercept = 1 - 0.05, linetype = 'dashed')

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

top_wflows <- find_equivilent_models(model_equivilance_model, 0.01, 0.05)

if(Sys.info()['sysname'] != 'Windows'){
  write_csv(top_wflows, '../Results/equivilant_top_models.csv.gz')
  
} else {
  top_disease_models <- top_wflows %>%
    rowwise %>%
    mutate(wflowset = list(extract_workflow_set_result(disease_model_tune, model)),
           # race_plot = list(plot_race(wflowset)),
           best_params = list(select_best(wflowset, 'mn_log_loss')),
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
  
}

#### Assess Model Accuracy on Test Dataset ####
if(Sys.info()['sysname'] == 'Windows'){
  top_disease_models %>%
    dplyr::select(model, bal_accuracy:mn_log_loss) %>%
    arrange(-bal_accuracy, mn_log_loss)
  
  top_disease_models  %>%
    arrange(-bal_accuracy) %>%
    dplyr::select(model, conf_mat) %$%
    set_names(conf_mat, model)
}


#### Write outputs ####
if(Sys.info()['sysname'] == 'Windows'){
  dplyr::select(top_disease_models, -mean:-upper, 
                -contains('race_plot'), -conf_mat) %>%
    ungroup %>%
    write_rds('../intermediate_files/top_disease_models.rds.gz')
  
}


#### Test-set metrics all models ####
if(Sys.info()['sysname'] == 'Windows'){
  test_set_metrics <- disease_model_tune %>%
    filter(!str_detect(wflow_id, 'null')) %>%
    rowwise %>%
    mutate(best_params = list(select_best(result, 'mn_log_loss')),
           final_wflow = list(extract_workflow(disease_model_tune, wflow_id) %>%
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
                         values_from = '.estimate')) %>%
    ungroup %>%
    dplyr::select(-where(is.list)) %>%
    arrange(mn_log_loss)
  write_csv(test_set_metrics, '../Results/test_set_metrics.csv')
  
}

#### Remove Directory ####
if(Sys.info()['sysname'] != 'Windows'){
  unlink(paste0(getwd(), '/../batch_files'), recursive = TRUE)
}
