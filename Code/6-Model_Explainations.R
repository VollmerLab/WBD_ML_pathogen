##TODO - individual shap plots https://stackoverflow.com/questions/71662140/create-shap-plots-for-tidymodel-objects
##TODO - summary shap plots https://christophm.github.io/interpretable-ml-book/shap.html#shap-feature-importance

#TODO - split shaps across more nodes regardless of train/test stuff
#TODO - Sort out how to make force/waterfall plots
#TODO - fix "explain" function to use %dorng%

if(!interactive()){
  rerun_shap <- TRUE
  NSIM <- 100
} else {
  rerun_shap <- FALSE
  NSIM <- 1
}

MAX_ASV <- 10
min_test_accuracy <- 0.97
maxZ <- 2

# To run on Discovery 
# sbatch --output=/scratch/j.selwyn/Panama_Tank_Field/slurm/shap_%j.out \
#   /work/vollmer/software/jds_scripts/runRscript.slurm \
#   /scratch/j.selwyn/Panama_Tank_Field/Code/6-Model_Explainations.R


library(batchtools)
library(doFuture)
library(doRNG)
library(tidyverse)
library(magrittr)
library(tidymodels)
library(workflows)
library(workflowsets)
library(desirability2)
library(discrim)
library(themis)
library(fastshap)
library(ggbump)
library(ggbeeswarm)
library(tidytext)


if(Sys.info()['sysname'] != 'Windows'){
  setwd('/scratch/j.selwyn/Panama_Tank_Field/Code')
} 


#### Read in Models & Data ####
coral_split <- read_rds('../intermediate_files/coral_split.rds.gz')
taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz', show_col_types = FALSE) %>%
  mutate(across(everything(), str_replace_na))

tank_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                      show_col_types = FALSE) %>%
  filter(geno != 'GE', time == '8_exp',
         anti == 'N', health == 'D') %>%
  select(asv_id, sample_id, resist, log2_cpm_norm) %>%
  rename(health = resist) %>%
  mutate(health = str_replace_all(health, c('R' = 'H', 'S' = 'D')) %>% factor()) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm)

all_models <- list.files('../Results/model_tuning/models', pattern = 'rds.gz$', full.names = TRUE) %>%
  R.utils::getAbsolutePath() %>%
  tibble(wflow_file = .) %>%
  mutate(wflow_id = str_extract(wflow_file, '[a-zA-Z0-9_]+.rds.gz$') %>% str_remove('.rds.gz$'))

model_stats <- read_rds('../Results/model_tuning_metrics.rds.gz') %>%
  group_by(wflow_id) %>%
  filter(.config == .config[train_test == 'test']) %>%
  ungroup %>%
  select(-where(is.list)) %>%
  pivot_longer(cols = -c(wflow_id:.iter),
               names_to = '.metric') %>%
  group_by(wflow_id, train_test, .metric) %>%
  summarise(mean = mean(value),
            se = sd(value) / sqrt(n()),
            .groups = 'drop')

equivilant_models <- read_csv('../Results/equivilant_top_models.csv.gz', show_col_types = FALSE) %>%
  # filter(pract_equiv == 1) %>%
  rename(wflow_id = model)

#### Join Models and arrange ####
top_disease_models <- inner_join(equivilant_models,
           all_models,
           by = 'wflow_id') %>%
  select(-mean:-pract_equiv) %>%
  left_join(filter(model_stats, train_test == 'test') %>%
              select(-train_test, -se) %>%
              pivot_wider(names_from = '.metric',
                          values_from = 'mean'),
            by = 'wflow_id') %>%
  
  mutate(bacc_d = d_max(bal_accuracy, use_data = TRUE),
         loss_d = d_min(mn_log_loss, use_data = TRUE),
         brier_d = d_min(brier_class, use_data = TRUE),
         overall = d_overall(across(ends_with("_d")))) %>%
  dplyr::select(-ends_with('_d')) %>%
  mutate(wflow_id = fct_reorder(wflow_id, overall)) %>%
  arrange(overall) 

#### Calculate SHAP Values ####
if(file.exists('../intermediate_files/model_shaps.csv.gz') & !rerun_shap){
  top_shaps <- read_csv('../intermediate_files/model_shaps.csv.gz', show_col_types = FALSE)
} else {
  predict_fun <- function(object, newdata) {
    # library(tidymodels); ; library(workflowsets)
    library(workflows)
    predict(object, new_data = newdata, type = 'prob') %>% 
      pull(.pred_D)
  }
  
  calc_shap <- function(wflow, data, sims = 1, features = NA, ...){
    model <- wflow
    
    
    if(all(is.na(features))){
      features <- colnames(select(data, where(is.numeric)))
    }
    `%dopar%` <- `%dorng%`
    model_shap <- explain(model, 
                          feature_names = features,
                          X = data,
                          pred_wrapper = predict_fun,
                          nsim = sims,
                          adjust = sims > 1,
                          ...)
    
    model_shap
  }
  
  slurm_calc_shap <- function(shap_sets){
    # registerDoRNG()
    registerDoFuture()
    plan(cluster, workers = 20)
    
    fit_wflow <- read_rds(shap_sets$wflow_file[1]) %>%
      fit(shap_sets$train_data[[1]])
    
    calc_shap(wflow = fit_wflow,
              data = shap_sets$new_data[[1]],
              sims = NSIM,
              parallel = TRUE) %>%
      as_tibble() %>%
      mutate(across(everything(), as.numeric)) %>%
      # bind_cols(select(new_data, sample_id, health), .) %>%
      rename_with(.cols = starts_with('ASV'),
                  ~str_c(., '_shap')) %>%
      bind_cols(rename_with(shap_sets$new_data[[1]], 
                            .cols = starts_with('ASV'),
                            ~str_c(., '_value')), 
                .) %>%
      mutate(wflow_id = shap_sets$wflow_id[[1]], 
             test_train = shap_sets$test_train[[1]], 
             .before = sample_id)
  }
  
  reg <- makeRegistry(file.dir = paste0(getwd(), '/../batch_files_shap'), 
                      packages = c('doFuture', 'doRNG',
                                   'tidyverse', 'tidymodels',
                                   'discrim', 'bonsai', 'baguette',
                                   'plsmod', 'rules', 'themis', 'finetune',
                                   'fastshap'))
  
  reg$cluster.functions <- makeClusterFunctionsSlurm(template = "~/slurm_template.tmpl",
                                                     array.jobs = TRUE,
                                                     nodename = "localhost",
                                                     scheduler.latency = 1,
                                                     fs.latency = 65)
  
  successful_jobs <- character()
  z <- 0
  newly_out <- vector('list', maxZ)
  # walk(paste0(getwd(), '/../Results/', c('plots', 'metrics', 'models')),
  #      ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))
  
  while(length(successful_jobs) < (nrow(top_disease_models) * 3) & z < maxZ){
    z <- z + 1
    
    incomplete_shaps <- top_disease_models %>%
      expand_grid(test_train = c('test', 'train', 'tank')) %>%
      # mutate(test_train = 'train') %>%
      mutate(job_id = str_c(wflow_id, test_train, sep = '_')) %>%
      filter(!job_id %in% successful_jobs) %>%
      
      mutate(train_data = list(analysis(coral_split)),
             new_data = case_when(test_train == 'test' ~ list(assessment(coral_split)),
                                  test_train == 'train' ~ list(analysis(coral_split)),
                                  test_train == 'tank' ~ list(tank_data))) %>%
      select(job_id, wflow_id, wflow_file, train_data, test_train, new_data)
    
    # incomplete_shaps <- sample_n(incomplete_shaps, 5)
    
    if(z > 1){
      message('Running Job IDs: ', str_c(incomplete_shaps$job_id, collapse = '; '))
    }
    
    run_shaps <- split(incomplete_shaps, 1:nrow(incomplete_shaps))
    
    
    batchMap(fun = slurm_calc_shap, shap_sets = run_shaps)
    batchExport(list(predict_fun = predict_fun,
                     calc_shap = calc_shap,
                     NSIM = NSIM))
    submitJobs(resources = list(max.concurrent.jobs = 30))
    waitForJobs()
    
    newly_finished <- dplyr::as_tibble(getJobTable()) %>%
      dplyr::filter(is.na(error) & !is.na(done)) %>%
      dplyr::select(job.pars) %>%
      unnest(job.pars) %>%
      unnest(job.pars)
    
    if(nrow(newly_finished) > 0){
      newly_finished <- pull(newly_finished, job_id)
      newly_out[[z]] <- purrr::map_dfr(which(incomplete_shaps$job_id %in% newly_finished), loadResult)
      successful_jobs <- c(successful_jobs, newly_finished)
    } 
    
    clearRegistry()
  }
  
  top_shaps <- bind_rows(newly_out)
  
  write_csv(top_shaps, '../intermediate_files/model_shaps.csv.gz')
}

#### Clean up SHAP results ####
individual_shap_values <- top_shaps %>%
  group_by(wflow_id, test_train) %>%
  mutate(across(ends_with('value'), ~scale(.)[,1])) %>%
  ungroup %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = c('asv_id', '.value'),
               names_pattern = '(.*)_(.*)')

model_important_asvs <- individual_shap_values %>%
  # filter(test_train == 'train') %>%
  group_by(wflow_id, asv_id) %>%
  summarise(shap_importance = mean(abs(shap)),
            .groups = 'drop') %>%
  group_by(wflow_id) %>%
  mutate(asv_rank = rank(-shap_importance)) %>%
  ungroup %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(higher_taxonomy = str_c(family, genus, sep = '; ') %>%
           str_remove('; NA')) %>%
  select(wflow_id, asv_id, shap_importance, asv_rank, higher_taxonomy)

#To use SHAPviz
# library(shapviz)
# 
# shap_mat <- top_shaps %>%
#   filter(model == 'pca_lda') %>%
#   select(sample_id, ends_with('shap')) %>%
#   rename_with(~str_remove(., '_shap')) %>%
#   column_to_rownames('sample_id') %>%
#   as.matrix()
# 
# value_mat <- top_shaps %>%
#   filter(model == 'pca_lda') %>%
#   select(sample_id, ends_with('value')) %>%
#   rename_with(~str_remove(., '_value')) %>%
#   column_to_rownames('sample_id') %>%
#   as.matrix()
# 
# shapviz(object = shap_mat, X = value_mat) %>%
#   sv_importance(kind = 'bar')

#### SHAP Importance Plots ####
shap_importance <- model_important_asvs %>% 
  filter(asv_rank <= MAX_ASV) %>%
  mutate(higher_taxonomy = str_c(asv_id, higher_taxonomy, sep = ': '),
         higher_taxonomy = reorder_within(higher_taxonomy, -asv_rank, wflow_id)) %>%
  ggplot(aes(x = shap_importance, y = higher_taxonomy)) +
  geom_vline(xintercept = 0) +
  geom_segment(xend = 0, aes(yend = higher_taxonomy), linetype = 'dashed') +
  geom_point() +
  facet_wrap(~wflow_id, scales = 'free_y') +
  scale_y_reordered()
ggsave('../Results/shap_importance_facets.png', plot = shap_importance, height = 21, width = 21)

#### SHAP Importance - sankey ####
shap_importance_bump <- model_important_asvs %>% #count(asv_rank)
  filter(asv_rank <= MAX_ASV) %>%
  mutate(higher_taxonomy = fct_reorder(higher_taxonomy, asv_rank)) %>% 
  rename(model = wflow_id) %>%
  
  ggplot(aes(x = model, y = shap_importance, 
             color = higher_taxonomy, 
             group = asv_id)) +
  
  
  
  geom_bump() +
  geom_point() +
  
  # geom_text(data = . %>% filter(model == levels(model)[1]),
  #           aes(label = asv_id, x = model),
  #           vjust = -1) +
  
  # coord_cartesian(ylim=c(MAX_ASV, 1)) +
  # scale_y_reverse(labels = 1:MAX_ASV, breaks = 1:MAX_ASV) +
  guides(colour = guide_legend(ncol = 3)) + 
  labs(x = NULL,
       y = 'SHAP Importance Ranking',
       colour = NULL) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')
ggsave('../Results/shap_importance_bump.png', plot = shap_importance_bump, height = 7, width = 10)

#### SHAP Family Importance ####
individual_shap_values %>%
  # filter(test_train == 'train') %>%
  # left_join(taxonomy, by = c('asv_id')) %>%
  group_by(wflow_id, asv_id) %>%
  summarise(shap_importance = mean(abs(shap)),
            .groups = 'drop') %>%
  left_join(taxonomy, by = 'asv_id') %>%
  group_by(wflow_id, family) %>%
  summarise(shap_importance = sum(shap_importance),
            n_asv = n_distinct(asv_id),
            .groups = 'drop_last') %>%
  mutate(taxa_rank = rank(-shap_importance)) %>%
  ungroup %>%
  select(wflow_id, family, shap_importance, taxa_rank, n_asv) %>%
  
  filter(taxa_rank <= 5) %>%
  mutate(family = fct_reorder(family, taxa_rank)) %>% 
  rename(model = wflow_id) %>%
  
  ggplot(aes(x = model, y = shap_importance, 
             color = family, 
             group = family)) +
  
  
  
  geom_bump() +
  geom_point(aes(size = n_asv)) +
  
  # geom_text(data = . %>% filter(model == levels(model)[1]),
  #           aes(label = asv_id, x = model),
  #           vjust = -1) +
  
  # coord_cartesian(ylim=c(MAX_ASV, 1)) +
  # scale_y_reverse(labels = 1:MAX_ASV, breaks = 1:MAX_ASV) +
  guides(colour = guide_legend(ncol = 3)) + 
  labs(x = NULL,
       y = 'SHAP Importance Ranking',
       colour = NULL) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

#### SHAP Summary ####
shap_summary_plot <- individual_shap_values %>% 
  inner_join(filter(model_important_asvs,
                    asv_rank <= MAX_ASV) %>%
               select(wflow_id, asv_id, 
                      higher_taxonomy, asv_rank), 
             by = c('wflow_id', 'asv_id')) %>%
  rename(model = wflow_id) %>%
  # filter(model == 'pca_lda') %>%
  mutate(higher_taxonomy = str_c(asv_id, higher_taxonomy, sep = ': '),
         higher_taxonomy = reorder_within(higher_taxonomy, -asv_rank, model)) %>%
  
  ggplot(aes(y = higher_taxonomy, x = shap, colour = abs(value))) +
  geom_vline(xintercept = 0) +
  geom_quasirandom(groupOnX = FALSE, dodge.width = 0.1) +
  facet_wrap(~model,
             scales = 'free') +
  scale_y_reordered() +
  scale_color_distiller(type = 'div', palette = 'Spectral') +
  labs(y = NULL)
ggsave('../Results/shap_summary_facets.png', plot = shap_summary_plot, height = 21, width = 21)

#### Directional Shap Importance ####
individual_shap_values %>%
  group_by(wflow_id, health, asv_id) %>%
  summarise(shap = mean(shap),
            .groups = 'drop') %>%
  inner_join(filter(model_important_asvs, asv_rank <= MAX_ASV),
             by = c('wflow_id', 'asv_id')) %>%
  mutate(higher_taxonomy = fct_reorder(higher_taxonomy, asv_rank)) %>% 
  rename(model = wflow_id) %>%
  
  ggplot(aes(x = model, y = shap, 
             color = higher_taxonomy, 
             group = interaction(asv_id, health))) +
  
  
  
  geom_bump() +
  geom_point() +
  
  # geom_text(data = . %>% filter(model == levels(model)[1]),
  #           aes(label = asv_id, x = model),
  #           vjust = -1) +
  
  # coord_cartesian(ylim=c(MAX_ASV, 1)) +
  # scale_y_reverse(labels = 1:MAX_ASV, breaks = 1:MAX_ASV) +
  guides(colour = guide_legend(ncol = 3)) + 
  labs(x = NULL,
       y = 'SHAP Importance Ranking',
       colour = NULL) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

#### SHAP Dependence Plots ####
# shapviz(object = shap_mat, X = value_mat) %>%
#   sv_dependence(v = 'ASV25')
# 
# shap_dependence_plot <- function(shaps, asv, colour_var = 'health'){
#   shaps %>%
#     filter(asv_id == asv) %>%
#     ggplot(aes(x = value, y = shap, colour = !!sym(colour_var))) +
#     geom_point() +
#     facet_wrap(~model) +
#     labs(x = str_c('log2cpm ', asv),
#          y = 'SHAP Value',
#          colour = str_to_sentence(colour_var)) +
#     theme_classic()
# }
# 
# shap_dependence_plot(individual_shap_values, asv = 'ASV30')

#### SHAP Force Plot ####
# shap_force_plot <- function(shaps, max_asv, plot_samples = 'test'){
#   individual_shap_values %>%
#     filter(test_train == 'test') %>%
#     filter(sample_id == '2016_S_CK14_D_314') %>%
#     filter(model == 'base_rda') %>%
#     arrange(-abs(shap))
# }
# 
# shapviz(object = shap_mat, X = value_mat) %>%
#   sv_force(row_id = 109L, max_display = 10)

#### Plot important ASVs by health ####
# training(coral_split) 
# testing(coral_split) %>%
#   ggplot(aes(x = health, y = ASV700)) +
#   geom_boxplot()
# 
# 
# model_important_asvs
# 
# 
# filter(model_important_asvs,
#        asv_rank <= MAX_ASV) %>%
#   arrange(asv_rank) %>%
#   filter(model == 'base_rda')

#### Remove Directory ####
if(Sys.info()['sysname'] != 'Windows'){
  unlink(paste0(getwd(), '/../batch_files_shap'), recursive = TRUE)
}
