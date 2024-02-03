# https://stackoverflow.com/questions/71662140/create-shap-plots-for-tidymodel-objects
# https://christophm.github.io/interpretable-ml-book/shap.html#shap-feature-importance

if(!interactive()){
  rerun_shap <- TRUE
  NSIM <- 1000
} else {
  rerun_shap <- FALSE
  NSIM <- 1
}

alpha <- 0.05
MAX_GROUP <- 1
min_test_accuracy <- 0.97
maxZ <- 2

# To run on Discovery 
# sbatch \
#   --dependency=afterany:40368753 \
#   --output=/scratch/j.selwyn/Panama_Tank_Field/slurm/shap_%j.out \
#   /work/vollmer/software/jds_scripts/runRscript.slurm \
#   /scratch/j.selwyn/Panama_Tank_Field/Code/4-Model_Explainations.R


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
library(patchwork)
library(shapviz)

if(Sys.info()['sysname'] != 'Windows'){
  setwd('/scratch/j.selwyn/Panama_Tank_Field/Code')
} else {
  setwd("~/Google Drive/Research/Vollmer Lab PostDoc/Panama_Tank_Field/Code")
}


#### Read in Models & Data ####
coral_split <- read_rds('../intermediate_files/coral_split.rds.gz')
taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz', 
                     show_col_types = FALSE) %>%
  mutate(across(everything(), str_replace_na))

tank_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                      show_col_types = FALSE) %>%
  filter(geno != 'GE', time == '8_exp',
         anti == 'N', health == 'D') %>%
  # mutate(sample_id = str_c(sample_id, time, sep = '_')) %>%
  select(asv_id, sample_id, health, log2_cpm_norm) %>%
  mutate(health = factor(health)) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm)

model_predictions <- list.files('../Results/model_tuning/metrics/',
           pattern = 'prediction.csv.gz$', full.names = TRUE) %>%
  tibble(pred_file = .) %>%
  mutate(wflow_id = str_extract(pred_file, '[a-zA-Z0-9_]+.csv.gz$') %>% 
           str_remove('_prediction.csv.gz$')) %>%
  rowwise(wflow_id) %>%
  reframe(read_csv(pred_file, show_col_types = FALSE)) %>%
  select(wflow_id, sample_id, .pred_class)


all_models <- list.files('../Results/model_tuning/models',
                         pattern = 'rds.gz$', full.names = TRUE) %>%
  R.utils::getAbsolutePath() %>%
  tibble(wflow_file = .) %>%
  mutate(wflow_id = str_extract(wflow_file, '[a-zA-Z0-9_]+.rds.gz$') %>% 
           str_remove('.rds.gz$'))

model_stats <- read_rds('../Results/model_tuning_metrics.rds.gz') %>%
  group_by(wflow_id) %>%
  filter(.config == .config[train_test == 'test']) %>%
  ungroup %>%
  select(-where(is.list)) %>%
  pivot_longer(cols = c(bal_accuracy:overall),
               names_to = '.metric') %>%
  group_by(wflow_id, train_test, .metric) %>%
  summarise(mean = mean(value),
            se = sd(value) / sqrt(n()),
            .groups = 'drop')

equivilant_models <- read_csv('../Results/equivilant_top_models.csv.gz', show_col_types = FALSE) %>%
  # filter(pract_equiv >= 0.95) %>%
  rename(wflow_id = model)

#### Join Models and arrange ####
top_disease_models <- inner_join(equivilant_models,
           all_models,
           by = 'wflow_id') %>%
  mutate(wflow_id = fct_reorder(wflow_id, mean, .desc = TRUE)) %>%
  select(-mean:-upper) %>%
  left_join(filter(model_stats, train_test == 'test') %>%
              select(-train_test, -se) %>%
              pivot_wider(names_from = '.metric',
                          values_from = 'mean'),
            by = 'wflow_id') %>%
  mutate(wflow_id = fct_reorder(wflow_id, overall, .desc = TRUE)) %>%
  mutate(wflow_id = fct_reorder(wflow_id, -pract_equiv)) %>%
  arrange(wflow_id) %>%
  arrange(-overall) %>%
  identity()

models_use <- top_disease_models %>%
  filter(pract_equiv >= 1 - 0.2) %>%
  pull(wflow_id)

#### Split up data into roughly equal sized chunks ####
full_data <- bind_rows(train = analysis(coral_split),
                       test = assessment(coral_split),
                       tank = tank_data,
                       .id = 'train_test') %>%
  mutate(grouping = row_number() %% MAX_GROUP, 
         .before = 'train_test') 

#### Make Model Predictions for everyone ####
# model_file <- 'C:/Users/jdsel/Documents/Google Drive/Research/Vollmer Lab PostDoc/Panama_Tank_Field/Results/model_tuning/models/pls_decision.rds.gz'
# get_pred_class <- function(model_file){
#   read_rds(model_file) %>%
#     fit(analysis(coral_split)) %>%
#     predict(new_data = full_data) %>%
#     bind_cols(select(full_data, sample_id), .)
# }
# 
# model_predictions <- top_disease_models %>%
#   rowwise(wflow_id) %>%
#   reframe(get_pred_class(wflow_file))

#### Calculate SHAP Values ####
if(file.exists('../intermediate_files/model_shaps.csv.gz') & !rerun_shap){
  top_shaps <- read_csv('../intermediate_files/model_shaps.csv.gz', show_col_types = FALSE) %>%
    mutate(wflow_id = factor(wflow_id, levels = levels(top_disease_models$wflow_id)))
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
              data = select(shap_sets$new_data[[1]], -grouping, -train_test),
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
             # test_train = shap_sets$test_train[[1]], 
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
  
  while(length(successful_jobs) < (nrow(top_disease_models) * MAX_GROUP) & z < maxZ){
    z <- z + 1
    
    # incomplete_shaps <- top_disease_models %>%
    #   expand_grid(test_train = c('test', 'train', 'tank')) %>%
    #   # mutate(test_train = 'train') %>%
    #   mutate(job_id = str_c(wflow_id, test_train, sep = '_')) %>%
    #   filter(!job_id %in% successful_jobs) %>%
    #   
    #   mutate(train_data = list(analysis(coral_split)),
    #          new_data = case_when(test_train == 'test' ~ list(assessment(coral_split)),
    #                               test_train == 'train' ~ list(analysis(coral_split)),
    #                               test_train == 'tank' ~ list(tank_data))) %>%
    #   select(job_id, wflow_id, wflow_file, train_data, test_train, new_data)
    
    # incomplete_shaps <- sample_n(incomplete_shaps, 5)
    
    incomplete_shaps <- top_disease_models %>%
      expand_grid(group_var = 0:(MAX_GROUP - 1)) %>% 
      # mutate(test_train = 'train') %>%
      mutate(job_id = str_c(wflow_id, group_var, sep = '_')) %>%
      filter(!job_id %in% successful_jobs) %>%
      rowwise %>%
      mutate(train_data = list(analysis(coral_split)),
             new_data = list(filter(full_data, grouping == group_var))) %>%
      select(job_id, wflow_id, wflow_file, train_data, group_var, new_data)
    
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
  
  top_shaps <- bind_rows(newly_out) %>%
    dplyr::select(-grouping)
  
  write_csv(top_shaps, '../intermediate_files/model_shaps.csv.gz')
}

#### Clean up SHAP results ####
individual_shap_values <- top_shaps %>%
  group_by(wflow_id, train_test) %>%
  mutate(across(ends_with('value'), ~scale(.)[,1])) %>%
  ungroup %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = c('asv_id', '.value'),
               names_pattern = '(.*)_(.*)')


#https://www.datarobot.com/blog/using-feature-importance-rank-ensembling-fire-for-advanced-feature-selection/
# shaps <- filter(individual_shap_values, wflow_id %in% models_use)
rank_asvs <- function(shaps, type = 'rank', adjust = 'fdr'){
  #Type must be 'rank' or 'shap' or 'model' which then models the rankings
  if(type == 'rank' | type == 'model'){
    direction <- 'less'
    wide_vals <- shaps %>%
      # filter(test_train == 'train') %>%
      group_by(wflow_id, asv_id) %>%
      summarise(shap_importance = mean(abs(shap)),
                .groups = 'drop') %>%
      pivot_wider(names_from = wflow_id, 
                  values_from = shap_importance) %>%
      mutate(across(where(is.numeric), ~rank(-.))) %>%
      rowwise %>%
      mutate(median_rank = median(c_across(where(is.numeric)))) %>%
      ungroup 
    
    filtered_wide_vals <- filter(wide_vals, median_rank != max(median_rank)) 
    
  } else {
    direction <- 'greater'
    wide_vals <- shaps %>%
      group_by(wflow_id, asv_id) %>%
      summarise(shap_importance = mean(abs(shap)),
                .groups = 'drop') %>%
      pivot_wider(names_from = wflow_id, 
                  values_from = shap_importance) %>%
      # mutate(across(where(is.numeric), ~rank(-.))) %>%
      rowwise %>%
      mutate(median_shap = median(c_across(where(is.numeric))),
             all_0 = all(c_across(where(is.numeric)) == 0)) %>%
      ungroup 
    
    filtered_wide_vals <- filter(wide_vals, !all_0) %>%
      select(-all_0)
  }
  
  if(type != 'model'){
    asv_tests <- filtered_wide_vals %>%
      pivot_longer(cols = -any_of(c('asv_id', 'median_rank', 'median_shap', 'all_0'))) %>%
      nest(data = -any_of(c('asv_id', 'median_rank', 'median_shap', 'all_0'))) %>%
      mutate(the_row = row_number()) %>%
      rowwise(any_of(c('asv_id', 'median_rank', 'median_shap', 'all_0'))) %>%
      mutate(other_data = list(bind_rows(.$data[-the_row]))) %>%
      
      summarise(tidy(wilcox.test(data$value, y = other_data$value, alternative = direction)),
                .groups = 'drop') %>%
      mutate(p_adjust = p.adjust(p.value, method = adjust), .after = 'p.value')
  } else {
    library(glmmTMB)
    library(emmeans)
    
    model_data <- filtered_wide_vals %>%
      select(-median_rank) %>%
      pivot_longer(cols = where(is.numeric),
                   names_to = 'model',
                   values_to = 'rank') 
    
    rank_model <- glmmTMB(rank ~ asv_id + (1 | model),
                          dispformula = ~asv_id,
                          data = model_data,
                          family = Gamma(link = 'log'),
                          REML = TRUE)
    
    asv_tests <- emmeans(rank_model, ~asv_id, type = 'response') %>%
      as_tibble() %>%
      left_join(emmeans(rank_model, ~asv_id) %>%
                  contrast(side = "<", adjust = adjust) %>%
                  as_tibble() %>%
                  mutate(asv_id = str_remove(contrast, ' effect')) %>%
                  select(asv_id, p.value),
                by = 'asv_id') %>%
      rename(p_adjust = p.value) %>%
      mutate(alternative = 'less')
  }
  
  
  full_join(select(wide_vals, -contains('all_0')),
            select(asv_tests, -any_of(c('median_rank', 'median_shap', 'all_0'))),
            by = 'asv_id') %>%
    arrange(p_adjust)
}

rank_plot <- function(asv_ranks){
  direction <- unique(na.omit(asv_ranks$alternative))
  
  asv_ranks %>%
    filter(!is.na(p_adjust)) %>%
    select(-any_of(c('response', 'SE', 'df', 'asymp.LCL', 'asymp.UCL'))) %>%
    pivot_longer(cols = -any_of(c('asv_id', 'median_rank', 'median_shap', 'statistic', 
                                  'p.value', 'p_adjust', 'method', 'alternative'))) %>%
    
    
    left_join(taxonomy, by = 'asv_id') %>%
    mutate(higher_taxonomy = str_c(family, genus, asv_id, sep = '; ') %>%
             str_remove('; NA')) %>%

    mutate(higher_taxonomy = fct_reorder(higher_taxonomy, value, 
                                         .desc = if_else(direction == 'less', TRUE, FALSE))) %>%
    
    # filter(median_rank <= 50) %>%
    
    ggplot(aes(x = value, y = higher_taxonomy, colour = p_adjust < alpha)) +
    # geom_point(size = 0.5) +
    # geom_boxplot(outlier.shape = NA) +
    # stat_summary(fun.data = median_hilow, show.legend = FALSE) +
    stat_summary(fun = median, fun.min = min, fun.max = max, show.legend = FALSE) +
    # stat_summary(fun.data = mean_se, show.legend = FALSE) +
    # stat_summary(fun.data = mean_cl_boot, show.legend = FALSE) +
    scale_colour_manual(values = c('FALSE' = 'black', 'TRUE' = 'red')) +
    labs(x = if_else(direction == 'less', 'ASV Rank', 'SHAP Importance'),
         y = NULL) +
    theme_classic() 
}

asv_ranks <- rank_asvs(filter(individual_shap_values, wflow_id %in% models_use), 
                       type = 'model', adjust = 'fdr')
write_csv(asv_ranks, '../Results/asv_importance.csv.gz')

rank_plots <- rank_plot(asv_ranks)
ggsave('../Results/asv_medianImportance.png', plot = rank_plots, height = 20, width = 7)

retain_asvs <- filter(asv_ranks, p_adjust < alpha) %>%
  pull(asv_id)

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
write_csv(model_important_asvs, '../intermediate_files/model_shaps_importance.csv.gz')

#### Overrepresentation of Families ####
row_fisher <- function(x, k, m, N, direction = 'two.sided'){
  ##https://dputhier.github.io/ASG/practicals/go_statistics_td/go_statistics_td_2015.html
  data.frame(significant = c(x, m - x),
             not_significant = c(k - x, N - m - (k - x))) %>%
    fisher.test(alternative = direction) %>%
    tidy
}

ora_asv <- function(asv_ranks, min_size = 5, lowest_level = 'family'){
  asv_ranks %>%
    mutate(p_adjust = if_else(is.na(p_adjust), 1, p_adjust)) %>%
    select(asv_id, p_adjust) %>%
    left_join(taxonomy, by = 'asv_id') %>%
    group_by(across(domain:!!sym(lowest_level))) %>%
    summarise(n_sig_fam = sum(p_adjust < alpha),
              total_fam = n(),
              .groups = 'rowwise') %>%
    mutate(total_sig = sum(.$n_sig_fam),
           total_asv = sum(.$total_fam)) %>%
    filter(total_fam >= min_size) %>%
    summarise(row_fisher(n_sig_fam, total_fam, total_sig, total_asv, direction = 'greater'),
              .groups = 'drop') %>%
    rename_with(.cols = -c(domain:!!sym(lowest_level)),
                ~str_c('ora', ., sep = '_'))
}

overrepresented_families <- ora_asv(asv_ranks, min_size = 5, lowest_level = 'family')

#### GSEA ####
gsea_asv <- function(asv_ranks, min_size = 5, lowest_level = 'family', the_test = wilcox.test){
  DIRECTION <-  unique(na.omit(asv_ranks$alternative))
  
  asv_ranks %>%
    select(asv_id, starts_with('median')) %>%
    left_join(taxonomy, by = 'asv_id') %>%
    rename(value = 2) %>%
    nest(data = -c(domain:!!sym(lowest_level))) %>%
    mutate(the_row = row_number()) %>%
    rowwise(domain:!!sym(lowest_level)) %>%
    mutate(other_data = list(bind_rows(.$data[-the_row]))) %>%
    filter(nrow(data) >= min_size) %>%
    summarise(the_test(data$value, other_data$value, alternative = DIRECTION) %>% tidy,
              .groups = 'drop') %>%
    rename_with(.cols = -c(domain:!!sym(lowest_level)),
                ~str_c('gsea', ., sep = '_'))
}

gsea_families <- gsea_asv(asv_ranks, min_size = 5, lowest_level = 'family', the_test = wilcox.test)

full_join(overrepresented_families,
          gsea_families) %>%
  write_csv('../Results/ASV_oraGSEA.csv.gz')

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
  # filter(asv_rank <= MAX_ASV) %>%
  filter(asv_id %in% retain_asvs,
         wflow_id %in% models_use) %>%
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
  # filter(asv_rank <= MAX_ASV) %>%
  filter(asv_id %in% retain_asvs,
         wflow_id %in% models_use) %>%
  mutate(higher_taxonomy = str_c(higher_taxonomy, asv_id, sep = '; '),
         higher_taxonomy = fct_reorder(higher_taxonomy, asv_rank)) %>% 
  rename(model = wflow_id) %>%
  
  ggplot(aes(x = model, y = shap_importance, 
             color = higher_taxonomy, 
             group = asv_id)) +
  
  
  
  geom_bump(show.legend = FALSE) +
  geom_point() +
  
  # geom_text(data = . %>% filter(model == levels(model)[1]),
  #           aes(label = asv_id, x = model),
  #           vjust = -1) +
  
  # coord_cartesian(ylim=c(MAX_ASV, 1)) +
  # scale_y_reverse(labels = 1:MAX_ASV, breaks = 1:MAX_ASV) +
  guides(colour = guide_legend(ncol = 3)) + 
  scale_x_discrete(labels = ~str_remove(., 'base_') %>% 
                     str_replace_all(c('forest' = 'Random\nForest',
                                       'knn' = 'KNN',
                                       'lasso' = 'LASSO',
                                       'mlp' = 'MLP',
                                       'null' = 'Null',
                                       'pls' = 'PLS',
                                       'svmLinear' = 'SVM'))) +
  labs(x = NULL,
       y = 'SHAP Importance',
       colour = NULL) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.text = element_text(colour = 'black', size = 8))
ggsave('../Results/Fig7_shap_importance_bump.png', plot = shap_importance_bump, height = 7, width = 11)

#### SHAP Summary ####
shap_summary_plot <- individual_shap_values %>% 
  inner_join(filter(model_important_asvs,
                    asv_id %in% retain_asvs,
                    wflow_id %in% models_use) %>%
               select(wflow_id, asv_id, 
                      higher_taxonomy, asv_rank), 
             by = c('wflow_id', 'asv_id')) %>%
  rename(model = wflow_id) %>%
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


#### SHAP Dependence Plots ####
# shapviz(object = shap_mat, X = value_mat) %>%
#   sv_dependence(v = 'ASV25')
# 
shap_dependence_plot <- function(shaps, asv, colour_var = 'health'){
  shaps %>%
    filter(asv_id == asv) %>%
    ggplot(aes(x = value, y = shap, colour = !!sym(colour_var), 
               shape = train_test)) +
    geom_point() +
    facet_grid( ~ wflow_id) +
    labs(x = str_c('log2cpm ', asv),
         y = 'SHAP Value',
         colour = str_to_sentence(colour_var)) +
    scale_size_manual(values = c('TRUE' = 1.5, 'FALSE' = 4), 
                      labels = c('TRUE' = 'Correct', 'FALSE' = 'Incorrect')) +
    scale_shape_manual(values = c('train' = 'circle', 'test' = 'square', 'tank' = 'triangle'),
                       labels = ~str_to_title(.),
                       breaks = c('train', 'test', 'tank')) +
    theme_classic() +
    theme(panel.background = element_rect(colour = 'black'),
          strip.background = element_blank())
}
# 
# shap_dependence_plot(individual_shap_values, asv = 'ASV25')

shap_dependence_plot <- individual_shap_values %>%
  inner_join(filter(asv_ranks, p_adjust < alpha) %>%
               select(asv_id, starts_with('median'), alternative) %>%
               rename(asv_rank = starts_with('median')) %>%
               mutate(asv_rank = if_else(alternative == 'greater', rank(-asv_rank), rank(asv_rank))) %>%
               select(-alternative),
             by = 'asv_id') %>%
  filter(!is.na(value)) %>%
  left_join(model_predictions,
            by = c('wflow_id', 'sample_id')) %>%
  filter(wflow_id %in% models_use) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(asv_name = str_c(order, family, genus, asv_id, sep = '\n'),
         asv_name = fct_reorder(asv_name, asv_rank)) %>%
  
  ggplot(aes(x = value, y = shap, colour = health, 
             shape = train_test, size = .pred_class == health)) +
  geom_point() +
  # facet_grid(wflow_id ~ asv_id, scales = 'free_x') +
  facet_grid(asv_name ~ wflow_id, scales = 'free_y') +
  labs(x = 'log2cpm',
       y = 'SHAP Value',
       colour = 'Disease\nState',
       size = 'Correct\nClassification',
       shape = 'Data\nSource') +
  scale_size_manual(values = c('TRUE' = 0.5, 'FALSE' = 2), 
                    labels = c('TRUE' = 'Correct', 'FALSE' = 'Incorrect')) +
  scale_shape_manual(values = c('train' = 'circle', 'test' = 'square', 'tank' = 'triangle'),
                     labels = ~str_to_title(.),
                     breaks = c('train', 'test', 'tank')) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank())
ggsave('../Results/shap_dependence.png', plot = shap_dependence_plot, height = 25, width = 25)

#### SHAP Force Plot ####
make_shap_force_plot <- function(shap_out, model_name, base_value, max_asvs = 15, collapse_to = NULL){
  #collapse by taxonomy
  if(!is.null(collapse_to)){
    shap_out <- shap_out %>%
      pivot_longer(cols = starts_with('ASV'),
                   names_to = c('asv_id', '.value'),
                   names_pattern = '(.*)_(.*)') %>%
      left_join(taxonomy, by = 'asv_id') %>%
      group_by(train_test, wflow_id, sample_id, health, !!sym(collapse_to)) %>%
      summarise(value = mean(value, na.rm = TRUE),
                shap = sum(shap, na.rm = TRUE),
                .groups = 'drop') %>%
      rename(taxa_var = collapse_to) %>%
      pivot_wider(names_from = taxa_var,
                  values_from = c('value', 'shap'),
                  names_glue = '{taxa_var}_{.value}')
  }
  
  #setup
  shap_mat <- shap_out %>%
    filter(wflow_id == model_name) %>%
    select(sample_id, ends_with('shap')) %>%
    rename_with(~str_remove(., '_shap')) %>%
    column_to_rownames('sample_id') %>%
    as.matrix()
  
  value_mat <- shap_out %>%
    filter(wflow_id == model_name) %>%
    select(sample_id, ends_with('value')) %>%
    rename_with(~str_remove(., '_value')) %>%
    column_to_rownames('sample_id') %>%
    as.matrix()
  
  shap_viz <- shapviz(object = shap_mat, X = value_mat, 
                      baseline = base_value)
  
  force_data <- shap_viz$S %>%
    as_tibble(rownames = 'sample_id') %>%
    # filter(sample_id == 'P6_Bin5_N_D_PI') %>%
    # filter(sample_id == 'P6_Bin2_N_D_YE') %>%
    # filter(sample_id == 'P6_Bin1_N_D_NE') %>%
    # filter(sample_id %in% c('P6_Bin5_N_D_PI', 'P6_Bin2_N_D_YE', 'P6_Bin1_N_D_NE')) %>%
    # filter(sample_id == 'P6_Bin2_N_D_NE') %>%
    left_join(filter(shap_out, wflow_id == model_name) %>%
                select(sample_id, health, train_test),
              by = 'sample_id')  %>%
    pivot_longer(cols = -c(sample_id, health, train_test),
                 names_to = 'asv_id',
                 values_to = 'shap_val') %>%
    
    # filter(health == 'H',
    #        train_test == 'tank') %>%
    
    group_by(sample_id, health, train_test) %>%
    arrange(-abs(shap_val), .by_group = TRUE) %>%
    mutate(a = row_number(-abs(shap_val))) %>%
    
    mutate(asv_id = case_when(a > max_asvs ~ "Other", 
                              TRUE ~ as.character(asv_id))) %>%
    group_by(sample_id, health, train_test, asv_id) %>%
    summarize(shap_val = sum(shap_val),
              .groups = 'drop') %>%
    mutate(direction = shap_val > 0) %>%
    ungroup %>%
    group_by(sample_id, health, train_test) %>%
    mutate(direction = case_when(sum(shap_val) > 0 ~ direction,
                                 TRUE ~ !direction)) %>%
    arrange(direction, -shap_val, .by_group = TRUE) %>%
    mutate(to = cumsum(shap_val),
           to = to + get_baseline(shap_viz),
           from = lag(to),
           from = if_else(is.na(from), to - shap_val, from)) %>%
    mutate(asv_id = str_remove(asv_id, '_asv_id')) %>%
    ungroup %>%
    mutate(across(c(to, from), ~case_when(. < 0 ~ 0,
                                          . > 1 ~ 1,
                                          TRUE ~ .)))
  
  force_data %>%
    # select(sample_id, health, train_test, asv_id, to, from) %>%
    pivot_longer(cols = c(to, from)) %>%
    ggplot(aes(y = sample_id, x = value)) +
    geom_path(aes(colour = asv_id,
                  group = interaction(direction, sample_id)),
              position = position_dodge(0.05)) +
    geom_point(data = force_data %>% #filter(sample_id == 'P6_Bin5_N_D_PI') %>%
                 group_by(train_test, sample_id, health) %>%
                 summarise(prob_d = sum(shap_val),
                           .groups = 'drop') %>%
                 mutate(prob_d = prob_d + get_baseline(shap_viz),
                        prob_d = case_when(prob_d < 0 ~ 0,
                                           prob_d > 1 ~ 1,
                                           TRUE ~ prob_d)),
               aes(x = prob_d)) +
    geom_vline(xintercept = shap_viz$baseline, linetype = 'dashed') +
    
    facet_grid(health + train_test ~ ., scales = 'free_y',
               labeller = labeller(health = c('D' = 'Disease', 'H' = 'Healthy'))) +
    # scale_x_continuous(labels = scales::percent_format(),
    #                    limits = c(0, 1)) +
    labs(y = NULL,
         x = 'Probability of Disease',
         colour = ifelse(is.null(collapse_to), 'ASV ID', str_to_title(collapse_to))) +
    guides(colour = guide_legend(title.position = 'top', title.hjust = 0.5, nrow = 5)) +
    theme_classic() +
    theme(legend.position = 'bottom',
          # axis.text = element_text(colour = 'black', size = 12),
          panel.background = element_rect(colour = 'black'),
          strip.background = element_blank()) 
  
}


force_plots <- select(top_disease_models, wflow_id, wflow_file) %>%
  filter(wflow_id %in% models_use) %>%
  expand_grid(correct_id = c(TRUE, FALSE)) %>%
  rename(model = wflow_id) %>% 
  # sample_n(2) %>%
  rowwise %>%
  mutate(data = list(left_join(top_shaps,
                               model_predictions,
                               by = c('wflow_id', 'sample_id')) %>%
                       filter(wflow_id == model, 
                              (health == .pred_class) == correct_id))) %>%
  left_join(group_by(model_predictions, wflow_id, .pred_class) %>%
              summarise(n = n(),
                        .groups = 'drop_last') %>%
              mutate(prop = n / sum(n)) %>%
              filter(.pred_class == 'D') %>%
              select(wflow_id, prop) %>%
              ungroup,
            by = c('model' = 'wflow_id')) %>%

  mutate(plot = list(make_shap_force_plot(data, model, prop, collapse_to = 'order') + 
                       labs(title = model)))

force_correct_plot <- force_plots %>%
  filter(correct_id) %>%
  pull(plot) %>%
  wrap_plots(nrow = 1) & 
  theme(legend.position = 'none') & 
  plot_annotation(title = 'Correct ID')
ggsave('../Results/all_force_plot_correctID.png', plot = force_correct_plot, width = 25, height = 10)

force_wrong_plot <- force_plots %>%
  filter(!correct_id) %>%
  pull(plot) %>%
  wrap_plots(nrow = 1) & 
  theme(legend.position = 'none') & 
  plot_annotation(title = 'Incorrect ID')
ggsave('../Results/all_force_plot_incorrectID.png', plot = force_wrong_plot, width = 25, height = 10)


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
