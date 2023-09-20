##TODO - individual shap plots https://stackoverflow.com/questions/71662140/create-shap-plots-for-tidymodel-objects
##TODO - summary shap plots https://christophm.github.io/interpretable-ml-book/shap.html#shap-feature-importance

#TODO - run with many bootstraps
#TODO - figure out problem with pls_svm

library(tidyverse)
library(magrittr)
library(tidymodels)
library(discrim)
library(themis)
library(fastshap)
library(ggbump)
library(ggbeeswarm)
library(tidytext)
library(multidplyr)

rerun_shap <- TRUE
NSIM <- 2
MAX_ASV <- 20

#### Read in Models & Data ####
coral_split <- read_rds('../intermediate_files/coral_split.rds.gz')
top_disease_models <- read_rds('../intermediate_files/top_disease_models.rds.gz') %>%
  mutate(model = fct_reorder(model, -bal_accuracy))
taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz', show_col_types = FALSE) %>%
  mutate(across(everything(), str_replace_na))

#### Calculate SHAP Values ####
if(file.exists('../intermediate_files/model_shaps.csv.gz') & !rerun_shap){
  top_shaps <- read_csv('../intermediate_files/model_shaps.csv.gz', show_col_types = FALSE)
} else {
  predict_fun <- function(object, newdata) {
    predict(object, new_data = newdata, type = 'prob') %>% 
      pull(.pred_D)
  }
  calc_shap <- function(wflow, data, sims = 1, features = NA, ...){
    model <- extract_workflow(wflow)
    
    if(all(is.na(features))){
      features <- colnames(select(data, where(is.numeric)))
    }
    
    model_shap <- explain(model, 
                          feature_names = features,
                          X = data,
                          pred_wrapper = predict_fun,
                          nsim = sims,
                          adjust = sims > 1,
                          ...)
    
    model_shap
  }
  
  cluster <- new_cluster(parallel::detectCores() - 1)
  cluster_library(cluster, c('tidymodels', 'fastshap', 'dplyr', 'stringr'))
  cluster_copy(cluster, c('predict_fun', 'calc_shap', 'NSIM'))
  
  top_shaps <- top_disease_models %>%
    expand_grid(test_train = c('test', 'train')) %>%
    
    mutate(new_data = case_when(test_train == 'test' ~ list(analysis(coral_split)),
                                test_train == 'train' ~ list(assessment(coral_split)))) %>%
    
    # filter(test_train == 'test') %>%
    # sample_n(4) %>%
    
    rowwise(test_train, model, bal_accuracy) %>%
    partition(cluster) %>%
    summarise(calc_shap(final_wflow,
                        data = new_data,
                        # features = c('ASV8', 'ASV25', 'ASV700'),
                        sims = NSIM) %>%
                as_tibble() %>%
                mutate(across(everything(), as.numeric)) %>%
                # bind_cols(select(new_data, sample_id, health), .) %>%
                rename_with(.cols = starts_with('ASV'),
                            ~str_c(., '_shap')) %>%
                bind_cols(rename_with(new_data, 
                                      .cols = starts_with('ASV'),
                                      ~str_c(., '_value')), 
                          .)) %>%
    
    # mutate(shap_out = list(calc_shap(final_wflow,
    #                                  data = new_data,
    #                                  sims = NSIM,
    #                                  shap_only = FALSE))) %>%
    collect() %>%
    # select(model, bal_accuracy, shap_out) %>%
    ungroup
  write_csv(top_shaps, '../intermediate_files/model_shaps.csv.gz')
  rm(cluster); gc()
}

#### Clean up SHAP results ####
individual_shap_values <- top_shaps %>%
  # filter(test_train == 'train') %>%
  mutate(across(ends_with('value'), ~scale(.)[,1])) %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = c('asv_id', '.value'),
               names_pattern = '(.*)_(.*)')

model_important_asvs <- individual_shap_values %>%
  group_by(model, asv_id) %>%
  summarise(shap_importance = mean(abs(shap)),
            .groups = 'drop') %>%
  group_by(model) %>%
  mutate(asv_rank = rank(-shap_importance)) %>%
  ungroup %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(higher_taxonomy = str_c(family, genus, sep = '; ') %>%
           str_remove('; NA')) %>%
  select(model, asv_id, shap_importance, asv_rank, higher_taxonomy)

#To use SHAPviz
library(shapviz)

shap_mat <- top_shaps %>%
  filter(model == 'pca_lda') %>%
  select(sample_id, ends_with('shap')) %>%
  rename_with(~str_remove(., '_shap')) %>%
  column_to_rownames('sample_id') %>%
  as.matrix()

value_mat <- top_shaps %>%
  filter(model == 'pca_lda') %>%
  select(sample_id, ends_with('value')) %>%
  rename_with(~str_remove(., '_value')) %>%
  column_to_rownames('sample_id') %>%
  as.matrix()

shapviz(object = shap_mat, X = value_mat) %>%
  sv_importance(kind = 'bar')

#### SHAP Importance Plots ####
model_important_asvs %>% 
  filter(asv_rank <= MAX_ASV) %>%
  mutate(higher_taxonomy = str_c(asv_id, higher_taxonomy, sep = ': '),
         higher_taxonomy = reorder_within(higher_taxonomy, -asv_rank, model)) %>%
  ggplot(aes(x = shap_importance, y = higher_taxonomy)) +
  geom_vline(xintercept = 0) +
  geom_segment(xend = 0, aes(yend = higher_taxonomy), linetype = 'dashed') +
  geom_point() +
  facet_wrap(~model, scales = 'free_y') +
  scale_y_reordered()

#### SHAP Importance - sankey ####
model_important_asvs %>% #count(asv_rank)
  filter(asv_rank <= MAX_ASV) %>%
  
  ggplot(aes(x = model, y = shap_importance, 
             color = higher_taxonomy, 
             group = asv_id)) +
  geom_bump() +
  geom_point() +
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
individual_shap_values %>% 
  inner_join(filter(model_important_asvs,
                    asv_rank <= MAX_ASV) %>%
               select(model, asv_id, 
                      higher_taxonomy, asv_rank), 
             by = c('model', 'asv_id')) %>%
  # filter(model == 'pca_lda') %>%
  mutate(higher_taxonomy = str_c(asv_id, higher_taxonomy, sep = ': '),
         higher_taxonomy = reorder_within(higher_taxonomy, -asv_rank, model)) %>%
  
  ggplot(aes(y = higher_taxonomy, x = shap, colour = abs(value))) +
  geom_vline(xintercept = 0) +
  geom_quasirandom(groupOnX = FALSE, dodge.width = 0.1) +
  facet_wrap(~model,
             scales = 'free_y') +
  scale_y_reordered() +
  scale_color_distiller(type = 'div', palette = 'Spectral') +
  labs(y = NULL)
  
#### SHAP Dependence Plots ####
shapviz(object = shap_mat, X = value_mat) %>%
  sv_dependence(v = 'ASV25')

#### SHAP Force Plot ####
shapviz(object = shap_mat, X = value_mat) %>%
  sv_force(row_id = 109L, max_display = 10)
