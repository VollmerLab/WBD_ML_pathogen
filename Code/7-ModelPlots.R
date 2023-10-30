##TODO - add prediction contours
##TODO - add classification probability into plots
##TODO - annotate plots with accuracy of each of the 3 groups
##TODO - why is lasso always equal quality...

#### Libraries ####
library(tidyverse)
library(tidymodels)
library(tidyposterior)
library(discrim)
library(plsmod)
library(ggrepel)
library(patchwork)

#### Data ####
model_dir <- '../Results/model_tuning/models'

coral_split <- read_rds('../intermediate_files/coral_split.rds.gz')
tank_data_wide <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                           show_col_types = FALSE) %>%
  filter(geno != 'GE', time == '8_exp',
         anti == 'N', health == 'D') %>%
  select(asv_id, sample_id, resist, log2_cpm_norm) %>%
  rename(health = resist) %>%
  mutate(health = str_replace_all(health, c('R' = 'H', 'S' = 'D')) %>% factor()) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm)

taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz')
model_list <- read_csv('../Results/equivilant_top_models.csv.gz', show_col_types = FALSE) %>%
  filter(pract_equiv >= 0.95) %>%
  rename(wflow_id = model) %>%
  pull(wflow_id)

shap_importance <- read_csv('../intermediate_files/model_shaps_importance.csv.gz', 
                            show_col_types = FALSE)

#Plot for pls_lda - want to plot decision boundary - https://stats.stackexchange.com/questions/92157/compute-and-graph-the-lda-decision-boundary
#https://stackoverflow.com/questions/30619616/how-to-plot-classification-borders-on-an-linear-discrimination-analysis-plot-in
plot_nonBase <- function(model_name, x = 'PLS1', y = 'PLS2', N_GRID = 100){
  fit_model <- str_c(model_dir, '/', model_name, '.rds.gz') %>%
    read_rds() %>%
    fit(analysis(coral_split))
  
  prepped_model <- extract_preprocessor(fit_model) %>%
    prep()
  
  if(str_detect(model_name, 'pca')){
    pct_exp <- round(100 * (prepped_model$steps[[6]]$res$sdev^2 / sum(prepped_model$steps[[6]]$res$sdev ^ 2)), 1)
    x_name <- str_c(x, ' (', pct_exp[str_extract(x, '[0-9]+') %>% as.numeric()], '%)')
    y_name = str_c(y, ' (', pct_exp[str_extract(y, '[0-9]+') %>% as.numeric()], '%)')
  } else {
    x_name <- x; y_name = y
  }
  
  prediction_data <- bind_rows(train = prepped_model %>%
                                 juice() %>%
                                 bind_cols(predict(fit_model, new_data = analysis(coral_split))),
                               
                               test = prepped_model %>%
                                 bake(new_data = assessment(coral_split)) %>%
                                 bind_cols(predict(fit_model, new_data = assessment(coral_split))), 
                               
                               tank = prepped_model %>%
                                 bake(new_data = tank_data_wide) %>%
                                 bind_cols(predict(fit_model, new_data = tank_data_wide)),
                               
                               .id = 'data_source') %>%
    mutate(match_mismatch = health == .pred_class)
  
  
  #WORKING HERE TO ADD CONTORS
  extra_prediction_columns <- str_subset(colnames(prediction_data), 
                                         str_remove_all(c(x, y), '[0-9]+') %>% 
                                           unique) %>%
    str_subset(pattern = str_c(x, y, sep = '|'), negate = TRUE)
  
  
  prediction_grid <- expand_grid(!!sym(x) := modelr::seq_range(prediction_data[[x]], n = N_GRID, expand = 0.1),
                                 !!sym(y) := modelr::seq_range(prediction_data[[y]], n = N_GRID, expand = 0.1))
  if(length(extra_prediction_columns) > 0){
    for(variable in extra_prediction_columns){
      prediction_grid <- mutate(prediction_grid, !!sym(variable) := median(prediction_data[[variable]]))
    }
  }
  
  #lda
  if(str_detect(model_name, 'lda')){
    pred_grid <- extract_fit_engine(fit_model) %>%
      predict(newdata = prediction_grid) %>%
      bind_cols %>%
      select(class) %>%
      bind_cols(prediction_grid)
  }
  #decision
  if(str_detect(model_name, 'decision')){
    pred_grid <- extract_fit_engine(fit_model) %>%
      predict(newdata = prediction_grid) %>%
      apply(1, which.max) %>%
      bind_cols(prediction_grid, class = .) %>%
      mutate(class = if_else(class == 1, 'D', 'H') %>% factor)
  }
  #forest
  if(str_detect(model_name, 'forest')){
    pred_grid <- extract_fit_engine(fit_model) %>%
      predict(data = prediction_grid) %>% 
      .$predictions %>% 
      apply(1, which.max) %>%
      bind_cols(prediction_grid, class = .) %>%
      mutate(class = if_else(class == 1, 'D', 'H') %>% factor)
  }
  #svm
  if(str_detect(model_name, 'svm')){
    library(kernlab)
    pred_grid <- extract_fit_engine(fit_model) %>%
      predict(newdata = prediction_grid) %>%
      bind_cols(prediction_grid, class = .)
  }
  
  out <- prepped_model %>%
    tidy(6) %>%
    select(-id) %>%
    mutate(value = 50 * value) %>%
    pivot_wider(names_from = component, 
                values_from = value) %>%
    left_join(filter(shap_importance, wflow_id == model_name) %>%
                select(-wflow_id, -higher_taxonomy),
              by = c('terms' = 'asv_id')) %>%
    left_join(taxonomy, by = c('terms' = 'asv_id')) %>% 
    mutate(genus = if_else(is.na(genus), family, genus)) %>%
    
    ggplot(aes(x = !!sym(x), y = !!sym(y))) +
    geom_raster(data = pred_grid, 
                aes(fill = factor(class)),
                alpha = 0.5,
                show.legend = FALSE) +
    geom_segment(xend = 0, yend = 0, 
                 aes(alpha = 1/asv_rank)) +
    geom_text(data = . %>% filter(asv_rank <= 5),
                    aes(label = str_c(genus, terms, sep = '; '),
                        alpha = 1/asv_rank),
                    show.legend = FALSE) +
    
    geom_point(data = prediction_data,
               aes(colour = health, size = match_mismatch, shape = data_source)) +
    
    scale_size_manual(values = c('TRUE' = 1.5, 'FALSE' = 4), 
                      labels = c('TRUE' = 'Correct', 'FALSE' = 'Incorrect')) +
    scale_alpha_continuous('ASV\nImportance', breaks = c(1/(1:5)), labels = c(1:5)) +
    scale_shape_manual(values = c('train' = 'circle', 'test' = 'square', 'tank' = 'triangle'),
                       labels = ~str_to_title(.),
                       breaks = c('train', 'test', 'tank')) +
    labs(x = x_name,
         y = y_name,
         size = 'Correct\nClassification',
         shape = 'Data\nSource',
         colour = 'Disease\nState') +
    theme_classic()
  
  if(str_detect(model_name, 'lda')){
    out <- out + 
      geom_point(data = fit_model$fit$fit$fit$means %>%
                   as_tibble(rownames = 'health'),
                 size = 10, aes(fill = health),
                 shape = 21, show.legend = FALSE)
  }
  out
}
plot_nonBase('pls_lda', x = 'PLS1', y = 'PLS2', N_GRID = 1000)
plot_nonBase('pca_lda', x = 'PC1', y = 'PC2')

model_list %>%
  str_subset('base', negate = TRUE) %>%
  tibble(model = .) %>%
  rowwise %>%
  mutate(plot = list(plot_nonBase(model, N_GRID = 1000) + labs(title = model))) %>%
  pull(plot) %>%
  wrap_plots() + 
  plot_layout(guides = 'collect')
ggsave('../Results/topPLSmodels.png', height = 15, width = 15)

## Lasso
model_name <- 'base_lasso'
fit_model <- str_c(model_dir, '/', model_name, '.rds.gz') %>%
  read_rds() %>%
  fit(analysis(coral_split))
lambda <- fit_model %>%
  extract_fit_parsnip() %>%
  tidy %>%
  pull(penalty) %>%
  unique

read_rds('../Results/model_tuning_metrics.rds.gz') %>%
  filter(wflow_id == 'base_lasso') %>%
  filter(train_test == 'train') %>%
  unnest(params) %>%
  group_by(penalty) %>%
  summarise(brier_class  = mean(brier_class )) %>%
  ggplot(aes(x = penalty, y = brier_class)) +
  geom_point()
  
  
  group_by(penalty, .config)
  select(.config, params) %>%
  distinct %>%
  unnest(params)
  filter(.config == .config[train_test == 'tank']) %>%
  select(params) %>%
  unnest(params)

extract_fit_engine(fit_model) %>% 
  coef %>% 
  as.matrix %>% 
  as_tibble(rownames = 'param') %>%
  filter(param != '(Intercept)') %>%
  pivot_longer(cols = starts_with('s'),
               names_to = 'penalty_cat',
               values_to = 'coef') %>%
  mutate(lambda = str_extract(penalty_cat, '[0-9]+') %>% as.integer(),
         lambda = extract_fit_engine(fit_model)$lambda[lambda + 1]) %>%
  
  ggplot(aes(x = lambda, y = coef, group = param)) +
  geom_vline(xintercept = 1.06065935783328e-10) +
  geom_line()

