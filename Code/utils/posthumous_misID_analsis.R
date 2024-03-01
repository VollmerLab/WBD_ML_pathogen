
library(tidyverse)
library(tidymodels)
library(magrittr)

coral_split <- read_rds('../../intermediate_files/coral_split.rds.gz')

field_asvs <- bind_rows(train = analysis(coral_split),
                        test = assessment(coral_split),
                        .id = 'train_test')

field_data <- read_csv('../../intermediate_files/normalized_field_asv_counts.csv', show_col_types = FALSE) %>%
  select(sample_id, health, year, season, site, dataset) %>%
  distinct %>%
  left_join(field_asvs %>%
              select(train_test, sample_id),
            by = 'sample_id')

model_predictions <- list.files('../../Results/model_tuning/metrics/',
                                pattern = 'prediction.csv.gz$', full.names = TRUE) %>%
  tibble(pred_file = .) %>%
  mutate(wflow_id = str_extract(pred_file, '[a-zA-Z0-9_]+.csv.gz$') %>% 
           str_remove('_prediction.csv.gz$')) %>%
  rowwise(wflow_id) %>%
  reframe(read_csv(pred_file, show_col_types = FALSE)) %>%
  select(wflow_id, sample_id, .pred_class)

model_misid_rates <- inner_join(field_data,
           model_predictions,
           by = 'sample_id') %>%
  mutate(correct_id = health == .pred_class) %>%
  group_by(wflow_id) %>%
  summarise(accuracy = mean(correct_id))

base_rate <- count(field_data, health) %>%
  mutate(prop = n/sum(n)) %>%
  filter(health == 'H') %>%
  pull(prop)

mis_id <- inner_join(field_data,
           model_predictions,
           by = 'sample_id') %>%
  mutate(correct_id = health == .pred_class) %>%
  select(-.pred_class) %>%
  pivot_wider(names_from = wflow_id,
              values_from = correct_id) %>%
  select(-base_null, -base_knn, -base_pls) %>%
  filter(if_any(starts_with('base'), ~!.)) %>%
  rowwise() %>%
  mutate(n_model = sum(!c_across(starts_with('base')))) %>%
  ungroup

top_asvs <- read_csv('../../Results/Table45_asv_table.csv', 
                     show_col_types = FALSE) %>%
  filter(likely_type != '-') %>%
  select(ID, likely_type)

#### Basic Stats ####
mis_id %>%
  count(n_model)

mis_id %>%
  filter(n_model > 1)

#### Breakdown who is mis ID'd ####
full_join(count(mis_id, health),
          count(field_data, health, name = 'total'),
          by = 'health') %>%
  mutate(not = total - n) %>%
  select(-total) %>%
  column_to_rownames('health') %T>%
  print %>%
  chisq.test()

prop.test(c(7,5), n = c(143,270))

full_join(count(mis_id, year),
          count(field_data, year, name = 'total'),
          by = 'year') %>%
  mutate(not = total - n) %>%
  select(-total) %>%
  column_to_rownames('year') %T>%
  print %>%
  chisq.test()

full_join(count(mis_id, season),
          count(field_data, season, name = 'total'),
          by = 'season') %>%
  mutate(not = total - n) %>%
  select(-total) %>%
  column_to_rownames('season') %T>%
  print %>%
  chisq.test()

full_join(count(mis_id, year, season),
          count(field_data, year, season, name = 'total'),
          by = c('year', 'season')) %>%
  mutate(timepoint = str_c(year, season, sep = '_'),
         across(where(is.integer), ~replace_na(., 0L)),
         .keep = 'unused') %>%
  mutate(not = total - n) %>%
  select(-total) %>%
  column_to_rownames('timepoint') %T>%
  print %>%
  # chisq.test()
  chisq.posthoc.test::chisq.posthoc.test()

full_join(count(mis_id, site),
          count(field_data, site, name = 'total'),
          by = 'site') %>%
  mutate(not = total - n) %>%
  select(-total) %>%
  column_to_rownames('site') %T>%
  print %>%
  chisq.test()

full_join(count(mis_id, train_test),
          count(field_data, train_test, name = 'total'),
          by = 'train_test') %>%
  mutate(not = total - n) %>%
  select(-total) %>%
  column_to_rownames('train_test') %T>%
  print %>%
  chisq.test()

#### Breakdown of Important ASVs and misID'd ####

select(mis_id, sample_id) %>%
  left_join(field_asvs, 
            by = 'sample_id') %>%
  select(sample_id, health, all_of(top_asvs$ID)) %>%
  select(-ASV38) %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = 'asv_id',
               values_to = 'log2cpm') %>%
  group_by(asv_id) %>%
  summarise(broom::tidy(t.test(log2cpm ~ health)))

tst <- select(mis_id, sample_id) %>%
  left_join(field_asvs, 
            by = 'sample_id') %>%
  select(sample_id, health, all_of(top_asvs$ID)) %>%
  select(-ASV38) %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = 'asv_id',
               values_to = 'log2cpm') %>%
  group_by(asv_id) %>%
  summarise(t_test = list(MKinfer::boot.t.test(log2cpm ~ health)))

select(mis_id, sample_id) %>%
  left_join(field_asvs, 
            by = 'sample_id') %>%
  select(sample_id, health, all_of(top_asvs$ID)) %>%
  select(-ASV38) %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = 'asv_id',
               values_to = 'log2cpm') %>%
  ggplot(aes(x = asv_id, y = log2cpm, colour = health)) +
  stat_summary(position = position_dodge(0.5))


select(mis_id, sample_id) %>%
  left_join(mutate(field_asvs, across(starts_with('ASV'), ~. > min(.))), 
            by = 'sample_id') %>%
  select(sample_id, health, all_of(top_asvs$ID)) %>%
  group_by(health) %>%
  summarise(across(starts_with('ASV'), sum),
            n = n()) %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = 'asv_id',
               values_to = 'present') %>%
  group_by(asv_id) %>%
  summarise(present = list(present),
            total = list(n),
            .groups = 'rowwise') %>%
  reframe(broom::tidy(prop.test(present, total)))

  
#### SHAP FORCE PLOTS ####
library(shapviz)
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

taxonomy <- read_csv('../../intermediate_files/taxonomy.csv.gz', show_col_types = FALSE)
top_shap <- read_csv('../../intermediate_files/model_shaps.csv.gz', show_col_types = FALSE) 


top_shap %>%
  filter(sample_id == '2016_W_SB_H_15',
         wflow_id == 'base_mlp') %>%
  make_shap_force_plot(model_name = 'base_mlp', 
                       base_value = base_rate, 
                       collapse_to = 'genus')


tst <- mis_id %>%
  filter(health == 'H') %>%
  select(sample_id, starts_with('base')) %>%
  pivot_longer(cols = starts_with('base'),
               names_to = 'model',
               values_to = 'correct_id') %>%
  filter(!correct_id) %>%
  select(-correct_id) %>%
  nest_by(model) %>%
  mutate(data = list(pull(data, sample_id))) %>%
  rename(samples = data) %>%
  left_join(model_misid_rates,
            by = c('model' = 'wflow_id')) %>%
  mutate(data = list(filter(top_shap, sample_id %in% samples))) %>%
  mutate(plot = list(make_shap_force_plot(data,
                                          model_name = model, 
                                          base_value = base_rate, 
                                          collapse_to = 'genus')))


mis_id %>%
  filter(health == 'H') %>%
  select(sample_id, starts_with('base')) %>%
  pivot_longer(cols = starts_with('base'),
               names_to = 'wflow_id',
               values_to = 'correct_id') %>%
  
  left_join(top_shap,
            by = c('sample_id', 'wflow_id')) %>%
  select(sample_id, correct_id, wflow_id, starts_with('ASV8_'), starts_with('ASV25_'))
