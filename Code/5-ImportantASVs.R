# cpm ~ (control-healthyEnd, disease-healthyEnd, disease-diseaseEnd) * dosed/not 

#### Libraries ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(lubridate)
library(emmeans)
library(broom)
library(patchwork)
library(multidplyr)
library(ggtext)

#### N ASVs ####
alpha <- 0.05
rerun_models <- FALSE

#### Data ####
model_list <- read_csv('../Results/equivilant_top_models.csv.gz',
                       show_col_types = FALSE) %>%
  filter(pract_equiv >= 0.95) %>%
  rename(wflow_id = model) %>%
  pull(wflow_id)

taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz',
                     show_col_types = FALSE) %>%
  mutate(across(where(is.character), 
                ~str_replace_na(., replacement = '')))


shap_importance <- read_csv('../Results/asv_importance.csv.gz',
                            show_col_types = FALSE) %>%
  filter(p_adjust < alpha) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(asv_id = fct_reorder(asv_id, median_rank)) %>%
  select(asv_id, domain:genus) %>%
  distinct %>%
  arrange(asv_id)

asv_associations <- read_csv('../intermediate_files/model_shaps.csv.gz', show_col_types = FALSE) %>%
  filter(wflow_id %in% str_c('base_', c('mlp', 'svmLinear', 'lasso', 'forest'))) %>%
  pivot_longer(starts_with('ASV'),
               names_to = c('asv_id', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  group_by(sample_id, asv_id, health, value) %>%
  summarise(shap = mean(shap),
            .groups = 'drop') %>%
  group_by(asv_id) %>%
  reframe(broom::tidy(lm(shap ~ value))) %>%
  filter(term != '(Intercept)') %>%
  arrange(desc(asv_id)) %>%
  mutate(`Disease Associated` = estimate > 0,
         `Healthy Associated` = estimate < 0) %>%
  select(asv_id, `Disease Associated`, `Healthy Associated`)

field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE)

whacky <- c('Bin1_N_D_YE', 'Bin2_N_D_YE', 'Bin5_N_D_OR', 'Bin5_N_D_GR')
tank_asv_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                          show_col_types = FALSE) %>% 
  # filter(asv_id == 'ASV25') %>%
  filter(anti == 'N',
         sample_id != 'P6_Bin2_N_N_GR') %>% 
  mutate(time = str_extract(time, '[0-9]') %>% as.integer(),
         fragment_id = str_c(geno, exposure, str_extract(tank, 'Bin[0-9]+'), sep = '_'),
         tank_id = str_c(exposure, str_extract(tank, 'Bin[0-9]+'), sep = '_'),
         treatment = if_else(exposure == 'pre', 'pre',
                             str_c(exposure, health, sep = '_')),
         time_treat = str_c(time, exposure, health, sep = '_')) %>%
  select(-dataset, -anti, -plate, -year:-site, -cpm_norm:-cpm, 
         -norm.factors, -domain:-species) %>%
  
  #Get rid of samples which shouldnt exist
  filter(!(str_detect(sample_id, str_c(whacky, collapse = '|')) & time == 8))


tank_asv_data %>%
  count(tank_id, time) %>%
  pivot_wider(names_from = time, values_from = n)

tank_asv_data %>%
  count(geno, time) %>%
  pivot_wider(names_from = time, values_from = n)


#### Model ASVs in Field ####
# model <- field_asv_models$model[[1]]
# data <- field_asv_models$data[[1]]

plot_field <- function(model, data, caption = FALSE, alpha = 0.05, alpha_filter = 1){
  if(all(class(model) == 'lm')){
    plot_data <- emmeans(model, ~health * year * season * site,
                         data = data) %>%
      as_tibble()
  }
  
  if(all(class(model) == 'lmerModLmerTest') | 'stanreg' %in% class(model)){
    plot_data <- emmeans(model, ~health * year * season,
                         data = data) %>%
      as_tibble()
  }
  
  the_caption <- NULL
  if(caption){
    if('stanreg' %in% class(model)){
      the_caption <- NULL
    } else {
      if(all(class(model) == 'lm')){
        the_aov <- anova(model)
        
        the_aov <- as_tibble(the_aov, rownames = 'effect') %>%
          filter(effect != 'Residuals') %>%
          rename(NumDF = Df) %>%
          mutate(DenDF = the_aov$Df[length(the_aov$Df)])
        
      } else {
        the_aov <- anova(model, ddf = 'Kenward-Roger') %>%
          as_tibble(rownames = 'effect') 
      }
      the_caption <- the_aov %>%
        filter(`Pr(>F)` < alpha_filter) %>%
        mutate(effect = str_replace_all(effect, c('_c' = '', '_' = ' ', ':' = ' x ')),
               effect = str_to_title(effect) %>% str_replace('X', 'x'),
               across(c(where(is.numeric), -`F value`, -`Pr(>F)`), ~round(., 1)),
               `F value` = round(`F value`, 3),
               output = str_c(effect, ': F(', NumDF, ', ', DenDF, ') = ', 
                              `F value`, '; p = ', scales::pvalue(`Pr(>F)`))) %>%
        pull(output) %>%
        str_c(collapse = '\n')
    }
  }
  
  out <- plot_data %>%
    rename_with(~str_replace(., 'HPD', 'CL')) %>%
    mutate(grouping = row_number()) %>%
    mutate(timepoint = if_else(season == 'S', 
                               ymd(str_c(year, 'July', 1, sep = '-')),
                               ymd(str_c(year, 'January', 1, sep = '-')))) %>%
    ggplot(aes(x = timepoint, y = emmean, colour = health,
               ymin = lower.CL, ymax = upper.CL, group = grouping)) +
    geom_pointrange(position = position_dodge(50)) +
    labs(caption = the_caption)
  
  if(all(class(model) == 'lmerModLmerTest')){
    timepoint_diff <- emmeans(model, ~health | year + season) %>%
      contrast('pairwise') %>%
      as_tibble() %>%
      mutate(grouping = row_number(),
             timepoint = if_else(season == 'S', 
                                 ymd(str_c(year, 'July', 1, sep = '-')),
                                 ymd(str_c(year, 'January', 1, sep = '-')))) 
    
    out <- out +
      geom_text(data = timepoint_diff, inherit.aes = FALSE,
                aes(x = timepoint, y = Inf, 
                    label = if_else(p.value < alpha, '*', '')),
                vjust = 1, size = 8)
  }
  
  out
}

extract_significance <- function(model, anova_summary = 'summary'){
  if('lmerModLmerTest' %in% class(model)){
    if(anova_summary == 'summary'){
      anova_table <- summary(model, ddf = 'Kenward-Roger')$coefficients
    } else if(anova_summary == 'anova'){
      anova_table <- anova(model, ddf = 'Kenward-Roger')
    }
    
    out <- as_tibble(anova_table, rownames = 'param') %>%
      select(param, starts_with('Pr')) %>%
      pivot_wider(names_from = param, 
                  values_from = starts_with('Pr'),
                  names_prefix = 'p.') %>%
      rename_with(~str_replace_all(., c(':' = 'X', '_' = '')) %>%
                    str_remove_all('\\(|\\)'))
    
    out$summary_table <- list(anova_table)
  }
  
  if('stanreg' %in% class(model)){
    anova_table <- summary(model, 
                           pars = names(fixef(model)),
                           probs = c(0.1, 0.9)) %>%
      as_tibble(rownames = 'param') %>%
      mutate(sig = (`10%` < 0 & `90%` < 0) | (`10%` > 0 & `90%` > 0))
    
    out <- select(anova_table, param, sig) %>%
      pivot_wider(names_from = param, 
                  values_from = sig,
                  names_prefix = 'p.') %>%
      rename_with(~str_replace_all(., c(':' = 'X', '_' = '')) %>%
                    str_remove_all('\\(|\\)'))
    
    out$summary_table <- list(anova_table)
  }
  
  if(all(class(model) == 'lm')){
    if(anova_summary == 'summary'){
      anova_table <- summary(model)$coefficients
    } else if(anova_summary == 'anova'){
      anova_table <- anova(model)
    }
    
    out <- as_tibble(anova_table, rownames = 'param') %>%
      select(param, starts_with('Pr')) %>%
      pivot_wider(names_from = param, 
                  values_from = starts_with('Pr'),
                  names_prefix = 'p.') %>%
      rename_with(~str_replace_all(., c(':' = 'X', '_' = '')) %>%
                    str_remove_all('\\(|\\)'))
    
    out$summary_table <- list(anova_table)
  }
  
  out
}


# model <- field_asv_models$model[[1]]; data <- field_asv_models$data[[1]]; form <- ~health | season * year
# model <- tank_asv_models$model[[1]]; data <- tank_asv_models$data[[1]]; form <- ~health | exp_dis + exp_hea

make_coef_table <- function(form, model, data){
  
  x_vars <- str_split(as.character(form)[2], ' [\\+\\*] ') %>% unlist %>% unique %>%
    str_remove_all('\\(|\\)')
  
  out <- emmeans(model, form, data = data) %>%
    as_tibble %>%
    select(all_of(x_vars), emmean)
  
  if(any(x_vars == 'exp_dis')){
    out <- out %>%
      filter(!(exp_dis == 1 & exp_hea == 1)) %>%
      mutate(treatment = case_when(exp_dis == 0 & exp_hea == 0 ~ str_c('pre', health, sep = '_'),
                                   exp_dis == 1 ~ str_c('postDisease', health, sep = '_'),
                                   exp_hea == 1 ~ str_c('postHealth', health, sep = '_'),)) %>%
      select(treatment, emmean)
    x_vars <- 'treatment'
  } 
  
  out %>%
    pivot_wider(names_from = all_of(x_vars),
                values_from = emmean,
                names_prefix = 'mu_')
}

make_fc_table <- function(form, model, data){
  x_vars <- str_split(as.character(form)[2], ' [\\+\\*\\|] ') %>% unlist %>% unique %>%
    str_remove_all('\\(|\\)')
  
  out <- emmeans(model, form, data = data) %>%
    contrast(if_else(any(x_vars == 'exp_dis'), 'revpairwise', 'pairwise')) %>%
    as_tibble %>%
    select(any_of(x_vars), estimate)
  
  if(any(x_vars == 'exp_dis')){
    out <- out %>%
      filter(!(exp_dis == 1 & exp_hea == 1)) %>%
      mutate(treatment = case_when(exp_dis == 0 & exp_hea == 0 ~ 'pre',
                                   exp_dis == 1 ~ 'postDisease',
                                   exp_hea == 1 ~ 'postHealth')) %>%
      select(treatment, estimate)
    x_vars <- 'treatment'
  } 
  
  out %>%
    pivot_wider(names_from = any_of(x_vars),
                values_from = estimate,
                names_prefix = 'fcDH_')
}

test_within_timepoint <- function(model, data){
  tp_contrasts <- emmeans(model, ~health | year * season, data = data) %>%
    contrast('pairwise') %>%
    as_tibble
  
  select(tp_contrasts, year, season, p.value) %>%
    pivot_wider(names_from = c(year, season),
                values_from = p.value,
                names_prefix = 'p.withinYear_') %>%
    mutate(tp_contrasts = list(tp_contrasts))
  
}

cluster <- new_cluster(parallel::detectCores() / 2)
cluster_library(cluster, c('lmerTest', 'emmeans', 'dplyr', 
                           'tidyr', 'tibble', 'ggplot2', 'stringr',
                           'lubridate', 'magrittr'))
cluster_copy(cluster, c('plot_field', 'extract_significance',
                        'make_coef_table', 'test_within_timepoint',
                        'make_fc_table'))

if(file.exists('../intermediate_files/field_asv_models.rds.gz') & !rerun_models){
  field_asv_models <- read_rds('../intermediate_files/field_asv_models.rds.gz')
} else {
  field_asv_models <- field_data %>%
    mutate(year = as.character(year)) %>%
    select(asv_id, sample_id, health, year, season, site, log2_cpm_norm, n_reads, lib.size) %>%
    nest(data = -c(asv_id)) %>%
    # filter(asv_id == 'ASV25') %>%
    rowwise %>%
    mutate(data = list(mutate(data, log2_cpm_norm = log2_cpm_norm + rnorm(nrow(data), 0, 6.11/100)))) %>%
    partition(cluster) %>%
    
    mutate(model = list(lmer(log2_cpm_norm ~ health * year * season + (1 | site),
                           data = data)),
         plot = list(plot_field(model, data, caption = TRUE) +
                       labs(y = 'log2(CPM)', x = NULL) +
                       theme_classic() +
                       theme(legend.position = 'none',
                             panel.background = element_rect(colour = 'black'),
                             axis.text = element_text(colour = 'black'))),
         extract_significance(model, anova_summary = 'anova'),
         
         average_health_fit = list(emmeans(model, ~health)),
         just_health_fit = list(emmeans(model, ~health, at = list(year = '2017', season = 'W'))),
         all_health_fit = list(emmeans(model, ~health | year * season)),
         
         test_within_timepoint(model, data),
         make_coef_table(~health * year * season, model, data),
         make_fc_table(~health | year * season, model, data)) %>%
    
    
    collect %>%
    ungroup 
  write_rds(field_asv_models, '../intermediate_files/field_asv_models.rds.gz')
}

field_asv_models %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>%
  
  left_join(taxonomy, by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(plot + labs(title = str_c(order, family, genus, asv_id, sep = '\n'),
                                 caption = NULL) +
                       theme(plot.title = element_text(size = 10)))) %>%
  ungroup %>%
  pull(plot) %>%
  wrap_plots(byrow = TRUE)

field_asv_models %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(across(contains('p.within'), ~p.adjust(., method = 'fdr'))) %>%
  
  rowwise %>%
  mutate(consistent = all(c_across(starts_with('p.within')) < alpha) & 
           (all(c_across(starts_with('fcDH')) > 0) | 
              all(c_across(starts_with('fcDH')) < 0)),
         .keep = 'unused') %>%
  ungroup %>%
  select(asv_id, consistent) %>%
  left_join(asv_associations, by = 'asv_id') %>%
  write_csv('../intermediate_files/ML_asv_groups.csv')


#### Model ASVs in tanks ####
# model<-tank_asv_models$model[[122]]; data <- tank_asv_models$data[[122]]; field_fitted <- tank_asv_models$all_health_fit[[122]]

plot_prePost <- function(model, data, field_fitted){
  real_zero <- min(data$log2_cpm_norm)
  initial_emmean <- emmeans(model, ~time_treat)
  treatment_levels <- levels(initial_emmean)$time_treat
  
  
  ref_grid(model) %>%
    add_grouping('treat_out', 'time_treat', 
                 str_remove(treatment_levels, '[0-9]+_')) %>%
    emmeans(~ treat_out) %>%
    broom::tidy(conf.int = TRUE) %>%
    mutate(prePost = if_else(str_detect(treat_out, 'pre'), 'Pre', 'Post'),
           prePost = factor(prePost, levels = c('Field', 'Pre', 'Post'), labels = c('Field', 'Pre-Exposure', 'Post-Exposure')),
           
           exposure = str_extract(treat_out, '^[DN]'),
           health = str_extract(treat_out, '[DH]$')) %>%
    bind_rows(broom::tidy(field_fitted, conf.int = TRUE) %>%
                mutate(prePost = factor('Field', levels = c('Field', 'Pre', 'Post')))) %>%
    mutate(health = factor(health, levels = c('H', 'D'), labels = c('Healthy', 'Diseased')),
           exposure = if_else(is.na(exposure), 'Field', exposure),
           exposure = factor(exposure, levels = c('Field', 'N', 'D'), labels = c('Field', 'Control', 'Disease')),
           season = case_when(is.na(season) ~ NA_character_,
                              season == 'S' ~ 'Jul',
                              season == 'W' ~ 'Jan'),
           timepoint = if_else(is.na(year), 'Tank', str_c(year, season, sep = '_'))) %>%
    
    # mutate(across(c(estimate, conf.low, conf.high), ~. - real_zero)) %>%
    
    
    ggplot(aes(x = prePost, y = estimate, colour = health, shape = exposure)) +
    # geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
    #                 position = position_dodge2(0.5)) +
    geom_linerange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge2(0.5), 
                  linewidth = 1, show.legend = FALSE) +
    geom_point(position = position_dodge2(0.5), size = 5) +
    # scale_shape_manual(values = c('D' = 'circle filled', 'N' = 'triangle filled'), 
    #                    na.value = 'square filled') +
    scale_shape_manual(values = c('Disease' = 'diamond', 'Control' = 'square', 'Field' = 'circle')) +
    scale_colour_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, type = "continuous"),
                                           c('Healthy', 'Diseased')),
                        breaks = c('Diseased', 'Healthy'))
  
}

tank_posthocs <- function(model){
  #Significance between healthy vs disease after exposure
  #significance of diseased compared to before exposure
  #significance of healthy compared to before exposure
  initial_emmean <- emmeans(model, ~time_treat,  
                            lmer.df = 'kenward-roger')
  treatment_levels <- levels(initial_emmean)$time_treat
  
  all_contrasts <- list('DvH' = case_when(str_detect(treatment_levels, '^.*_D$') ~ 1/2, #diseased vs healthy
                                          str_detect(treatment_levels, '^[28].*_H$') ~ -1/4,
                                          TRUE ~ 0),
                        
                        'DvN' = case_when(str_detect(treatment_levels, '^.*_D.*$') ~ 1/4, #diseased vs healthy exposure
                                          str_detect(treatment_levels, '^.*_N.*$') ~ -1/2,
                                          TRUE ~ 0),
                        
                        'DDvDH' = case_when(str_detect(treatment_levels, '^.*_D_D$') ~ 1/2,
                                            str_detect(treatment_levels, '^.*_D_H$') ~ -1/2,
                                            TRUE ~ 0),
                        
                        'DDvNH' = case_when(str_detect(treatment_levels, '^.*_D_D$') ~ 1/2,
                                            str_detect(treatment_levels, '^.*_N_H$') ~ -1/2,
                                            TRUE ~ 0),
                        
                        'DHvNH' = case_when(str_detect(treatment_levels, '^.*_D_H$') ~ 1/2,
                                            str_detect(treatment_levels, '^.*_N_H$') ~ -1/2,
                                            TRUE ~ 0),
                        
                        'PostvPreD' = case_when(str_detect(treatment_levels, '0_pre_H') ~ -1,
                                                str_detect(treatment_levels, '^.*_D$') ~ 1/2,
                                                TRUE ~ 0),
                        
                        'PostvPreH' = case_when(str_detect(treatment_levels, '0_pre_H') ~ -1,
                                                str_detect(treatment_levels, '^.*_H$') ~ 1/4,
                                                TRUE ~ 0))
  
  posthoc_sided <- contrast(initial_emmean, all_contrasts, side = '>', adjust = 'none')
  base_posthoc <- contrast(initial_emmean, all_contrasts, side = '=', adjust = 'none')
  
  full_posthoc <- base_posthoc %>%
    broom::tidy(conf.int = TRUE) %>%
    select(contrast, estimate, std.error, conf.low, conf.high) %>%
    left_join(broom::tidy(posthoc_sided) %>%
                select(contrast, df, statistic, p.value),
              by = 'contrast')
  
  select(full_posthoc, contrast, estimate, p.value) %>%
    rename(fc = estimate, pvalue = p.value) %>%
    pivot_wider(names_from = contrast, values_from = c(fc, pvalue)) %>%
    mutate(tank_posthoc = list(full_posthoc),
           tank_emmeans = list(base_posthoc),
           tank_directional_emmeans = list(posthoc_sided))
}

cluster_copy(cluster, c('plot_prePost', 'tank_posthocs'))

if(file.exists('../intermediate_files/tank_asv_models.rds.gz') & !rerun_models){
  tank_asv_models <- read_rds('../intermediate_files/tank_asv_models.rds.gz')
} else {
  tank_asv_models <- tank_asv_data %>%
    nest_by(asv_id) %>%
    left_join(select(field_asv_models, asv_id, all_health_fit, starts_with('p.withinYear')),
              by = 'asv_id') %>%
    arrange(asv_id) %>%
    rowwise %>%
    
    #Add a tiny bit of noise to help with fitting
    mutate(data = list(mutate(data, log2_cpm_norm = jitter(log2_cpm_norm)))) %>%
    partition(cluster) %>%
    mutate(model = list(lmer(log2_cpm_norm ~ time_treat + (1 | geno / fragment_id) + (1 | tank_id),
                             data = data))) %>%
    
    mutate(plot = list(plot_prePost(model, data, field_fitted = all_health_fit)),
           
           extract_significance(model, anova_summary = 'anova'),
           tank_posthocs(model)) %>%
    collect() %>%
    ungroup
  
  write_rds(tank_asv_models, '../intermediate_files/tank_asv_models.rds.gz')
}


consistent_asvs <- tank_asv_models %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(across(contains('p.within'), ~p.adjust(., method = 'fdr'))) %>%
  filter_at(vars(contains('p.within')), all_vars(. < alpha)) 

classified_asvs <- consistent_asvs %>%
  mutate(across(contains('pvalue_'), ~p.adjust(., method = 'fdr'))) %>%
  
  mutate(likely_type = case_when(pvalue_DDvDH < alpha & 
                                   #pvalue_DDvNH < alpha & 
                                   pvalue_DvH < alpha &
                                   pvalue_PostvPreD < alpha ~ 'Pathogen',
                                 
                                 pvalue_DvN < alpha &
                                   pvalue_PostvPreD < alpha ~ 'Opportunist',
                                 
                                 TRUE ~ 'Commensalist'),
         likely_type = factor(likely_type, levels = c('Pathogen', 'Opportunist', 'Commensalist')),
         asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>% 
  arrange(asv_id) 

classified_plots <- classified_asvs %>%
  filter(likely_type %in% c('Pathogen', 'Opportunist') | asv_id == 'ASV40') %>%
  
  left_join(taxonomy, by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = if_else(likely_type == 'Opportunist',
                        list(plot + scale_x_discrete(labels = ~str_replace(., '-', '-\n'))),
                        list(plot + scale_x_discrete(labels = ~str_replace(., '-', '-\n'))))) %>%
  mutate(plot = list(plot + 
                       geom_vline(xintercept = 1.5, linewidth = 0.5) +
                       scale_y_continuous(limits = c(4.75, 17.5),
                                          breaks = c(5, 7.5, 10, 12.5, 15, 17.5)) +
                       labs(y = 'log<sub>2</sub>(CPM)', 
                            x = NULL,
                            colour = 'Disease\nState',
                            shape = 'Experimental\nCondition') +
                       theme_classic() +
                       theme(panel.background = element_rect(colour = 'black'),
                             axis.text = element_text(colour = 'black', size = 12),
                             legend.text = element_text(colour = 'black', size = 12),
                             legend.key = element_blank(),
                             axis.title.y = element_markdown(colour = 'black', size = 18),
                             legend.title = element_text(colour = 'black', size = 14),
                             axis.text.x = element_text(colour = 'black', size = 10)))) %>%
  
  mutate(name = case_when(asv_id == 'ASV25' ~ '<i>Cysteiniphilum litorale</i>',
                          asv_id == 'ASV8' ~ '<i>Vibrio sp.</i>',
                          asv_id == 'ASV26' ~ 'Oceanospirillaceae',
                          asv_id == 'ASV30' ~ '<i>Thalassotalea sp.</i>',
                          asv_id == 'ASV361' ~ '<i>Endozoicomonas atrinae</i>',
                          asv_id == 'ASV38' ~ '<i>Neptuniibacter sp.</i>',
                          asv_id == 'ASV51' ~ '<i>Neptuniibacter sp.</i>',
                          asv_id == 'ASV40' ~ '<i>Qipengyuania sp.</i>',
                          TRUE ~ 'Unknown'),
         plot = list(plot + 
                       labs(title = str_c(asv_id, name, sep = ' - ')) +
                       theme(plot.title = element_markdown(colour = 'black', size = 10)))) %>%
  
  ungroup %>%
  filter(asv_id != 'ASV40') %>%
  arrange(likely_type) %>%
  mutate(asv_id = factor(asv_id, levels = str_c('ASV', c(25, 8, 30, 26, 361, 38, 51)))) %>%
  arrange(asv_id) %>%
  group_by(likely_type) %>%
  summarise(plot_set = list(wrap_plots(plot, nrow = 1) + plot_layout(axis_titles = 'collect_y', axes = 'collect_y'))) %>%
  pull(plot_set) %>%
  wrap_plots(nrow = 2) +
  plot_layout(guides = 'collect') +
  
  plot_annotation(tag_levels = list(c('A', '', '', 'B', '', '', '', 'C'))) &
  theme(plot.tag = element_text(face = 'bold', size = 22)) 
  
classified_plots
ggsave('../Results/Fig6_important_asvs.png', height = 12, width = 12)
ggsave('../Results/Fig6.svg', height = 12, width = 12)
ggsave('../Results/Fig6.pdf', height = 12, width = 12)
