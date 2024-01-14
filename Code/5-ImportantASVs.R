# cpm ~ (control-healthyEnd, disease-healthyEnd, disease-diseaseEnd) * dosed/not 

#### Libraries ####
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
# library(rstanarm)
library(lubridate)
library(emmeans)
library(broom)
library(patchwork)
library(multidplyr)

#### N ASVs ####
alpha <- 0.05
rerun_models <- FALSE

#### Data ####
the_nmds <- read_rds('../intermediate_files/field_tank_nmds.rds.gz')
asv_fit <- read_rds('../intermediate_files/field_tank_asvArrows.rds.gz')

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

field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE)

tank_asv_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                     show_col_types = FALSE) %>%
  filter(geno != 'GE') %>%
  select(-health) %>%
  # group_by(sample_id, time, tank, anti_health, geno, resist) %>%
  separate(anti_health, into = c('with_anti', 'disease_exposure')) %>%
  rename(health = resist) %>%
  filter(with_anti == 'N') %>%
  mutate(time = str_extract(time, '[0-9]') %>% as.integer(),
         fragment_id = str_c(geno, disease_exposure, tank, sep = '_')) %>%
  mutate(exposure = as.integer(time == 0),
         exp_dis = as.integer(time != 0) * as.integer(disease_exposure == 'D'),
         exp_hea = as.integer(time != 0) * as.integer(disease_exposure == 'N'),
         tank_id = str_c(tank, disease_exposure, sep = '_')) %>%
  select(asv_id, sample_id, log2_cpm_norm, exposure, exp_dis, exp_hea, health,
         time, geno, fragment_id, tank_id) %>%
  mutate(dose = case_when(exp_dis == 1 ~ 'D',
                          exp_hea == 1 ~ 'H',
                          TRUE ~ NA_character_),
         time_c = time + 1) %>%
  nest(data = -c(asv_id)) 


read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
         show_col_types = FALSE) %>%
  filter(geno != 'GE') %>%
  select(-health) %>%
  # group_by(sample_id, time, tank, anti_health, geno, resist) %>%
  separate(anti_health, into = c('with_anti', 'disease_exposure')) %>%
  rename(health = resist) %>%
  filter(with_anti == 'N') %>%
  select(sample_id, tank, with_anti, disease_exposure, geno, health, time) %>%
  distinct %>%
  pivot_wider(names_from = 'time',
              values_from = sample_id) %>%
  count(geno, disease_exposure)

#### NMDS ####
species_fits <- scores(the_nmds)$species %>%
  as_tibble(rownames = 'asv_id') %>%
  left_join(tibble(asv_id = names(asv_fit$vectors$r),
                   r2 = asv_fit$vectors$r,
                   p = asv_fit$vectors$pvals),
            by = 'asv_id') %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(useful = str_c(family, genus, asv_id, sep = '_') %>%
           str_replace_all('__', '_')) %>% 
  filter(family %in% unique(shap_importance$family)) %>%
  mutate(important_asv = asv_id %in% unique(shap_importance$asv_id)) %>%
  identity()
  

species_fits %>%
  arrange(important_asv) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_segment(aes(colour = important_asv),
               xend = 0, yend = 0, show.legend = FALSE) +
  geom_point(aes(colour = important_asv, size = important_asv)) +
  geom_point(data = scores(the_nmds)$sites %>%
               as_tibble(rownames = 'sample_id') %>%
               left_join(select(field_data, sample_id, health) %>%
                           distinct(),
                         by = 'sample_id'),
             aes(fill = health),
             size = 2, shape = 21) +
  geom_text(data = . %>% filter(important_asv),
            aes(label = useful), show.legend = FALSE) +
  facet_wrap(~family) +
  scale_colour_manual(values = c('TRUE' = 'black', 'FALSE' = 'grey50')) +
  scale_size_manual(values = c('TRUE' = 5, 'FALSE' = 1)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = 'black')),
         size = 'none') +
  labs(colour = 'ASV Importance',
       fill = 'Health State') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'))
ggsave('../Results/nmds_important_asvs.png', height = 15, width = 15)

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
cluster_library(cluster, c('rstanarm', 'lmerTest', 'emmeans', 'dplyr', 
                           'tidyr', 'tibble', 'ggplot2', 'stringr',
                           'lubridate'))
cluster_copy(cluster, c('plot_field', 'extract_significance',
                        'make_coef_table', 'test_within_timepoint',
                        'make_fc_table'))

if(file.exists('../intermediate_files/field_asv_models.rds.gz') & !rerun_models){
  field_asv_models <- read_rds('../intermediate_files/field_asv_models.rds.gz')
} else {
  field_asv_models <- field_data %>%
    mutate(year = as.character(year)) %>%
    select(asv_id, sample_id, health, year, season, site, log2_cpm_norm) %>%
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
  # mutate(across(starts_with('p.'), ~p.adjust(., 'fdr'))) %>%
  # filter(if_any(c(starts_with('p.'), contains('health')), ~. < alpha)) %>%
  
  left_join(taxonomy, by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(plot + labs(title = str_c(order, family, genus, asv_id, sep = '\n')))) %>%
  ungroup %>%
  pull(plot) %>%
  wrap_plots(byrow = TRUE)
ggsave('../Results/important_asvs_field.png', height = 15, width = 15)


#### Model ASVs in tanks ####
# model<-tank_asv_models$model[[1]]; data <- tank_asv_models$data[[1]]; field_fitted <- tank_asv_models$just_health_fit[[1]]
# field_fitted <- field_asv_models$all_health_fit[[1]]

plot_tank <- function(model, data, field_fitted, caption = FALSE, alpha = 1){
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
        mutate(effect = str_replace_all(effect, c('health' = 'Resistance',
                                                  'exp_dis' = 'Disease Treatment',
                                                  'exp_hea' = 'Healthy Treatment'))) %>%
        filter(`Pr(>F)` < alpha) %>%
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
  
  out <- emmeans(model, ~health * (exp_dis + exp_hea),
          data = data) %>%
    as_tibble() %>%
    rename_with(~str_replace(., 'HPD', 'CL')) %>%
    filter(!(exp_dis == 1 & exp_hea == 1)) %>%
    mutate(grouping = str_c('g', row_number()),
           exposure = as.integer(!(exp_dis == 0 & exp_hea == 0)),
           dose = case_when(exp_dis == 1 ~ 'D',
                            exp_hea == 1 ~ 'H',
                            TRUE ~ NA_character_)) %>%
    rowwise %>%
    mutate(time_c = list(time_c = seq(1, 9, length.out = 100))) %>% 
    unnest(time_c) %>% 
    filter(!(exposure == 1 & time_c < 3 + 0.25),
           !(exposure == 0 & time_c > 3 - 0.25)) %>% 
    mutate(health = as.character(health)) %>%
    
    ggplot(aes(x = time_c - 1, y = emmean, colour = dose,
               group = grouping)) +
    geom_vline(data = tibble(health = c('R', 'S'), treat_tank = 1.75), 
               aes(xintercept = treat_tank)) +
    
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = dose),
                colour = NA, alpha = 0.5) +
    geom_line() +
    geom_pointrange(data = data %>%
                      group_by(time_c, health, dose) %>%
                      summarise(mean_log2 = mean(log2_cpm_norm),
                                se_log2 = sd(log2_cpm_norm) / sqrt(n()),
                                n = n(),
                                .groups = 'drop'),
                    aes(y = mean_log2,
                        ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2,
                        group = 1),
                    # size = 1.5, linewidth = 1,
                    position = position_dodge2(0.5)) +
    
    
    
    scale_x_continuous(breaks = c(0, 2, 8)) +
    facet_grid(~health, scales = 'free_x', space = 'free_x') +
    labs(caption = the_caption)
  
  if(!is.null(field_fitted)){
    out <- out +
      geom_pointrange(data = as_tibble(field_fitted) %>%
                        rename(dose = health) %>%
                        mutate(health = 'field'),
                      inherit.aes = FALSE, 
                      aes(x = 1, y = emmean, colour = dose,
                          ymin = lower.CL, ymax = upper.CL),
                      position = position_dodge(1.5))
  }
  out
}

plot_tank_foldChange <- function(model, data, field_fitted, caption = FALSE, alpha = 1){
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
        mutate(effect = str_replace_all(effect, c('health' = 'Resistance',
                                                  'exp_dis' = 'Disease Treatment',
                                                  'exp_hea' = 'Healthy Treatment'))) %>%
        filter(`Pr(>F)` < alpha) %>%
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

  out <- emmeans(model, ~health | (exp_dis + exp_hea),
          data = data) %>%
    contrast('revpairwise') %>% 
    broom::tidy(conf.int = TRUE) %>%
    filter(!(exp_dis == 1 & exp_hea == 1)) %>%
    mutate(time = if_else(exp_dis == 0 & exp_hea == 0, 'Before', 'After'),
           dose = case_when(exp_dis == 1 ~ 'Disease',
                            exp_hea == 1 ~ 'Healthy',
                            TRUE ~ NA_character_)) %>%
    bind_rows(contrast(field_fitted, 'pairwise') %>% 
                broom::tidy(conf.int = TRUE) %>%
                mutate(time = 'Field')) %>%
    mutate(time = factor(time, levels = c('Field', 'Before', 'After')),
           dose = if_else(is.na(dose), 'None', dose),
           group = str_c('g', row_number())) %>%
    ggplot(aes(x = time, y = estimate, colour = dose, group = group)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = 0.1, position = position_dodge(0.5))  +
    geom_point(position = position_dodge(0.5)) +
    labs(caption = the_caption)
  
  out
}

cluster_copy(cluster, c('plot_tank', 'plot_tank_foldChange'))

if(file.exists('../intermediate_files/tank_asv_models.rds.gz') & !rerun_models){
  tank_asv_models <- read_rds('../intermediate_files/tank_asv_models.rds.gz')
} else {
  tank_asv_models <- tank_asv_data %>%
    left_join(select(field_asv_models, asv_id, just_health_fit, all_health_fit),
              by = 'asv_id') %>%
    # filter(asv_id == 'ASV25') %>%
    arrange(asv_id) %>%
    rowwise %>%
    mutate(data = list(mutate(data, log2_cpm_norm = log2_cpm_norm + rnorm(nrow(data), 0, 6.11/100)))) %>%
    partition(cluster) %>%
    mutate(model = list(lmer(log2_cpm_norm ~ health * (exp_dis + exp_hea) + 
                               (1 | fragment_id) + (1 | tank_id) + (1 | geno),
                             data = data)),
           
           plot = list(plot_tank(model, data, field_fitted = just_health_fit, caption = TRUE) +
                         labs(y = 'log2(CPM)', x = 'Days') +
                         theme_classic() +
                         theme(legend.position = 'none',
                               panel.background = element_rect(colour = 'black'),
                               axis.text = element_text(colour = 'black'),
                               strip.background = element_blank())),
           
           plot_fc = list(plot_tank_foldChange(model, data, field_fitted = all_health_fit, caption = TRUE) +
                            labs(y = 'D log2(CPM) - H log2(CPM)', x = NULL) +
                            theme_classic() +
                            theme(legend.position = 'none',
                                  panel.background = element_rect(colour = 'black'),
                                  axis.text = element_text(colour = 'black'),
                                  strip.background = element_blank())),
           
           extract_significance(model),
           make_coef_table(~health * (exp_dis + exp_hea), model, data),
           make_fc_table(~health | (exp_dis + exp_hea), model, data)) %>%
    collect() %>%
    ungroup
  
  write_rds(tank_asv_models, '../intermediate_files/tank_asv_models.rds.gz')
}


tank_asv_models %>%
  
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>%
  
  # inner_join(filter(field_asv_models,
  #                   if_all(starts_with('p.withinYear'), ~. < alpha),
  #                   asv_id %in% shap_importance$asv_id) %>%
  #              select(asv_id),
  #            by = 'asv_id') %>%
  
  # select(-contains('Intercept')) %>%
  # mutate(across(starts_with('p.'), ~p.adjust(., 'fdr'))) %>%
  # filter(if_any(c(starts_with('p.') & contains('SX')), ~. < alpha)) %>%
  # filter(if_any(c(starts_with('p.') & contains('SX')))) %>%
  # filter(p.healthSXexpdis < alpha) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(plot + labs(title = str_c(order, family, genus, asv_id, sep = '\n')))) %>%
  ungroup %>%
  pull(plot) %>%
  wrap_plots()
# ggsave('../Results/important_asvs_tank.png', height = 15, width = 15)


tank_asv_models %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>%
  
  # inner_join(filter(field_asv_models,
  #                   if_all(starts_with('p.withinYear'), ~. < alpha),
  #                   asv_id %in% shap_importance$asv_id) %>%
  #              select(asv_id),
  #            by = 'asv_id') %>%
  # select(-contains('Intercept')) %>%
  # mutate(across(starts_with('p.'), ~p.adjust(., 'fdr'))) %>%
  # filter(if_any(c(starts_with('p.') & contains('SX')), ~. < alpha)) %>%
  # filter(if_any(c(starts_with('p.') & contains('SX')))) %>%
  # filter(p.healthSXexpdis < alpha) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  rowwise %>%
  mutate(plot_fc = list(plot_fc + labs(title = str_c(order, family, genus, asv_id, sep = '\n')))) %>%
  ungroup %>%
  pull(plot_fc) %>%
  wrap_plots()
ggsave('../Results/Fig5_important_asvs_FC.png', height = 15, width = 15)

#### Output asv coefficients ####
asv_coef <- full_join(select(field_asv_models, asv_id, starts_with('mu'), starts_with('fcDH')),
                      select(tank_asv_models, asv_id, starts_with('mu'), starts_with('fcDH')),
                      by = 'asv_id') %>%
  select(asv_id, starts_with('fcDH'), starts_with('mu_'))

write_csv(asv_coef, '../Results/asv_coefficients.csv.gz') 
