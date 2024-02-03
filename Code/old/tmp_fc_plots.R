model <- filter(tank_asv_models, asv_id == 'ASV25') %>%
  select(model) %>%
  pull(model) %>%
  pluck(1)

data <- filter(tank_asv_models, asv_id == 'ASV25') %>%
  select(data) %>%
  pull(data) %>%
  pluck(1)

field_fitted <- filter(tank_asv_models, asv_id == 'ASV25') %>%
  select(all_health_fit) %>%
  pull(all_health_fit) %>%
  pluck(1)


data %>%
  filter(health == 'D')


tank_posthocs <- function(model){
  #Significance between healthy vs disease after exposure
  #significance of diseased compared to before exposure
  #significance of healthy compared to before exposure
  initial_emmean <- emmeans(model, ~time_treat,  
                            lmer.df = 'kenward-roger')
  treatment_levels <- levels(initial_emmean)$time_treat
  
  all_contrasts <- list('DvH' = case_when(str_detect(treatment_levels, '^.*_D$') ~ 1/2,
                                          str_detect(treatment_levels, '^.*_H$') ~ -1/4,
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
  
 posthoc_sided <- contrast(initial_emmean, all_contrasts, side = '>', adjust = 'fdr')
 full_posthoc <- contrast(initial_emmean, all_contrasts, side = '=', adjust = 'fdr') %>%
   broom::tidy(conf.int = TRUE) %>%
   select(contrast, estimate, std.error, conf.low, conf.high) %>%
   left_join(broom::tidy(posthoc_sided) %>%
               select(contrast, df, statistic, adj.p.value),
             by = 'contrast')
 
 select(full_posthoc, contrast, estimate, adj.p.value) %>%
   rename(fc = estimate, p = adj.p.value) %>%
   pivot_wider(names_from = contrast, values_from = c(fc, p)) %>%
   mutate(posthoc = list(full_posthoc))
}


# Pathogen test - 
#high in disease compared to healthy
#higher after exposure than before - only in diseased

# in disease exposure up in diseased but not healthy






#opportunist is up in disease exposure 


plot_tank <- function(model, data, field_fitted){
  
  out <- emmeans(model, ~ time_treat) %>%
    as_tibble() %>%
    separate(time_treat, into = c('time', 'exposure', 'health'), 
             sep = '_', fill = 'right') %>%
    mutate(time = as.integer(time)) %>%
    rename_with(~str_replace(., 'HPD', 'CL')) %>%
    
    ggplot(aes(x = time, y = emmean, colour = health,
               shape = exposure)) +
    geom_vline(xintercept = 1.5, linetype = 'dashed', linewidth = 0.1) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  position = position_dodge(1), width = 0.1, linewidth = 1.5) +
    geom_point(position = position_dodge(1), size = 3) +
    # geom_point(data = data, aes(y = log2_cpm_norm),
    #            # size = 1.5, linewidth = 1,
    #            position = position_dodge2(0.5)) +
    scale_x_continuous(breaks = c(0, 2, 8))
  
  if(!is.null(field_fitted)){
    out <- out +
      geom_vline(xintercept = -0.5, linetype = 'solid', linewidth = 0.1) +
      geom_pointrange(data = as_tibble(field_fitted) %>%
                        mutate(exposure = 'field'),
                      inherit.aes = FALSE, size = 0.5,
                      aes(x = -1.5, y = emmean, colour = health,
                          shape = exposure,
                          ymin = lower.CL, ymax = upper.CL),
                      position = position_dodge2(1))
  }
  out
}


plot_tank(tank_asv_models$model[[1]], tank_asv_models$data[[1]], tank_asv_models$all_health_fit[[1]])




plot_tank_fc <- function(model, data, field_fitted){
  
  ## FC preVpost (d & h)
  ## FC DvH post
  tank_contasts <- emmeans(model, ~ time_treat) %>%
    contrast(pathogen_test) %>%
    broom::tidy(conf.int = TRUE) %>%
    mutate(time = if_else(contrast == 'DvH', 'Post', 'Pre.v.Post'),
           health = if_else(contrast == 'DvH', NA_character_, str_extract(contrast, 'D$|H$'))) %>%
    select(health, time, estimate, starts_with('conf'), p.value)
  
  field_contrasts <- contrast(field_fitted, 'pairwise') %>%
    broom::tidy(conf.int = TRUE) %>%
    mutate(time = 'Field') %>%
    select(time, estimate, starts_with('conf'), p.value)
  
  
  out <- bind_rows(field_contrasts, 
            tank_contasts) %>%
    mutate(time = factor(time, levels = c('Field', 'Pre.v.Post', 'Post'))) %>%
    ggplot(aes(x = time, y = estimate, colour = health)) +
    geom_hline(yintercept = 0) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge2(0.5), 
                  linewidth = 1.5)

  out
}


plot_prePost <- function(model, data, field_fitted){
  initial_emmean <- emmeans(model, ~time_treat)
  treatment_levels <- levels(initial_emmean)$time_treat
  
  
  ref_grid(model) %>%
    add_grouping('treat_out', 'time_treat', 
                 str_remove(treatment_levels, '[0-9]+_')) %>%
    emmeans(~ treat_out) %>%
    broom::tidy(conf.int = TRUE) %>%
    mutate(prePost = if_else(str_detect(treat_out, 'pre'), 'Pre', 'Post'),
           prePost = factor(prePost, levels = c('Field', 'Pre', 'Post')),
           
           exposure = str_extract(treat_out, '^[DN]'),
           health = str_extract(treat_out, '[DH]$')) %>%
    bind_rows(broom::tidy(field_fitted, conf.int = TRUE) %>%
                mutate(prePost = factor('Field', levels = c('Field', 'Pre', 'Post')))) %>%
    ggplot(aes(x = prePost, y = estimate, colour = health, shape = exposure)) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    position = position_dodge2(0.5)) +
    scale_shape_manual(values = c('D' = 'circle', 'N' = 'triangle'), na.value = 'square')
    
}




plot_prePost(tank_asv_models$model[[1]], tank_asv_models$data[[1]], tank_asv_models$all_health_fit[[1]])


tank_asv_models <- tank_asv_data %>%
  nest_by(asv_id) %>%
  left_join(select(field_asv_models, asv_id, just_health_fit, all_health_fit),
            by = 'asv_id') %>%
  # filter(asv_id == 'ASV25') %>%
  
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  
  rowwise %>%
  mutate(data = list(mutate(data, log2_cpm_norm = log2_cpm_norm + rnorm(nrow(data), 0, 6.11/100)))) %>%
  mutate(model = list(lmer(log2_cpm_norm ~ time_treat +  
                             (1 | geno) + (1 | tank_id),
                           data = data))) %>%
  mutate(plot_prePost = list(plot_prePost(model, data, field_fitted = all_health_fit) +
                               labs(y = 'log2(CPM)', x = NULL,
                                    title = str_c(order, family, genus, asv_id, sep = '\n')) +
                               theme_classic() +
                               theme(legend.position = 'left',
                                     panel.background = element_rect(colour = 'black'),
                                     axis.text = element_text(colour = 'black'),
                                     strip.background = element_blank())))


wrap_plots(tank_asv_models$plot)
ggsave('../Results/testing_important_asvs_tank_1.png', height = 15, width = 15)


wrap_plots(tank_asv_models$plot_fc)
ggsave('../Results/testing_important_asvs_tank_2.png', height = 15, width = 15)

wrap_plots(tank_asv_models$plot_prePost, guides = 'collect')
ggsave('../Results/testing_important_asvs_tank_3.png', height = 15, width = 15)


select(field_asv_models, asv_id, starts_with('p.within')) %>%
  mutate(across(starts_with('p.within'), ~p.adjust(., method = 'fdr'))) %>%
  rowwise %>%
  filter(all(c_across(starts_with('p.within')) < 0.05)) %>%
  ungroup %>%
  select(-starts_with('p.within')) %>%
  inner_join(tank_asv_models,
             by = 'asv_id') %>%
  filter(asv_id == 'ASV51') %>%
  pull(plot_prePost) %>%
  wrap_plots(guides = 'collect')


extract_significance(model, anova_summary = 'anova')
