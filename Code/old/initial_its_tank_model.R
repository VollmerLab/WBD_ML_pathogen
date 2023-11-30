#### Model ASV change through time ####
#Just use data with no antibiotic exposure and in samples exposed to the disease
test_data <- tank_raw %>%
  # filter(n_reads > 0) %>%
  filter(with_anti == 'N', (disease_exposure == 'D' | time == 'T0')) %>%
  select(asv_id, geno, time, health, tank, fragment_id, where(is.numeric)) %>%
  
  filter(asv_id == 'ASV25') %>%
  mutate(time_c = str_extract(time, '[0-9]+') %>% as.integer,
         time_c = time_c + 1,
         disease_exposure = as.integer(time != 'T0'),
         
         post_exposure_time = time_c * disease_exposure,
         post_exposure_time = if_else(post_exposure_time == 0, 0, 
                                      post_exposure_time - min(post_exposure_time[post_exposure_time > 0]) + 1))

count(test_data, geno, tank, time) %>%
  pivot_wider(names_from = time,
              values_from = n)

test_data %>%
  group_by(geno, time) %>%
  summarise(the_sd = sd(log2_cpm_norm))


its_model <- lm(log2_cpm_norm ~ resist * (time_c + disease_exposure),
                data = test_data)
summary(its_model)

its_model <- lmer(log2_cpm_norm ~ resist * (time_c + disease_exposure) + (1 | tank),
                  data = test_data)
summary(its_model)
anova(its_model)

emmeans(its_model, ~resist * (disease_exposure)) %>%
  as_tibble() %>%
  rowwise %>%
  mutate(time_c = list(time_c = seq(1, 9, length.out = 100))) %>%
  unnest(time_c) %>%
  filter(!(disease_exposure == 1 & time_c < 3 + 0.25),
         !(disease_exposure == 0 & time_c > 3 - 0.25)) %>%
  ggplot(aes(x = time_c - 1, y = emmean, colour = resist, group = interaction(resist, disease_exposure))) +
  
  geom_vline(xintercept = 1.75) +
  
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = resist),
              colour = NA, alpha = 0.5) +
  geom_line() +
  # geom_point(data = test_data, aes(y = log2_cpm_norm), size = 3) +
  geom_pointrange(data = test_data %>%
                    group_by(time_c, resist) %>%
                    summarise(mean_log2 = mean(log2_cpm_norm),
                              se_log2 = sd(log2_cpm_norm) / sqrt(n()),
                              .groups = 'drop'),
                  aes(y = mean_log2,
                      ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2,
                      group = 1),
                  size = 3, linewidth = 2,
                  position = position_dodge2(0.5)) +
  scale_x_continuous(breaks = c(0, 2, 8))



emmeans(its_model, ~resist * (time_c + disease_exposure), at = list(time_c = seq(1, 9, length.out = 100))) %>%
  as_tibble() %>%
  # rowwise %>%
  # mutate(time_c = if_else(disease_exposure == 0, list(1:3), list(3:9))) %>%
  # unnest(time_c) %>%
  filter(!(disease_exposure == 1 & time_c < 3 + 0.25),
         !(disease_exposure == 0 & time_c > 3 - 0.25)) %>%
  
  ggplot(aes(x = time_c - 1, y = emmean, colour = resist, group = interaction(resist, disease_exposure))) +
  
  geom_vline(xintercept = 1.75) +
  
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = resist),
              colour = NA, alpha = 0.5) +
  geom_line() +
  # geom_point(data = test_data, aes(y = log2_cpm_norm), size = 3) +
  geom_pointrange(data = test_data %>%
                    group_by(time_c, resist) %>%
                    summarise(mean_log2 = mean(log2_cpm_norm),
                              se_log2 = sd(log2_cpm_norm) / sqrt(n()),
                              .groups = 'drop'),
                  aes(y = mean_log2,
                      ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2,
                      group = 1),
                  size = 3, linewidth = 2,
                  position = position_dodge(0.5)) +
  scale_x_continuous(breaks = c(0, 2, 8))



#### Full Models ####
plot_results <- function(model, data, variable_time, buffer = 0.25){
  if(variable_time){
    plot_data <- emmeans(model, ~health * (time_c + disease_exposure), 
                         at = list(time_c = seq(1, 9, length.out = 100))) %>%
      as_tibble() %>%
      filter(!(disease_exposure == 1 & time_c < 3 + buffer),
             !(disease_exposure == 0 & time_c > 3 - buffer))
  } else {
    plot_data <- emmeans(model, ~health * (disease_exposure)) %>%
      as_tibble() %>%
      rowwise %>%
      mutate(time_c = list(time_c = seq(1, 9, length.out = 100))) %>%
      unnest(time_c) %>%
      filter(!(disease_exposure == 1 & time_c < 3 + buffer),
             !(disease_exposure == 0 & time_c > 3 - buffer))
  }
  
  output_stats <- anova(model, ddf = "Kenward-Roger") %>%
    as_tibble(rownames = 'effect') %>%
    mutate(effect = str_replace_all(effect, c('_c' = '', '_' = ' ', ':' = ' x ')),
           effect = str_to_title(effect) %>% str_replace('X', 'x'),
           across(c(where(is.numeric), -`F value`, -`Pr(>F)`), ~round(., 1)),
           `F value` = round(`F value`, 3),
           output = str_c(effect, ': F(', NumDF, ', ', DenDF, ') = ', 
                          `F value`, '; p = ', scales::pvalue(`Pr(>F)`))) %>%
    pull(output) %>%
    str_c(collapse = '\n')
  
  plot_data %>%
    mutate(health = if_else(health == 'H', 'Healthy', 'Diseased')) %>%
    ggplot(aes(x = time_c - 1, y = emmean, colour = health, group = interaction(health, disease_exposure))) +
    
    geom_vline(xintercept = 1.75) +
    
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = health),
                colour = NA, alpha = 0.5, show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    # geom_point(data = test_data, aes(y = log2_cpm_norm), size = 3) +
    geom_pointrange(data = data %>%
                      mutate(health = if_else(health == 'H', 'Healthy', 'Diseased')) %>%
                      group_by(time_c, health) %>%
                      summarise(mean_log2 = mean(log2_cpm_norm),
                                se_log2 = sd(log2_cpm_norm) / sqrt(n()),
                                .groups = 'drop'),
                    aes(y = mean_log2,
                        ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2,
                        group = 1),
                    size = 3, linewidth = 2,
                    position = position_dodge(0.5)) +
    scale_x_continuous(breaks = c(0, 2, 8)) +
    guides(colour = guide_legend(override.aes = list(size = 1, linewidth = 0))) +
    labs(colour = NULL,
         fill = NULL,
         subtitle = output_stats)
}

select(tank_raw, geno, tank, time, fragment_id) %>%
  distinct %>%
  count(fragment_id)

select(tank_raw, geno, tank, time) %>%
  distinct %>%
  count(geno, tank, time) %>%
  pivot_wider(names_from = time,
              values_from = n) %T>%
  print() %>%
  group_by(geno) %>%
  summarise(across(where(is.numeric), ~sum(., na.rm = TRUE)))

its_model_analysis <- tank_raw %>%
  # filter(n_reads > 0) %>%
  # filter(with_anti == 'N', (disease_exposure == 'D' | time == 'T0')) %>%
  select(asv_id, geno, time, health, tank, fragment_id, where(is.numeric)) %>%
  mutate(time_c = str_extract(time, '[0-9]+') %>% as.integer,
         time_c = time_c + 1,
         disease_exposure = as.integer(time != 'T0'),
         
         post_exposure_time = time_c * disease_exposure,
         post_exposure_time = if_else(post_exposure_time == 0, 0, 
                                      post_exposure_time - min(post_exposure_time[post_exposure_time > 0]) + 1)) %>%
  nest_by(asv_id) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  filter(family %in% shap_importance$family) %>%
  mutate(important_asv = asv_id %in% shap_importance$asv_id) %>%
  filter(!(family == 'Endozoicomonadaceae' & !important_asv),
         !(family == 'Rhodobacteraceae' & !important_asv)) %>% #ungroup %>% count(family)
  filter(important_asv) %>%
  # inner_join(shap_importance, by = 'asv_id') %>%
  mutate(useful = str_c(family, genus, asv_id, sep = '_') %>%
           str_replace_all('__', '_')) %>% 
  # mutate(model = list(lmer(log2_cpm_norm ~ health * (time_c + disease_exposure) + (1 | tank) + (1 | geno),
  #                          data = data))) %>%
  # mutate(model = list(lmer(log2_cpm_norm ~ health * (time_c + disease_exposure) + (1 | fragment_id),
  #                          data = data))) %>%
  
  mutate(model = list(lmer(log2_cpm_norm ~ health * (disease_exposure) + (1 | fragment_id),
                           data = data)),
         pvalue_time_c = 1,
         pvalue_healthHXtime_c = 1) %>%
  
  mutate(summary(model, ddf = "Kenward-Roger") %>% 
           coef %>% 
           as_tibble(rownames = 'param') %>% 
           rename(beta = Estimate,
                  pvalue = `Pr(>|t|)`) %>%
           select(param, beta, pvalue) %>%
           mutate(param = str_remove_all(param, '\\(|\\)') %>% str_replace_all(':', 'X')) %>%
           pivot_wider(names_from = param, values_from = c(beta, pvalue)),
         
         plot = list(plot_results(model, data, 
                                  variable_time = (pvalue_time_c < alpha |
                                                     pvalue_healthHXtime_c < alpha)) +
                       labs(x = 'Time', 
                            y = 'log2(CPM)',
                            title = useful) +
                       theme_classic())) %>%
  ungroup

its_model_analysis %>%
  # mutate(across(starts_with('pvalue'), ~p.adjust(., method = 'fdr'))) %>%
  # select(starts_with('pvalue'))
  filter(pvalue_healthHXdisease_exposure < alpha | 
           pvalue_healthH < alpha | 
           pvalue_healthHXtime_c < alpha) %>%
  # select(starts_with('pvalue'))
  pull(plot) %>%
  wrap_plots() + 
  plot_layout(guides = 'collect')
ggsave('../Results/significant_associated_asvs.png', height = 15, width = 15)


its_model_analysis %>%
  filter(!(pvalue_healthHXdisease_exposure < alpha | 
             pvalue_healthH< alpha | 
             pvalue_healthHXtime_c < alpha)) %>%
  pull(plot) %>%
  wrap_plots() + 
  plot_layout(guides = 'collect')
ggsave('../Results/nonsignificant_associated_asvs.png', height = 21, width = 21)


#### Test Clustering of ASVs ####
library(tidymodels)

asv_pca <- select(its_model_analysis, asv_id, starts_with('beta'), starts_with('pvalue')) %>%
  ungroup %>%
  select(-contains('Intercept')) %>%
  pivot_longer(cols = -asv_id,
               names_to = c('.value', 'coef'),
               names_pattern = c('(beta|pvalue)_(.*)')) %>%
  mutate(beta = if_else(pvalue > alpha, 0, beta)) %>%
  select(-pvalue) %>%
  pivot_wider(names_from = coef, 
              values_from = beta) %>%
  recipe(~.) %>%
  update_role(asv_id, new_role = "ID") %>%
  step_pca(all_predictors()) %>%
  prep


tidy(asv_pca, 1) %>%
  select(-id) %>%
  pivot_wider(names_from = component, 
              values_from = value) %>%
  mutate(terms = str_remove(terms, '^beta_')) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_segment(xend = 0, yend = 0) +
  geom_text(aes(label = terms)) +
  geom_point(data = left_join(juice(asv_pca),
                              select(its_model_analysis, asv_id, starts_with('pvalue'))),
             aes(colour = pvalue_resistSXdisease_exposure < alpha))



tank_data