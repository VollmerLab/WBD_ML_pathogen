#### Libraries ####
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(patchwork)

#### N ASVs ####
N_ASV <- 10
alpha <- 0.05

#### Data ####
model_list <- read_csv('../Results/equivilant_top_models.csv.gz', show_col_types = FALSE) %>%
  filter(pract_equiv >= 0.95) %>%
  rename(wflow_id = model) %>%
  pull(wflow_id)

taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz', show_col_types = FALSE)
shap_importance <- read_csv('../intermediate_files/model_shaps_importance.csv.gz', 
                            show_col_types = FALSE) %>%
  filter(wflow_id %in% model_list,
         asv_rank <= N_ASV) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(asv_id = fct_reorder(asv_id, asv_rank)) %>%
  select(asv_id, domain:genus) %>%
  distinct %>%
  arrange(asv_id)

tank_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                      show_col_types = FALSE) %>%
  mutate(across(domain:species, str_replace_na)) %>%
  select(-site, -dataset, -year, -season, -anti, -health) %>%
  rename(time_treat = time) %>%
  mutate(time = str_extract(time_treat, '[0-9]') %>% str_c('T', .)) %>%
  filter(geno != 'GE') %>%
  separate(anti_health, into = c('with_anti', 'disease_exposure')) %>%
  mutate(disease_exposure = if_else(time == 'T0', 'N', disease_exposure)) %>%
  mutate(fragment_id = str_c(geno, tank, with_anti, disease_exposure, sep = '_'),
         tank_id = str_c(tank, with_anti, disease_exposure, sep = '_')) %>%
  
  # filter(n_reads > 0) %>%
  
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(asv_id = factor(asv_id, levels = levels(shap_importance$asv_id)))


#### Plots ####
tank_data %>%
  group_by(asv_id, time, with_anti, disease_exposure, geno, resist) %>%
  # filter(asv_id == 'ASV25') %>%
  summarise(mean_count = mean(log2_cpm_norm),
            se_count = sd(log2_cpm_norm) / sqrt(n()),
            .groups = 'drop') %>%
  filter(with_anti == 'N', (disease_exposure == 'D' | time == 'T0')) %>%
  
  ggplot(aes(x = time, y = mean_count, ymin = mean_count - se_count, ymax = mean_count + se_count, 
             colour = resist, 
             group = interaction(geno, resist))) +
  geom_path() +
  geom_pointrange() +
  facet_wrap(~asv_id, scales = 'free_y')


#### Model ASV change through time ####
#Just use data with no antibiotic exposure and in samples exposed to the disease
test_data <- tank_data %>%
  # filter(n_reads > 0) %>%
  filter(with_anti == 'N', (disease_exposure == 'D' | time == 'T0')) %>%
  select(asv_id, geno, time, resist, tank, fragment_id, where(is.numeric)) %>%
  
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

library(lmerTest)
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
    plot_data <- emmeans(model, ~resist * (time_c + disease_exposure), 
                         at = list(time_c = seq(1, 9, length.out = 100))) %>%
      as_tibble() %>%
      filter(!(disease_exposure == 1 & time_c < 3 + buffer),
             !(disease_exposure == 0 & time_c > 3 - buffer))
  } else {
    plot_data <- emmeans(model, ~resist * (disease_exposure)) %>%
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
    mutate(health = if_else(resist == 'R', 'Healthy', 'Diseased')) %>%
    ggplot(aes(x = time_c - 1, y = emmean, colour = health, group = interaction(health, disease_exposure))) +
    
    geom_vline(xintercept = 1.75) +
    
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = health),
                colour = NA, alpha = 0.5, show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    # geom_point(data = test_data, aes(y = log2_cpm_norm), size = 3) +
    geom_pointrange(data = data %>%
                      mutate(health = if_else(resist == 'R', 'Healthy', 'Diseased')) %>%
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

its_model_analysis <- tank_data %>%
  # filter(n_reads > 0) %>%
  filter(with_anti == 'N', (disease_exposure == 'D' | time == 'T0')) %>%
  select(asv_id, geno, time, resist, tank, fragment_id, where(is.numeric)) %>%
  mutate(time_c = str_extract(time, '[0-9]+') %>% as.integer,
         time_c = time_c + 1,
         disease_exposure = as.integer(time != 'T0'),
         
         post_exposure_time = time_c * disease_exposure,
         post_exposure_time = if_else(post_exposure_time == 0, 0, 
                                      post_exposure_time - min(post_exposure_time[post_exposure_time > 0]) + 1)) %>%
  nest_by(asv_id) %>%
  left_join(shap_importance, by = 'asv_id') %>%
  mutate(taxa = if_else(is.na(genus), str_c(family, asv_id, sep = '; '), 
                        str_c(family, genus, asv_id, sep = '; ')),
         taxa = fct_reorder(taxa, as.integer(asv_id))) %>% 
  mutate(model = list(lmer(log2_cpm_norm ~ resist * (time_c + disease_exposure) + (1 | tank) + (1 | geno),
                           data = data)),
         
         summary(model, ddf = "Kenward-Roger") %>% 
           coef %>% 
           as_tibble(rownames = 'param') %>% 
           rename(beta = Estimate,
                  pvalue = `Pr(>|t|)`) %>%
           select(param, beta, pvalue) %>%
           mutate(param = str_remove_all(param, '\\(|\\)') %>% str_replace_all(':', 'X')) %>%
           pivot_wider(names_from = param, values_from = c(beta, pvalue)),
         
         plot = list(plot_results(model, data, 
                                  variable_time = (pvalue_time_c < alpha |
                                                     pvalue_resistSXtime_c < alpha)) +
                       labs(x = 'Time', 
                            y = 'log2(CPM)',
                            title = taxa) +
                       theme_classic()))


its_model_analysis %>%
  filter(pvalue_resistSXdisease_exposure < alpha | 
           pvalue_resistS < alpha | 
           pvalue_resistSXtime_c < alpha) %>%
  pull(plot) %>%
  wrap_plots() + 
  plot_layout(guides = 'collect')
ggsave('../Results/significant_associated_asvs.png', height = 15, width = 15)


its_model_analysis %>%
  filter(!(pvalue_resistSXdisease_exposure < alpha | 
           pvalue_resistS < alpha | 
           pvalue_resistSXtime_c < alpha)) %>%
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

