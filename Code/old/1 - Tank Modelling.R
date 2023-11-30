#### Libraries ####
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(patchwork)

#### Data ####
tank_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                      show_col_types = FALSE) %>%
  mutate(across(domain:species, str_replace_na)) %>%
  select(-site, -dataset, -year, -season, -anti, -health) %>%
  rename(time_treat = time) %>%
  mutate(time = str_extract(time_treat, '[0-9]') %>% str_c('T', .)) %>%
  filter(geno != 'GE') %>%
  separate(anti_health, into = c('with_anti', 'disease_exposure')) %>%
  mutate(fragment_id = str_c(geno, tank, with_anti, disease_exposure, sep = '_'),
         tank_id = str_c(tank, with_anti, disease_exposure, sep = '_'))
  # mutate(across(c(with_anti, disease_exposure), ~. != 'N'))

### Investigate metadata ####
tank_data %>%
  select(disease_exposure, time, tank, with_anti, geno, resist, plate) %>%
  distinct %>%
  mutate(across(where(is.logical), as.character)) %>%
  pivot_longer(cols = everything(),
               names_to = 'variable',
               values_to = 'value') %>%
  nest_by(variable) %>%
  reframe(count(data, value)) %>%
  group_split(variable)

#### Sort out experimental design ####
tank_data %>%
  group_by(sample_id, across(time_treat:plate), time, fragment_id, tank_id) %>%
  summarise(mean_l2cpm = mean(log2_cpm_norm),
            .groups = 'drop') %>%
  select(sample_id, tank_id, geno, fragment_id, time, with_anti, disease_exposure) %>%
  # muta
  pivot_wider(names_from = 'time',
              values_from = 'sample_id') %>%
  arrange(geno) %>%
  select(geno, with_anti, disease_exposure, tank_id, fragment_id, starts_with('T')) %>%
  View

#### Simple Plot ####
tank_data %>%
  group_by(sample_id, across(time_treat:plate), time) %>%
  summarise(mean_l2cpm = mean(log2_cpm_norm),
            .groups = 'drop') %>%
  ggplot(aes(x = time, y = mean_l2cpm, colour = resist, shape = geno)) +
  geom_point() +
  facet_grid(with_anti + disease_exposure ~ tank)

#### Compare Anti-biotic effect (T0) - no disease exposure has happened ####
antibiotic_effects <- tank_data %>%
  filter(time == 'T0') %>%
  nest_by(asv_id, across(domain:species)) %>%
  mutate(data2 = list(select(data, tank, geno, log2_cpm_norm, with_anti) %>%
                       pivot_wider(names_from = with_anti,
                                   values_from = 'log2_cpm_norm') %>%
                       filter(!is.na(A), !is.na(N)) %>%
                       mutate(log_fc = A - N)),
         
         full_model = list(lmer(log_fc ~ (1 | geno), data = data2)),
         
         summary(full_model)$coef %>% 
           as_tibble() %>% 
           rename(mean_fc = Estimate,
                  se_fc = `Std. Error`,
                  t_value = `t value`,
                  p_value = `Pr(>|t|)`)) %>%
  ungroup %>%
  mutate(fdr = p.adjust(p_value, method = 'fdr'))

filter(antibiotic_effects, fdr < 0.1)

#### Look at antibiotic effect over time only in disease exposed samples ####
tmp$data[[1]]

tmp <- tank_data %>%
  filter(disease_exposure == 'D') %>%
  nest_by(asv_id, across(domain:species)) %>%
  mutate(full_model = list(lmer(log2_cpm_norm ~ with_anti * time + (1 | geno / fragment_id), data = data)),
         
         full_model %>%
           anova(ddf = 'Kenward-Roger') %>%
           as_tibble(rownames = 'term') %>%
           select(term, `Pr(>F)`) %>%
           pivot_wider(names_from = 'term',
                       values_from = `Pr(>F)`) %>%
           rename_with(~str_replace_all(., ':', 'X'))) %>%
  ungroup %>%
  mutate(across(c(with_anti, time, with_antiXtime), ~p.adjust(., method = 'fdr')))

anti_time_interaction <- filter(tmp, with_antiXtime < 0.05 & with_anti < 0.05) %>%
  mutate(across(domain:species, str_replace_na)) %>%
  rowwise %>%
  mutate(plot = list(emmeans(full_model, ~time * with_anti) %>%
                       tidy(conf.int = TRUE) %>%
                       ggplot(aes(x = str_extract(time, '[0-9]') %>% as.integer(), 
                                  y = estimate, 
                                  colour = with_anti)) +
                       geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                                       position = position_dodge(0.5)) +
                       labs(x = time, y = 'log2CPM',
                            title = str_c(family, genus, asv_id, sep = '_')))) %>%
  pull(plot) %>%
  wrap_plots()

filter(tmp, with_anti < 0.05) %>%
  sample_n(12) %>%
  mutate(across(domain:species, str_replace_na)) %>%
  rowwise %>%
  mutate(plot = list(emmeans(full_model, ~time * with_anti) %>%
                       tidy(conf.int = TRUE) %>%
                       ggplot(aes(x = str_extract(time, '[0-9]') %>% as.integer(), 
                                  y = estimate, 
                                  colour = with_anti)) +
                       geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                                       position = position_dodge(0.5)) +
                       labs(x = time, y = 'log2CPM',
                            title = str_c(family, genus, asv_id, sep = '_')))) %>%
  ungroup %>%
  pull(plot) %>%
  wrap_plots()

tmp$full_model[[1]] %>%
  emmeans(~time * with_anti) %>%
  tidy(conf.int = TRUE) %>%
  ggplot(aes(x = str_extract(time, '[0-9]') %>% as.integer(), y = estimate, 
             colour = with_anti)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(0.5)) 

#### Compare fully crossed antibiotic & exposure effects @ T2 & T8 ####
anti_disease <- tank_data %>%
  filter(time != 'T0') %>%
  nest_by(asv_id, across(domain:species)) %>%
  mutate(full_model = list(lmer(log2_cpm_norm ~ disease_exposure * with_anti * time + (1 | geno / fragment_id), data = data)),
         
         full_model %>%
           anova(ddf = 'Kenward-Roger') %>%
           as_tibble(rownames = 'term') %>%
           select(term, `Pr(>F)`) %>%
           pivot_wider(names_from = 'term',
                       values_from = `Pr(>F)`) %>%
           rename_with(~str_replace_all(., ':', 'X'))) %>%
  ungroup 

anti_disease %>%
  mutate(across(disease_exposure:disease_exposureXwith_antiXtime, ~p.adjust(., method = 'fdr'))) %>%
  filter(if_any(c(contains('disease_exposure'), contains('with_anti')), ~. < 0.05)) %>%
  
  filter(disease_exposureXwith_antiXtime < 0.05) %>%
  
  rowwise %>%
  mutate(plot = list(emmeans(full_model, ~disease_exposure * with_anti * time) %>%
                       tidy(conf.int = TRUE) %>%
                       ggplot(aes(x = time, y = estimate, colour = disease_exposure, shape = with_anti)) +
                       geom_pointrange(aes(ymin = conf.low, ymax = conf.high), 
                                       position = position_dodge(0.5)) +
                       labs(x = 'time', y = 'log2CPM',
                            title = str_c(family, genus, asv_id, sep = '_')))) %>%
  pull(plot) %>%
  wrap_plots()



#### Model ASVs ####
tmp <- tank_data %>%
  mutate(time = str_extract(time, '[0-9]+') %>% as.integer) %>%
  # mutate(time = factor(time, ordered = TRUE)) %>%
  nest_by(asv_id, across(domain:species))
tmp$data[[1]]$time

tst <- lmer(log2_cpm_norm ~ time * with_anti * disease_exposure + 
       (1 | geno / fragment_id),
     data = tmp$data[[1]]) 

anova(tst, ddf = "Kenward-Roger", type = '3')
emmeans(tst, ~time * with_anti * disease_exposure, 
        at = list(time = c(0, 2, 8))) %>%
  tidy(conf.int = TRUE) %>%
  filter(!(time == 0 & disease_exposure == 'N')) %>%
   
  ggplot(aes(x = time, y = estimate, colour = disease_exposure, shape = with_anti)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(1), size = 1) +
  geom_point(data = tmp$data[[1]], aes(y = log2_cpm_norm),
             position = position_dodge(1), size = 1)


