read_csv('~/../Downloads/normalized_tank_asv_counts.csv', 
         show_col_types = FALSE) %>%
  filter(asv_id == 'ASV25') %>%
  count(health)


tmp_data <- read_csv('~/../Downloads/normalized_tank_asv_counts.csv', 
                     show_col_types = FALSE) %>%
  filter(asv_id == 'ASV25') %>%
  filter(geno != 'GE') %>%
  select(-health) %>%
  # group_by(sample_id, time, tank, anti_health, geno, resist) %>%
  separate(anti_health, into = c('with_anti', 'disease_exposure')) %>%
  rename(health = resist) %>%
  filter(with_anti == 'N') %>%
  mutate(time = str_extract(time, '[0-9]') %>% as.integer(),
         fragment_id = str_c(geno, disease_exposure, tank, sep = '_')) %>%
  select(asv_id, log2_cpm_norm, sample_id, fragment_id, time, disease_exposure,
         geno, health)

tmp_data %>%
  count(fragment_id, time) %>%
  pivot_wider(names_from = time, values_from = n)


whacky <- c('Bin1_N_D_YE', 'Bin2_N_D_YE', 'Bin5_N_D_OR', 'Bin5_N_D_GR')
tank_asv_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                          show_col_types = FALSE) %>% 
  filter(asv_id == 'ASV25') %>%
  filter(anti == 'N') %>% 
  mutate(time = str_extract(time, '[0-9]') %>% as.integer(),
         fragment_id = str_c(geno, exposure, str_extract(tank, 'Bin[0-9]+'), sep = '_'),
         tank_id = str_c(exposure, str_extract(tank, 'Bin[0-9]+'), sep = '_'),
         treatment = if_else(exposure == 'pre', 'pre',
                             str_c(exposure, health, sep = '_'))) %>%
  select(-dataset, -anti, -plate, -year:-site, -cpm_norm:-n_reads, 
         -lib.size, -norm.factors, -domain:-species) %>%
  
  #Get rid of samples which shouldnt exist
  filter(!(str_detect(sample_id, str_c(whacky, collapse = '|')) & time == 8))


tank_asv_data %>%
  count(fragment_id, time) %>%
  arrange(fragment_id, time) %>%
  pivot_wider(names_from = time, values_from = n)

tank_asv_data %>%
  count(tank_id, time) %>%
  pivot_wider(names_from = time, values_from = n)

tank_asv_data %>%
  count(geno, time) %>%
  pivot_wider(names_from = time, values_from = n)

count(tank_asv_data, exposure)

lmer(log2_cpm_norm ~ (1 | geno) + (1 | tank_id),
     data = tank_asv_data)


tmp_data2 <- tank_asv_data %>% 
  rename(disease_exposure = exposure) %>%
  mutate(exposure = as.integer(time == 0),
         exp_dis = as.integer(time != 0) * as.integer(disease_exposure == 'D'),
         exp_hea = as.integer(time != 0) * as.integer(disease_exposure == 'N')) %>%
  select(asv_id, sample_id, log2_cpm_norm, exposure, exp_dis, exp_hea, health,
         time, geno, tank_id) %>%
  mutate(dose = case_when(exp_dis == 1 ~ 'D',
                          exp_hea == 1 ~ 'H',
                          TRUE ~ NA_character_),
         time_c = time + 1)

tmp_data2 %>%
  filter(time == 2,
         health == 'D')


tst_model <- lmer(log2_cpm_norm ~ health + (exp_dis + exp_hea) + 
                    (1 | geno) + (1 | tank_id),
                  data = tmp_data2)
summary(tst_model, ddf = "Kenward-Roger")
anova(tst_model, ddf = "Kenward-Roger")

emmeans(tst_model, ~ health + (exp_dis + exp_hea)) %>%
  as_tibble() %>%
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
  
  ggplot(aes(x = time_c - 1, y = emmean, colour = dose,
             group = grouping)) +
  geom_vline(xintercept = 1.75) +
  
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = dose),
              colour = NA, alpha = 0.5) +
  geom_line() +
  # geom_point(data = test_data, aes(y = log2_cpm_norm), size = 3) +
  geom_pointrange(data = tmp_data2 %>%
                    group_by(time_c, health, dose) %>%
                    summarise(mean_log2 = mean(log2_cpm_norm),
                              se_log2 = sd(log2_cpm_norm) / sqrt(n()),
                              n = n(),
                              .groups = 'drop'),
                  aes(y = mean_log2,
                      ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2,
                      group = 1),
                  size = 3, linewidth = 2,
                  position = position_dodge2(0.5)) +
  scale_x_continuous(breaks = c(0, 2, 8)) +
  facet_wrap(~health)




#### Try Again ####
tank_asv_data

tst_model <- lmer(log2_cpm_norm ~ treatment + 
                    (1 | geno) + (1 | tank_id),
                  data = tank_asv_data)
summary(tst_model, ddf = "Kenward-Roger")
anova(tst_model, ddf = "Kenward-Roger")

emmeans(tst_model, ~ treatment) %>%
  as_tibble() %>%
  mutate(time = if_else(treatment == 'pre', 'pre', 'post'),
         time = factor(time, levels = c('pre', 'post'))) %>%
  separate(treatment, into = c('exposure', 'health'), sep = '_', fill = 'right') %>%
  mutate(health = if_else(is.na(health), 'H', health)) %>%
  ggplot(aes(x = time, y = emmean, colour = health, shape = exposure)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, position = position_dodge(0.5)) +
  geom_point(size = 2.5, position = position_dodge(0.5))



#### One more Try ####
tst_data <- mutate(tank_asv_data, time_treat = str_c(time, exposure, health, sep = '_'))
count(tst_data, time_treat)

tst_model <- lmer(log2_cpm_norm ~ time_treat +  
                    (1 | geno) + (1 | tank_id),
                  data = tst_data)
summary(tank_asv_models$model[[2]], ddf = "Kenward-Roger")
  

filter(tank_asv_models, asv_id == 'ASV25') %>%
  select(model) %>%
  pull(model) %>%
  pluck(1) %>%
  emmeans(~ time_treat) %>%
  as_tibble() %>%
  separate(time_treat, into = c('time', 'exposure', 'health'), 
           sep = '_', fill = 'right') %>%
  mutate(time = as.integer(time)) %>%
  ggplot(aes(x = time, y = emmean, colour = health, shape = exposure)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, position = position_dodge(0.5)) +
  geom_point(size = 2.5, position = position_dodge(0.5))



filter(tank_asv_models, asv_id == 'ASV10') %>%
  select(model) %>%
  pull(model) %>%
  pluck(1) %>%
  ref_grid(~ time_treat) %>%
  add_grouping('time_disease', 'time_treat', 
               c('0_H', '2_D', '2_H', '2_H', 
                 '8_D', '8_H', '8_D', '8_H')) %>%
  emmeans(~ time_disease) %>%
  as_tibble() %>%
  separate(time_disease, into = c('time', 'health'), 
           sep = '_', fill = 'right') %>%
  mutate(time = as.integer(time)) %>%
  ggplot(aes(x = time, y = emmean, colour = health)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, position = position_dodge(0.5)) +
  geom_point(size = 2.5, position = position_dodge(0.5))


filter(tank_asv_models, asv_id == 'ASV10') %>%
  select(model) %>%
  pull(model) %>%
  pluck(1) %>%
  ref_grid(~ time_treat) %>%
  add_grouping('time_disease', 'time_treat', 
               c('0_H', '2_D', '2_H', '2_H', 
                 '8_D', '8_H', '8_D', '8_H')) %>%
  add_grouping('healthy_disease', 'time_disease', 
               c('H', 'D', 'H', 'D', 'H')) %>%
  emmeans(~ healthy_disease) %>%
  contrast('pairwise')
