tmp_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
         show_col_types = FALSE) %>%
  filter(geno != 'GE') %>%
  select(-health) %>%
  # group_by(sample_id, time, tank, anti_health, geno, resist) %>%
  separate(anti_health, into = c('with_anti', 'disease_exposure')) %>%
  rename(health = resist) %>%
  filter(with_anti == 'N') %>%
  mutate(time = str_extract(time, '[0-9]') %>% as.integer(),
         fragment_id = str_c(geno, disease_exposure, tank, sep = '_')) %>%
  filter(asv_id == 'ASV30')


tmp_data %>%
  count(fragment_id, time) %>%
  pivot_wider(names_from = time, values_from = n)
  


tmp_data %>%
  ggplot(aes(x = time, y = log2_cpm_norm, colour = health)) +
  geom_point() +
  facet_grid(~ disease_exposure)


tmp_data2 <- tmp_data %>% 
  mutate(exposure = as.integer(time == 0),
         exp_dis = as.integer(time != 0) * as.integer(disease_exposure == 'D'),
         exp_hea = as.integer(time != 0) * as.integer(disease_exposure == 'N')) %>%
  select(asv_id, sample_id, log2_cpm_norm, exposure, exp_dis, exp_hea, health,
         time, fragment_id) %>%
  mutate(dose = case_when(exp_dis == 1 ~ 'D',
                          exp_hea == 1 ~ 'H',
                          TRUE ~ NA_character_),
         time_c = time + 1)

count(tmp_data2, health, exp_dis, exp_hea, dose)

tst_model <- lmer(log2_cpm_norm ~ health * (exp_dis + exp_hea) + (1 | fragment_id),
     data = tmp_data2)
summary(tst_model, ddf = "Kenward-Roger")
anova(tst_model, ddf = "Kenward-Roger")


emmeans(tst_model, ~health * (exp_dis + exp_hea)) %>%
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



#### All asvs ####
asv_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
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
         exp_hea = as.integer(time != 0) * as.integer(disease_exposure == 'N')) %>%
  select(asv_id, sample_id, log2_cpm_norm, exposure, exp_dis, exp_hea, health,
         time, fragment_id) %>%
  mutate(dose = case_when(exp_dis == 1 ~ 'D',
                          exp_hea == 1 ~ 'H',
                          TRUE ~ NA_character_),
         time_c = time + 1) %>%
  nest_by(asv_id) 


asv_data <- ungroup(asv_data) %>%
  sample_n(10) %>%
  bind_rows(filter(asv_data, asv_id == 'ASV25')) %>%
  distinct %>%
  rowwise(asv_id)

asv_data <- filter(asv_data, asv_id %in% shap_importance$asv_id)

library(rstanarm)
test <- asv_data %>%
  mutate(model = list(stan_lmer(log2_cpm_norm ~ health * (exp_dis + exp_hea) + (1 | fragment_id),
                                data = data,
                                chains = 4, cores = 4)),
         plot = list(emmeans(model, ~health * (exp_dis + exp_hea),
                             data = data) %>%
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
                       
                       geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD, fill = dose),
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
                                       size = 3, linewidth = 2,
                                       position = position_dodge2(0.5)) +
                       scale_x_continuous(breaks = c(0, 2, 8)) +
                       facet_wrap(~health)))


blat <- summary(test$model[[1]])
coef(test$model[[1]])
fixef(test$model[[1]])
unname(diff(summary(test$model[[1]])['healthS',c('10%', '90%')] > 0) == 0)

test %>%
  mutate(health_alone = unname(diff(summary(model)['healthS',c('10%', '90%')] > 0) == 0),
         dis_alone = unname(diff(summary(model)['exp_dis',c('10%', '90%')] > 0) == 0),
         dis_interaction = unname(diff(summary(model)['healthS:exp_dis',c('10%', '90%')] > 0) == 0)) %>%
  filter(health_alone | dis_alone | dis_interaction) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(plot = list(plot + labs(title = str_c(order, family, genus, asv_id, sep = '\n')))) %>%
  pull(plot) %>%
  wrap_plots()
ggsave('../Results/test_its_with_control.png', height = 7, width = 10)
