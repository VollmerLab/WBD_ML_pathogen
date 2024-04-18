library(tidyverse)
library(chisq.posthoc.test)
library(ggtext)

pathogens <- str_c('ASV', c(8, 25, 38))


#### Data ####
field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  filter(asv_id %in% pathogens) %>%
  select(asv_id, sample_id, n_reads, health, site, season, year)

#### Binomial Model ####
library(blme)
library(emmeans)

binomial_data <- field_data %>%
  group_by(asv_id, site, season, year, health) %>%
  summarise(n = sum(n_reads > 0),
            total = n(),
            .groups = 'drop') %>%
  mutate(time = str_c(season, year, sep = '_'))

pathogen_binomial_model <- bglmer(cbind(n, total - n) ~ asv_id * health * time + (1 + asv_id | site),
       data = binomial_data,
       family = 'binomial',
       fixef.prior = normal(cov = diag(9, 24)), 
       control = glmerControl(optimizer = 'bobyqa'))

car::Anova(pathogen_binomial_model)

emmeans(pathogen_binomial_model, ~time | asv_id,
        at = list(health = 'D')) %>%
  contrast('poly')

emmeans(pathogen_binomial_model, ~time | asv_id,
        at = list(health = 'H')) %>%
  contrast('pairwise')

emmeans(pathogen_binomial_model, ~asv_id,
        at = list(health = 'D')) %>%
  contrast('pairwise')

emmeans(pathogen_binomial_model, ~asv_id,
        at = list(health = 'H')) %>%
  contrast('pairwise')

emmeans(pathogen_binomial_model, ~asv_id,
        at = list(health = 'D')) %>%
  contrast('pairwise')

emmeans(pathogen_binomial_model, ~health * time * asv_id, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  mutate(health = if_else(health == 'D', 'Diseased', 'Healthy')) %>%
  ggplot(aes(x = asv_id, y = prob, ymin = conf.low, ymax = conf.high, 
             fill = health, group = time, shape = time)) +
  geom_errorbar(width = 0.05, position = position_dodge(0.5), 
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), size = 3) +
  
  guides(fill = guide_legend(override.aes = list(shape = 'circle filled', size = 4)),
         shape = guide_legend(override.aes = list(size = 4, fill = 'black'))) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous")),
                                       c('Healthy', 'Diseased')),
                    breaks = c('Diseased', 'Healthy'), drop = FALSE) +
  scale_shape_manual(values = c('S_2016' = 'square filled',
                                'S_2017' = 'triangle filled',
                                'W_2016' = 'diamond filled',
                                'W_2017' = 'triangle down filled'),
                     breaks = c('W_2016', 'S_2016', 'W_2017', 'S_2017'),
                     labels = c('S_2016' = '16Jul',
                                'S_2017' = '17Jul',
                                'W_2016' = '16Jan',
                                'W_2017' = '17Jan'),
                     drop = FALSE) +
  
  labs(x = NULL,
       y = 'Detectable Quantity',
       shape = 'Time',
       fill = 'Disease\nState') +
  theme_classic() +
  theme(legend.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        legend.key = element_blank(),
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 14))


emmeans(pathogen_binomial_model, ~health * time * asv_id, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  mutate(asv_id = factor(asv_id, levels = str_c('ASV', c(25, 8, 38)))) %>%
  arrange(asv_id) %>%
  mutate(health = if_else(health == 'D', 'Diseased', 'Healthy'),
         time = factor(time, levels = c('W_2016', 'S_2016', 'W_2017', 'S_2017')),
         name = case_when(asv_id == 'ASV25' ~ '<i>Cysteiniphilum litorale</i>',
                          asv_id == 'ASV8' ~ '<i>Vibrio sp.</i>',
                          asv_id == 'ASV26' ~ 'Oceanospirillaceae',
                          asv_id == 'ASV30' ~ '<i>Thalassotalea sp.</i>',
                          asv_id == 'ASV361' ~ '<i>Endozoicomonas atrinae</i>',
                          asv_id == 'ASV38' ~ '<i>Neptuniibacter sp.</i>',
                          asv_id == 'ASV51' ~ '<i>Neptuniibacter sp.</i>',
                          asv_id == 'ASV40' ~ '<i>Qipengyuania sp.</i>',
                          TRUE ~ 'Unknown'),
         name = str_c(asv_id, name, sep = ' - '),
         name = fct_reorder(name, as.integer(asv_id))) %>%
  ggplot(aes(x = time, y = prob, ymin = conf.low, ymax = conf.high, 
             colour = health)) +
  geom_errorbar(width = 0.05, position = position_dodge(0.5), 
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), size = 3, shape = 'circle') +
  
  guides(colour = guide_legend(override.aes = list(shape = 'circle', size = 4))) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(labels = c('S_2016' = '16Jul', 'S_2017' = '17Jul', 'W_2016' = '16Jan', 'W_2017' = '17Jan')) +
  scale_colour_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous")),
                                       c('Healthy', 'Diseased')),
                    breaks = c('Diseased', 'Healthy'), drop = FALSE) +
  facet_wrap(~name) +
  labs(x = NULL,
       y = 'ASV Prevalence',
       colour = 'Disease\nState') +
  theme_classic() +
  theme(legend.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        legend.key = element_blank(),
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 14),
        strip.background = element_blank(),
        strip.text = element_markdown(hjust = 0, size = 16, colour = 'black'))
ggsave('../Results/Fig7.png', height = 7, width = 12)



#### Independent Counts ####
field_data %>%
  # filter(health == 'D') %>%
  group_by(asv_id) %>%
  summarise(present = sum(n_reads > 0),
            total = n(),
            prop = present / total)

field_data %>%
  # filter(health == 'D') %>%
  group_by(asv_id, present = n_reads > 0) %>%
  summarise(n = n(),
            .groups = 'drop') %>%
  pivot_wider(names_from = present, values_from = n) %>%
  column_to_rownames('asv_id') %T>%
  print %>%
  chisq.test()


#### Pathogen Counting ####
pathogen_counts <- pivot_wider(field_data,
                               names_from = asv_id, 
                               values_from = n_reads) %>%
  mutate(across(starts_with('ASV'), ~. > 0),
         combination = case_when(ASV8 & !ASV25 & !ASV38 ~ 'only8',
                                 !ASV8 & ASV25 & !ASV38 ~ 'only25',
                                 !ASV8 & !ASV25 & ASV38 ~ 'only38',
                                 
                                 ASV8 & ASV25 & !ASV38 ~ 'both_8_25',
                                 ASV8 & !ASV25 & ASV38 ~ 'both_8_38',
                                 !ASV8 & ASV25 & ASV38 ~ 'both_25_38',
                                 
                                 !ASV8 & !ASV25 & !ASV38 ~ 'none',
                                 ASV8 & ASV25 & ASV38 ~ 'all',
                                 
                                 TRUE ~ 'other')) %>%
  count(health, combination) %>%
  pivot_wider(names_from = combination,
              values_from = n,
              values_fill = 0) 


pathogen_counts %>%
  pivot_longer(cols = -health) %>%
  mutate(prop = value / sum(value), .by = name) %>%
  select(-value) %>%
  filter(health == 'D') %>%
  pivot_wider(values_from = prop)

select(pathogen_counts, -starts_with('both'))

column_to_rownames(pathogen_counts, 'health') %>%
  apply(1, sum)

select(pathogen_counts, health, all, none) %>%
  column_to_rownames('health') %T>%
  print %>%
  chisq.test(simulate.p.value = FALSE)

select(pathogen_counts, health, starts_with('only'), starts_with('both')) %>%
  pivot_longer(cols = -health,
               names_to = 'asv',
               values_to = 'number') %>%
  mutate(asv = if_else(asv == 'all', '25_8_38', asv)) %>%
  mutate(asv = str_remove(asv, 'only'),
         asv = str_remove(asv, 'both_'),
         asv = str_split(asv, '_')) %>%
  unnest(asv) %>%
  group_by(health, asv) %>%
  summarise(n = sum(number),
            .groups = 'drop') %>%
  pivot_wider(names_from = asv,
              values_from = n,
              names_prefix = 'asv') %>%
  column_to_rownames('health') %T>%
  print %>%
  chisq.test(simulate.p.value = FALSE)
  # chisq.posthoc.test()


#### Others ####
select(pathogen_counts, health, starts_with('only')) %>%
  column_to_rownames('health') %T>%
  print %>%
  chisq.test()

#only25
select(pathogen_counts, health, starts_with('only')) %>%
  mutate(other = only38 + only8, .keep = 'unused') %>%
  column_to_rownames('health') %T>%
  print %>%
  fisher.test()

#only38
select(pathogen_counts, health, starts_with('only')) %>%
  mutate(other = only25 + only8, .keep = 'unused') %>%
  column_to_rownames('health') %T>%
  print %>%
  fisher.test()

#only8
select(pathogen_counts, health, starts_with('only')) %>%
  mutate(other = only25 + only38, .keep = 'unused') %>%
  column_to_rownames('health') %T>%
  print %>%
  fisher.test()


select(pathogen_counts, health, all, none, starts_with('only')) %>%
  column_to_rownames('health') %T>%
  print %>%
  # fisher.test(simulate.p.value = FALSE)
  chisq.posthoc.test()


select(pathogen_counts, health, starts_with('both')) %>%
  column_to_rownames('health') %T>%
  print %>%
  fisher.test(simulate.p.value = FALSE)

## 25 vs 8 ##
tibble(asv_id = str_c('ASV', c(25,8)),
       disease = c(12, 21),
       healthy = c(7, 60)) %>%
  rowwise %>%
  mutate(broom::tidy(prop.test(x = disease, n = disease + healthy)),
         se = sqrt(estimate * (1 - estimate) / (disease + healthy))) %>%
  select(asv_id, estimate, se)

#### Experimental ####
library(emmeans)
pathogen_counts %>%
  pivot_longer(cols = -health) %>%
  pivot_wider(names_from = health, values_from = value) %>%
  mutate(ratio = (D + 1) / (H + 1),
         total = D + H,
         prop = D / total) %>%
  mutate(has_8 = name == 'all' | (str_detect(name, '8') & !name %in% c('both_25_38', 'only38')),
         has_25 = name == 'all' | str_detect(name, '25'),
         has_38 = name == 'all' | str_detect(name, '38')) %>%
  
  # summarise(across(c(D, H), sum), .by = has_38) %>% mutate(prop = D / (D + H))
  
  glm(cbind(D, H) ~ has_8 + has_25 + has_38, data = ., family = binomial) %>%
  emmeans(~has_8, type = 'response') %>%
  contrast('revpairwise', null = 143 / 270)
# anova(test = 'LRT')
# lm(ratio ~ has_8 + has_25 + has_38, data = .) %>%
# anova


tst <- pathogen_counts %>%
  pivot_longer(cols = -health) %>%
  pivot_wider(names_from = health, values_from = value) %>%
  mutate(ratio = (D + 1) / (H + 1),
         total = D + H,
         prop = D / total,
         name = fct_relevel(name, 'none')) %>%
  glm(cbind(D, H) ~ name, data = ., family = binomial)
car::Anova(tst)
summary(tst)


pathogen_counts %>%
  pivot_longer(cols = -health) %>%
  pivot_wider(names_from = health, values_from = value) %>%
  mutate(has_8 = name == 'all' | (str_detect(name, '8') & !name %in% c('both_25_38', 'only38')),
         has_25 = name == 'all' | str_detect(name, '25'),
         has_38 = name == 'all' | str_detect(name, '38')) %>%
  select(-name) %>%
  pivot_longer(cols = starts_with('has')) %>%
  filter(value) %>%
  group_by(name) %>%
  summarise(across(c(D, H), sum),
            .groups = 'rowwise') %>%
  mutate(broom::tidy(binom.test(x = D, n = D + H, p = 0.34)))


pathogen_counts %>%
  pivot_longer(cols = -health) %>%
  filter(str_detect(name, 'only')) %>%
  pivot_wider(names_from = health,
              values_from = value) %>%
  rowwise  %>%
  mutate(broom::tidy(binom.test(x = D, n = D + H, p = 143 / (143 + 270))))



pathogen_counts %>%
  pivot_longer(cols = -health) %>%
  pivot_wider(names_from = health, values_from = value) %>%
  mutate(has_8 = name == 'all' | (str_detect(name, '8') & !name %in% c('both_25_38', 'only38')),
         has_25 = name == 'all' | str_detect(name, '25'),
         has_38 = name == 'all' | str_detect(name, '38')) %>%
  select(-name) %>%
  pivot_longer(cols = starts_with('has')) %>%
  group_by(name, value) %>%
  summarise(across(c(D, H), sum)) %>%
  mutate(D = D / 143,
         H = H / 270)

