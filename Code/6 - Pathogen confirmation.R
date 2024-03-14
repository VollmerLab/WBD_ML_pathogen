library(tidyverse)
library(chisq.posthoc.test)

pathogens <- str_c('ASV', c(8, 25, 38))


#### Data ####
field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  filter(asv_id %in% pathogens) %>%
  select(asv_id, sample_id, n_reads, health)

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
  chisq.test(simulate.p.value = FALSE)

select(pathogen_counts, health, all, none, starts_with('only')) %>%
  column_to_rownames('health') %T>%
  print %>%
  # fisher.test(simulate.p.value = FALSE)
  chisq.posthoc.test()


select(pathogen_counts, health, starts_with('both')) %>%
  column_to_rownames('health') %T>%
  print %>%
  fisher.test(simulate.p.value = FALSE)


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

