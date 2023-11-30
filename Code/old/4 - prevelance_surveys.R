
library(multcomp)
library(tidyverse)
library(lubridate)
library(lmerTest)
library(emmeans)

#### Data ####
prevelance_data <- read_csv('C:/Users/jdsel/Dropbox/1_AC_field_abundprev/Acer5sites_HS_Added.csv', 
                            show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  select(-x1) %>%
  mutate(date = str_replace(date, '\\?', '15'),
         date = mdy(date),
         season = str_extract(time, 'S|W'),
         year = year(date)) %>%
  select(-time, -month_year, -timepoint1, -transect_number) %>%
  mutate(timepoint = factor(timepoint, ordered = TRUE),
         year = factor(year))

#### Map Raw Data ####
count(prevelance_data, season, year, site) %>%
  pivot_wider(names_from = 'site',
              values_from = 'n')

prevelance_data %>%
  ggplot(aes(x = meter, y = if_else(side == 'R', -1, 1), fill = as.character(wbd))) +
  geom_tile() +
  facet_grid(timepoint ~ site) +
  scale_fill_manual(values = c('0' = 'black', '1' = 'red')) +
  labs(x = 'meter',
       y = 'meter',
       fill = 'WBD Presence')
ggsave('../Results/WBD_transects.png', height = 5, width = 10)

prevelance_data %>%
  ggplot(aes(x = meter, y = if_else(side == 'R', -1, 1), fill = as.character(acerv))) +
  geom_tile() +
  facet_grid(timepoint ~ site) +
  scale_fill_manual(values = c('0' = 'black', '1' = 'red')) +
  labs(x = 'meter',
       y = 'meter',
       fill = 'Acropora Presence')
ggsave('../Results/Acropora_transects.png', height = 5, width = 10)


prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            .groups = 'drop') %>%
  
  mutate(date = ymd(str_c(year, if_else(season == 'S', 7, 1), 15, sep = '-'))) %>%
  rowwise %>%
  mutate(broom::tidy(prop.test(x = n_wbd, n = n_acerv))) %>%
  
  ggplot(aes(x = date, y = n_wbd / n_acerv)) +
  geom_line() +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~site)
  
#### Acerv Analysis ####
acerv_model <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            .groups = 'drop') %>%
  glm(cbind(n_acerv, total_meters - n_acerv) ~ timepoint * site, 
      family = 'binomial',
      data = .)

car::Anova(acerv_model)
  
emmeans(acerv_model, ~timepoint * site, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  left_join(count(prevelance_data, timepoint, date, site, season, year),
            by = c('timepoint', 'site')) %>%
  
  ggplot(aes(x = date, y = prob, colour = site)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_x_date(date_labels = '%b\n%y') +
  labs(x = NULL, 
       y = 'Acerv Abundance (%)',
       colour = NULL) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black', size = 12))
  

emmeans(acerv_model, ~site) %>%
  contrast('pairwise')

emmeans(acerv_model, ~site, type = 'response')

emmeans(acerv_model, ~timepoint) %>%
  contrast('poly')

emmeans(acerv_model, ~timepoint, type = 'response') %>% 
  broom::tidy(conf.int = TRUE) %>%
  left_join(count(prevelance_data, timepoint, date, site, season, year) %>%
              group_by(timepoint, season, year) %>%
              summarise(date = median(date),
                        .groups = 'drop'),
            by = c('timepoint')) %>%
  
  ggplot(aes(x = date, y = prob)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_x_date(date_labels = '%b\n%y') +
  labs(x = NULL, 
       y = 'Acerv Abundance (%)',
       colour = NULL) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black', size = 12))

#### WBD Analysis ####
full_model <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            .groups = 'drop') %>%
  glm(cbind(n_wbd, n_acerv - n_wbd) ~ timepoint * site, 
      family = 'binomial',
      data = .)

null_model <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            .groups = 'drop') %>%
  glm(cbind(n_wbd, n_acerv - n_wbd) ~ 1, 
      family = 'binomial',
      data = .)

1-logLik(full_model)/logLik(null_model)

count(prevelance_data, timepoint, season)
count(prevelance_data, timepoint, year)
count(prevelance_data, timepoint, date, site)

car::Anova(full_model, test.statistic = 'LR')

model_grid <- ref_grid(full_model)

add_grouping(model_grid, 'season', 'timepoint', 
             factor(c('S', 'W', 'S', 'W', 'S'))) %>%
  emmeans(~season) %>%
  contrast('pairwise')

add_grouping(model_grid, 'season', 'timepoint', 
             factor(c('S', 'W', 'S', 'W', 'S'))) %>%
  emmeans(~season | site, type = 'response') %>%
  contrast('pairwise')


add_grouping(model_grid, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year, type = 'response') %>%
  contrast('poly')

add_grouping(model_grid, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year, type = 'response') %>%
  broom::tidy() %>%
  ggplot(aes(x = year, y = prob)) +
  geom_point()

add_grouping(model_grid, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year | site, type = 'response') %>%
  contrast('poly')


emmeans(model_grid, ~timepoint, type = 'response') %>%
  cld(Letters = LETTERS) %>%
  broom::tidy() %>%
  left_join(count(prevelance_data, timepoint, date, site, season, year) %>%
              group_by(timepoint, season, year) %>%
              summarise(date = median(date),
                        .groups = 'drop'),
            by = c('timepoint')) %>%
  
  ggplot(aes(x = date, y = prob)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(x = NULL, 
       y = 'White Band Disease Prevelance (%)') +
  theme_classic()

emmeans(model_grid, ~timepoint * site, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  left_join(count(prevelance_data, timepoint, date, site, season, year),
            by = c('timepoint', 'site')) %>%
  
  ggplot(aes(x = date, y = prob, colour = site)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_x_date(date_labels = '%b\n%y') +
  labs(x = NULL, 
       y = 'White Band Disease Prevelance (%)',
       colour = NULL) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black', size = 12))
ggsave('../Results/WBD_prevelance_acrossSites.png', height = 7, width = 7)

#### Site as Random Effect ####
full_model_random <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            .groups = 'drop') %>%
  glmer(cbind(n_wbd, n_acerv - n_wbd) ~ timepoint + (1 + timepoint | site), 
      family = 'binomial',
      data = .)

full_model_random_noInt <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            .groups = 'drop') %>%
  glmer(cbind(n_wbd, n_acerv - n_wbd) ~ timepoint + (0 + timepoint | site), 
        family = 'binomial',
        data = .)

full_model_random_noSlope <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            .groups = 'drop') %>%
  glmer(cbind(n_wbd, n_acerv - n_wbd) ~ timepoint + (1 | site), 
        family = 'binomial',
        data = .)

anova(full_model_random_noInt, full_model_random)
anova(full_model_random_noSlope, full_model_random)

model_data <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            date = median(date),
            .groups = 'drop')
full_model_random <-  glmer(cbind(n_wbd, n_acerv - n_wbd) ~ 
                              timepoint + (1 | site), 
                            family = 'binomial',
                            data = model_data)

car::Anova(full_model_random, method = 'LRT')
summary(full_model_random)

significance_groupings <- glht(full_model_random, 
                               linfct = mcp(timepoint = "Tukey")) %>%
  summary(test = adjusted("Westfall")) %>% 
  cld() %$%
  mcletters %$% 
  Letters %>% 
  enframe(name = 'timepoint', value = '.group') %>%
  mutate(.group = str_to_upper(.group))

 
emmeans(full_model_random, ~timepoint, type = 'response') %>%
  # cld(Letters = LETTERS, adjust = 'none') %>%
  broom::tidy(conf.int = TRUE) %>%
  left_join(significance_groupings, by = 'timepoint') %>%
  mutate(.group = str_trim(.group)) %>%
  left_join(count(prevelance_data, timepoint, date, site, season, year) %>%
              group_by(timepoint, season, year) %>%
              summarise(date = median(date),
                        .groups = 'drop'),
            by = c('timepoint')) %>%
  
  ggplot(aes(x = date, y = prob)) +
  # geom_point(data = model_data, size = 0.5,
  #            aes(y = n_wbd / n_acerv)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_text(aes(y = Inf, label = .group), vjust = 1.5) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(x = NULL, 
       y = 'White Band Disease Prevelance (%)') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        panel.background = element_rect(colour = 'black'))
ggsave('../Manuscript/Fig1.png', height = 7, width = 7)

model_grid_random <- ref_grid(full_model_random)

add_grouping(model_grid_random, 'low_high', 'timepoint', 
             factor(c('low', 'high', 'other', 'low', 'other'))) %>%
  emmeans(~low_high, type = 'response')


add_grouping(model_grid_random, 'season', 'timepoint', 
             factor(c('S', 'W', 'S', 'W', 'S'))) %>%
  emmeans(~season) %>%
  contrast('pairwise')

add_grouping(model_grid_random, 'season', 'timepoint', 
             factor(c('S', 'W', 'S', 'W', 'S'))) %>%
  emmeans(~season, type = 'response') 

add_grouping(model_grid_random, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year) %>%
  contrast('pairwise')


add_grouping(model_grid_random, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year, type = 'link') %>%
  # contrast('pairwise', adjust = 'mvt')
  pairs() %>%
  as.glht() %>%
  summary(test = adjusted("Westfall"))
  
add_grouping(model_grid_random, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year, type = 'response')

#### Model WBD ####
library(rstanarm)
full_model3 <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            .groups = 'drop') %>%
  stan_glmer(cbind(n_wbd, n_acerv - n_wbd) ~ timepoint + (1 + timepoint | site), 
        family = 'binomial',
        data = .,
        chains = 4, 
        cores = 4)


emmeans(full_model3, ~timepoint, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  left_join(count(prevelance_data, timepoint, date, site, season, year) %>%
              group_by(timepoint, season, year) %>%
              summarise(date = median(date),
                        .groups = 'drop'),
            by = c('timepoint')) %>%
  
  ggplot(aes(x = date, y = prob)) +
  geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(x = NULL, 
       y = 'White Band Disease Prevelance (%)') +
  theme_classic()

model_grid3 <- ref_grid(full_model3)

add_grouping(model_grid3, 'season', 'timepoint', 
             factor(c('S', 'W', 'S', 'W', 'S'))) %>%
  emmeans(~season) %>%
  contrast('pairwise')

add_grouping(model_grid3, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year) %>%
  contrast('pairwise')

add_grouping(model_grid3, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  mutate(grouping = c('A', 'B', 'A'))

