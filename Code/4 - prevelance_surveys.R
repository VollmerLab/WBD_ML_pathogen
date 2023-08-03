
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
  

#### Analysis ####
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

car::Anova(full_model)

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
full_model2 <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            .groups = 'drop') %>%
  glmer(cbind(n_wbd, n_acerv - n_wbd) ~ timepoint + (1 | site), 
      family = 'binomial',
      data = .)

AIC(full_model2, full_model)


emmeans(full_model2, ~timepoint, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
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

summary(full_model2)


model_grid2 <- ref_grid(full_model2)

add_grouping(model_grid2, 'season', 'timepoint', 
             factor(c('S', 'W', 'S', 'W', 'S'))) %>%
  emmeans(~season) %>%
  contrast('pairwise')

add_grouping(model_grid2, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year) %>%
  contrast('pairwise')



add_grouping(model_grid2, 'year', 'timepoint', 
             factor(c('2015', '2016', '2016', '2017', '2017'))) %>%
  emmeans(~year, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  mutate(grouping = c('A', 'B', 'A'))
