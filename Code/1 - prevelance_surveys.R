
library(multcomp)
library(tidyverse)
library(lubridate)
library(lmerTest)
library(emmeans)
library(patchwork)

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

# water temp data downloaded from: @1/9/24
#https://biogeodb.stri.si.edu/physical_monitoring/research/bocas#download
# Degree Heating Weeks
#https://www.pacioos.hawaii.edu/voyager/info/coral_bleaching_degree_heating_week.html
n_weeks <- 12
bleach_threshold <- 30

bocas_temp <- read_csv('../Data/bocas_temp/bocas_tower_wt_elect.csv', 
                       show_col_types = FALSE) %>%
  filter(chk_note == 'good') %>%
  mutate(datetime = mdy_hms(datetime),
         date = mdy(date)) %>%
  select(-raw:-chk_fail) %>%
  rename(water_temp = wt) %>%
  group_by(date) %>%
  summarise(water_temp = mean(water_temp)) %>%
  mutate(dhw = zoo::rollapply(water_temp, width = 7 * n_weeks, 
                              function(x) sum(x > bleach_threshold), 
                              align = "right", fill = 0),
         dhw = dhw / 7)

#### Temperature Plot ####
# temp_plot <- bocas_temp %>%
#   # sample_n(500) %>%
#   ggplot(aes(x = datetime, y = water_temp)) +
#   geom_hline(yintercept = 30, linetype = 'dashed') +
#   geom_path() +
#   scale_x_datetime(breaks = ymd_h(c('2015-01-01 0', '2015-07-01 0', 
#                                     '2016-01-01 0', '2016-07-01 0', 
#                                     '2017-01-01 0', '2017-07-01 0')), 
#                    date_labels = '%b %Y') 

temp_plot <- bocas_temp %>%
  filter(date >= min(prevelance_data$date) - months(6),
         date <= max(prevelance_data$date)) %>%
  # sample_n(500) %>%
  ggplot(aes(x = date, y = water_temp)) +
  geom_hline(yintercept = 30, linetype = 'dashed') +
  geom_path() +
  scale_x_date(breaks = ymd(c('2015-01-01', '2015-07-01', 
                              '2016-01-01', '2016-07-01', 
                              '2017-01-01', '2017-07-01')), 
               date_labels = '%b %Y') +
  labs(x = NULL, 
       y = 'Mean Daily Temperature (C)') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        panel.background = element_rect(colour = 'black'))

dhw_plot <- bocas_temp %>%
  filter(date >= min(prevelance_data$date) - months(6),
         date <= max(prevelance_data$date)) %>%
  # sample_n(500) %>%
  ggplot(aes(x = date, y = dhw)) +
  geom_hline(yintercept = 4, linetype = 'dashed') + #Sensitive Species Bleaching
  geom_hline(yintercept = 8, linetype = 'dotted') + #Widespread Bleaching
  geom_path() +
  scale_x_date(breaks = ymd(c('2015-01-01', '2015-07-01', 
                              '2016-01-01', '2016-07-01', 
                              '2017-01-01', '2017-07-01')), 
               date_labels = '%b %Y') +
  labs(x = NULL, 
       y = 'Degree Heating Weeks (C-weeks)') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        panel.background = element_rect(colour = 'black'))

#### Map Raw Data ####
count(prevelance_data, season, year, site) %>%
  pivot_wider(names_from = 'site',
              values_from = 'n')

model_data <- prevelance_data %>%
  group_by(site, timepoint, season, year) %>%
  summarise(n_acerv = sum(acerv),
            n_wbd = sum(wbd),
            total_meters = n(),
            date = median(date),
            .groups = 'drop')

#### Analysis of Acerv Density ####
acerv_model_random <-  glmer(cbind(n_acerv, total_meters - n_acerv) ~ 
                              timepoint + (1 | site), 
                            family = 'binomial',
                            data = model_data)
car::Anova(acerv_model_random, method = 'LRT')
summary(acerv_model_random)

acerv_significance_groupings <- glht(acerv_model_random, 
                               linfct = mcp(timepoint = "Tukey")) %>%
  summary(test = adjusted("Westfall")) %>% 
  cld() %$%
  mcletters %$% 
  Letters %>% 
  enframe(name = 'timepoint', value = '.group') %>%
  mutate(.group = str_to_upper(.group))


abundance_plot <- emmeans(acerv_model_random, ~timepoint, 
                          type = 'response') %>%
  # cld(Letters = LETTERS, adjust = 'none') %>%
  broom::tidy(conf.int = TRUE) %>%
  left_join(acerv_significance_groupings, by = 'timepoint') %>%
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
  scale_x_date(breaks = ymd(c('2015-01-01', '2015-07-01', 
                              '2016-01-01', '2016-07-01', 
                              '2017-01-01', '2017-07-01')), 
               date_labels = '%b %Y',
               limits = ymd(c('2015-01-01', '2017-08-01'))) +
  labs(x = NULL, 
       y = 'Acropora cervicornis prevelance (%)') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        panel.background = element_rect(colour = 'black'))

model_grid_acer <- ref_grid(acerv_model_random)

add_grouping(model_grid_acer, 'low_high', 'timepoint', 
             factor(c('high', 'other', 'low', 'low', 'low'))) %>%
  emmeans(~low_high, type = 'response')


add_grouping(model_grid_random, 'season', 'timepoint', 
             factor(c('S', 'W', 'S', 'W', 'S'))) %>%
  emmeans(~season) %>%
  contrast('pairwise')

#### WBD Analysis ####
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


prevelance_plot <- emmeans(full_model_random, ~timepoint, 
                           type = 'response') %>%
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
  scale_x_date(breaks = ymd(c('2015-01-01', '2015-07-01', 
                              '2016-01-01', '2016-07-01', 
                              '2017-01-01', '2017-07-01')), 
               date_labels = '%b %Y',
               limits = ymd(c('2015-01-01', '2017-08-01'))) +
  labs(x = NULL, 
       y = 'White Band Disease Prevelance (%)') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        panel.background = element_rect(colour = 'black'))
# ggsave('../Results/Fig1_wbd_prevelance.png', height = 7, width = 7)

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

#### Make Merged Plot ####
(temp_plot + theme(axis.text.x = element_blank())) / 
  (dhw_plot + theme(axis.text.x = element_blank())) /
  (abundance_plot + theme(axis.text.x = element_blank())) / 
  prevelance_plot
