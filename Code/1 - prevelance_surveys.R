
library(multcomp)
library(tidyverse)
library(lubridate)
library(lmerTest)
library(emmeans)
library(patchwork)
library(ggtext)

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
  filter(date >= min(prevelance_data$date) - months(3),
         date <= ymd('2017-08-01')) %>%
  # sample_n(500) %>%
  ggplot(aes(x = date, y = water_temp)) +
  geom_hline(yintercept = 30, linetype = 'dashed') +
  geom_path() +
  scale_x_date(breaks = ymd(c('2015-07-01', 
                              '2016-01-01', '2016-07-01', 
                              '2017-01-01', '2017-07-01')), 
               date_labels = '%b %Y',
               limits = ymd(c('2015-03-23', '2017-08-01'))) +
  labs(x = NULL, 
       y = 'Mean Daily Temperature (C)') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        panel.background = element_rect(colour = 'black'))

dhw_plot <- bocas_temp %>%
  filter(date >= min(prevelance_data$date) - months(3),
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

#### Cyclic Temperature Plot ####
cyclical_temp_plot <- bocas_temp %>%
  mutate(month = month(date),
         year = year(date),
         day = day(date)) %>%
  filter(year != 2005) %>%
  ggplot(aes(x = ymd(str_c('2000', month, day, sep = '-')), y = water_temp, group = year)) +
  geom_hline(yintercept = 30, linetype = 'dashed') +
  geom_line() +
  # geom_vline(data = mutate(the_dates, fake_date = ymd(str_c('2000',  month(date),day(date),sep='-'))), 
  #            aes(xintercept = fake_date)) +
  # geom_vline(xintercept = ymd(c('2000-04-01', '2000-07-01')), colour = 'blue') +
  # geom_vline(xintercept = ymd(c('2000-09-01', '2000-12-01')), colour = 'red') +
  scale_x_date(date_labels = '%b') +
  labs(x = NULL,
       y = 'Mean Daily Temperature (°C)') +
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
            .groups = 'drop') %>%
  left_join(bocas_temp,
            by = c('date'))

#get days above 30 in preceeding heat wave (april-june or sept-nov) & max temp & number of days since peak
model_data_2 <- bocas_temp %>%
  mutate(month = month(date),
         year = year(date),
         day = day(date)) %>%
  mutate(hw = case_when(month %in% c(4:6) ~ 'early',
                        month %in% c(9:11) ~ 'late',
                        TRUE ~ 'non')) %>%
  group_by(year, hw) %>%
  summarise(max_temp = max(water_temp),
            mean_temp = mean(water_temp),
            hot_day = sum(water_temp >= 30),
            n_day = n(),
            peak_date = date[water_temp == max_temp],
            .groups = 'drop') %>%
  filter(hw != 'non') %>%
  mutate(year = if_else(hw == 'late', year + 1, year)) %>%
  right_join(mutate(model_data, 
                    hw = if_else(season == 'W', 'late', 'early'),
                    year_f = year,
                    year = as.numeric(as.character(year))),
             by = c('year', 'hw')) %>%
  mutate(days_since = as.numeric(date - peak_date),
         timepoint = as.character(timepoint))

#### Analysis of Acerv Density ####
acerv_model_random <-  glmer(cbind(n_acerv, total_meters - n_acerv) ~ 
                              timepoint + (1 | site), 
                            family = 'binomial',
                            data = model_data_2)
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
  scale_x_date(breaks = ymd(c('2015-07-01', 
                              '2016-01-01', '2016-07-01', 
                              '2017-01-01', '2017-07-01')), 
               date_labels = '%b\n%Y',
               limits = ymd(c('2015-03-23', '2017-08-01'))) +
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

#### WBD Analysis ####
full_model_random <-  glmer(cbind(n_wbd, n_acerv - n_wbd) ~ 
                              timepoint + (1 | site), 
                            family = 'binomial',
                            data = model_data_2)

car::Anova(full_model_random, method = 'LRT')
summary(full_model_random)

## Significance of site
no_site_model <- glmmTMB(cbind(n_wbd, n_acerv - n_wbd) ~ 
                           timepoint, 
                         family = 'binomial',
                         data = model_data)

full_model_random <-  glmmTMB(cbind(n_wbd, n_acerv - n_wbd) ~ 
                                timepoint + (1 | site), 
                              family = 'binomial',
                              data = model_data)

anova(no_site_model, full_model_random)
##

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
               date_labels = '%b\n%Y',
               limits = ymd(c('2015-03-23', '2017-08-01'))) +
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


#### Associations with Temperature/WBD ####
library(nlme)
#model <- temp_association$model[[1]]
process_gls <- function(model){
  tst <- anova(model)
  as_tibble(tst, rownames = 'effect') %>%
    filter(effect == 'temp_metric') %>%
    select(-effect) %>%
    mutate(denDF = str_extract(attr(tst, 'label'), '[0-9]([0-9\\.]+)?') %>% 
             as.numeric(),
           .after = numDF) %>%
    janitor::clean_names()
}

temp_association <- emmeans(full_model_random, ~timepoint, 
                            type = 'link') %>%
  as_tibble() %>%
  left_join(select(model_data_2,  timepoint, hot_day, max_temp, mean_temp) %>%
              distinct,
            by = 'timepoint') %>%
  mutate(var = SE^2) %>%
  pivot_longer(cols = c(hot_day, max_temp, mean_temp),
               names_to = 'metric',
               values_to = 'temp_metric') %>%
  nest_by(metric) %>%
  mutate(model = list(gls(emmean ~ temp_metric, 
                          weight=varIdent(var), 
                          data = data)),
         process_gls(model))

temp_association

temp_wbd_plot <- ref_grid(temp_association$model[[1]], at = list(temp_metric = seq(0, 90, length.out = 50)),
         data = temp_association$data[[1]]) %>%
  update(tran = binomial(link = 'logit')) %>%
  emmeans(~temp_metric, 
          type = 'response') %>%
  as_tibble() %>%
  rename(hot_day = temp_metric) %>%
  ggplot(aes(x = hot_day, y = response, ymin = lower.CL, ymax = upper.CL)) +
  geom_ribbon(alpha = 0.5) +
  geom_pointrange(data = emmeans(full_model_random, ~timepoint, 
                                 type = 'response') %>%
                    as_tibble() %>%
                    left_join(select(model_data_2,  timepoint, hot_day) %>%
                                distinct,
                              by = 'timepoint'),
                  inherit.aes = FALSE,
                  aes(x = hot_day, y = prob, ymin = asymp.LCL, ymax = asymp.UCL)) +
  geom_line() +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(x = 'Days ≥ 30°C', 
       y = 'White Band Disease Prevelance (%)') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        panel.background = element_rect(colour = 'black'))



#### Make Merged Plot ####


(temp_plot + labs(y = 'Mean Daily Temperature (°C)') + 
   theme(axis.text.x = element_blank())) / 
  # (dhw_plot + theme(axis.text.x = element_blank())) /
  (prevelance_plot + labs(y = 'WBD Prevelance (%)') + 
     theme(axis.text.x = element_blank())) / 
  (abundance_plot + labs(y = '<i>A. cervicornis</i> Abundance (%)') +
     theme(axis.text.x = element_text(size = 12),
           axis.title.y = element_markdown())) &
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 16),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10, 
                                   colour = 'black')) 
# ggsave('../Results/Fig1_temp_acer_wbd.png', height = 10, width = 5)





(((temp_plot + labs(y = 'Mean Daily Temperature (°C)') + 
    theme(axis.text.x = element_blank())) / 
  # (dhw_plot + theme(axis.text.x = element_blank())) /
  (prevelance_plot + labs(y = 'WBD Prevelance (%)') + 
     theme(axis.text.x = element_blank())) / 
  (abundance_plot + labs(y = '<i>A. cervicornis</i> Abundance (%)') +
     theme(axis.text.x = element_text(size = 12),
           axis.title.y = element_markdown()))) | 
  ((cyclical_temp_plot +
      theme(axis.text.x = element_text(size = 12))) /
     (temp_wbd_plot + labs(y = 'WBD Prevelance (%)') +
        theme(axis.text.x = element_text(size = 12))))) &
  
  plot_annotation(tag_levels = list(c('A', 'C', 'D', 'B', 'E'))) &
  theme(plot.tag = element_text(face = 'bold', size = 16),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10, 
                                   colour = 'black')) 
ggsave('../Results/Fig1_temp_acer_wbd.png', height = 10, width = 7)
