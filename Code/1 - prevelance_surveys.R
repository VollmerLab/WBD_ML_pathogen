
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

#### WBD Analysis ####
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
ggsave('../Results/Fig1_wbd_prevelance.png', height = 7, width = 7)

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
