
the_data <- select(filtered_wide_vals, -median_rank) %>%
  pivot_longer(cols = starts_with('base'),
               names_to = 'model',
               values_to = 'rank') %>%
  mutate(rank = floor(rank))

library(glmmTMB)
library(emmeans)

filter(the_data, asv_id == 'ASV191')

tst <- glmmTMB(rank ~ asv_id + (1 | model), 
               dispformula = ~asv_id,
               family = genpois(link = 'log'),
               data = the_data)
summary(tst)

emmeans(tst, ~asv_id, type = 'response') %>%
  as_tibble() %>%
  arrange(response) %>%
  mutate(asv_id = fct_reorder(asv_id, response)) %>%
  ggplot(aes(x = response, xmin = asymp.LCL, xmax = asymp.UCL,
             y = asv_id, colour = asv_id %in% keep$asv_id)) +
  geom_pointrange()


keep <- emmeans(tst, ~asv_id, type = 'response') %>%
  contrast(method = 'eff', side = "<") %>%
  as_tibble() %>%
  filter(p.value < 0.05) %>%
  mutate(asv_id = str_remove(contrast, ' effect')) %>%
  arrange(ratio) 

direction <- 'less'

asv_ranks %>%
  filter(!is.na(p_adjust)) %>%
  pivot_longer(cols = -any_of(c('asv_id', 'median_rank', 'median_shap', 'statistic', 
                                'p.value', 'p_adjust', 'method', 'alternative'))) %>%
  
  
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(higher_taxonomy = str_c(family, genus, asv_id, sep = '; ') %>%
           str_remove('; NA')) %>%
  
  mutate(higher_taxonomy = fct_reorder(higher_taxonomy, value, 
                                       .desc = if_else(direction == 'less', TRUE, FALSE))) %>%
  
  # filter(median_rank <= 50) %>%
  
  ggplot(aes(x = value, y = higher_taxonomy, colour = asv_id %in% keep$asv_id)) +
  # geom_point(size = 0.5) +
  # geom_boxplot(outlier.shape = NA) +
  # stat_summary(fun.data = median_hilow, show.legend = FALSE) +
  stat_summary(fun = median, fun.min = min, fun.max = max, show.legend = FALSE) +
  # stat_summary(fun.data = mean_se, show.legend = FALSE) +
  # stat_summary(fun.data = mean_cl_boot, show.legend = FALSE) +
  scale_colour_manual(values = c('FALSE' = 'black', 'TRUE' = 'red')) +
  labs(x = if_else(direction == 'less', 'ASV Rank', 'SHAP Importance'),
       y = NULL) +
  theme_classic() 
