# cpm ~ (control-healthyEnd, disease-healthyEnd, disease-diseaseEnd) * dosed/not 

#### Libraries ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
# library(rstanarm)
library(lubridate)
library(emmeans)
library(broom)
library(patchwork)
library(multidplyr)
library(ComplexUpset)
library(cowplot)

#### N ASVs ####
alpha <- 0.05


cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('lmerTest', 'emmeans', 'tibble', 'nlme',
                           'dplyr', 'tidyr', 'stringr'))

#### Data ####
field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  mutate(timepoint = str_c(season, year, sep = '_'),
         across(domain:genus, str_replace_na)) %>%
  nest_by(across(c(domain:genus, asv_id)))

shap_importance <- read_csv('../Results/asv_importance.csv.gz',
                            show_col_types = FALSE) %>%
  filter(p_adjust < alpha) %>%
  mutate(asv_id = fct_reorder(asv_id, median_rank)) %>%
  select(asv_id, median_rank) %>%
  distinct %>%
  arrange(asv_id)

ml_classifications <- read_csv('../Results/Table45_asv_table_supplement.csv',
         show_col_types = FALSE) %>%
  select(ID, likely_type) %>%
  mutate(type = case_when(is.na(likely_type) ~ "ML ID'd",
                          likely_type == 'Commensalist' ~ 'Field Consistent',
                          TRUE ~ likely_type)) %>%
  select(-likely_type) %>%
  rename(asv_id = ID)

#### All Random Model ####
random_models <- field_data %>%
  partition(cluster) %>%
  mutate(model = list(lmer(log2_cpm_norm ~ health + (1 | timepoint) + (1 | site),
                           data = data)),
         anova_table = list(anova(model, ddf = "Kenward-Roger")),
         d_h = as_tibble(emmeans(model, 'pairwise'~health)$contrast)$estimate, #Diseased - Healthy
         nDF = anova_table$NumDF,
         dDF = anova_table$DenDF,
         fvalue = anova_table$`F value`,
         pvalue = anova_table$`Pr(>F)`) %>%
  collect %>%
  ungroup


#### Get Significance ####
random_models %>%
  mutate(padjust = p.adjust(pvalue, method = 'fdr')) %>%
  filter(padjust < alpha) %>%
  group_by(more_disease = d_h > 0) %>%
  summarise(n = n(),
            across(c(domain:asv_id), n_distinct))

random_models %>%
  mutate(padjust = p.adjust(pvalue, method = 'fdr')) %>%
  filter(padjust < alpha) %>%
  group_by(more_disease = d_h > 0, across(domain:genus)) %>%
  summarise(n = n()) %>%
  ungroup %>%
  mutate(higher_taxa = str_c(class, order, family, genus, sep = '; '),
         n = if_else(more_disease, n, -1 * n),
         higher_taxa = fct_reorder(higher_taxa, abs(n), .fun = sum)) %>%
  ggplot(aes(y = higher_taxa, x = n, fill = more_disease)) +
  geom_col() 
  


random_models %>%
  mutate(padjust = p.adjust(pvalue, method = 'fdr')) %>%
  filter(padjust < alpha) %>%
  group_by(more_disease = d_h > 0, across(domain:genus)) %>%
  summarise(n = n()) %>%
  ungroup %>%
  mutate(higher_taxa = str_c(class, order, family, genus, sep = '; '),
         higher_taxa = fct_reorder(higher_taxa, abs(n), .fun = sum)) %>%
  ggplot(aes(y = higher_taxa, x = n, fill = more_disease)) +
  geom_linerange(aes(xmin = 0, xmax = n), position = position_dodge(0.5),
                 linetype = 'dashed') +
  geom_point(shape = 21, position = position_dodge(0.5), size = 5) +
  scale_x_continuous(trans = scales::log10_trans())

select(random_models, -where(is.list)) %>%
  mutate(fdr = p.adjust(pvalue, method = 'fdr'),
         fdr = fdr < alpha,
         positive = d_h > 0) %>%
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        annotations = list(
          'Order' = ggplot(mapping = aes(fill = order)) +
            geom_bar(stat = 'count', position = 'fill') +
            scale_y_continuous(labels=scales::percent_format())
        )
  )


#### Fixed Model ####
fixed_models <- field_data %>%
  partition(cluster) %>%
  mutate(model = list(gls(log2_cpm_norm ~ health * timepoint * site,
                           data = data)),
         
         anova(model, type = 'marginal') %>%
           as_tibble(rownames = 'term') %>%
           select(term, `p-value`) %>%
           filter(term != '(Intercept)') %>%
           pivot_wider(names_from = term, values_from = `p-value`) %>%
           rename_with(~str_replace_all(., ':', 'X'))) %>%
  collect %>%
  ungroup


## Upset Plot ##


?upset

select(fixed_models, -where(is.list)) %>%
  mutate(across(health:healthXtimepointXsite, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepointXsite, ~. < alpha)) %>%
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        annotations = list(
          'Order' = ggplot(mapping = aes(fill = order)) +
            geom_bar(stat = 'count', position = 'fill') +
            scale_y_continuous(labels=scales::percent_format())
        )
  )

select(fixed_models, -where(is.list)) %>%
  mutate(across(health:healthXtimepointXsite, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepointXsite, ~. < alpha)) %>%
  filter(if_any(contains('health'))) 


#### Match Consistency Model - fixed timepoint, random site ####
fixed_models <- field_data %>%
  partition(cluster) %>%
  mutate(model = list(lmer(log2_cpm_norm ~ health * timepoint + (1 | site),
                           data = data)),
         
         anova(model, ddf = "Kenward-Roger") %>%
           as_tibble(rownames = 'term') %>%
           select(term, `Pr(>F)`) %>%
           pivot_wider(names_from = term, values_from = `Pr(>F)`) %>%
           rename_with(~str_replace_all(., ':', 'X')),
         
         up_disease = fixef(model)[['healthH']] < 0) %>%
  collect %>%
  ungroup



colour_options <- read_rds('../intermediate_files/asv_colors.rds')
microbe_colors <- set_names(colour_options$color_palette$hex,
          colour_options$color_palette$group)
levels(colour_options$asv_clumping$Top_order)

base_plot <- select(fixed_models, -where(is.list)) %>%
  left_join(colour_options$asv_clumping,
            by = c('order', 'genus')) %>%
  mutate(group = case_when(!is.na(group) ~ group,
                           is.na(group) & order == 'Oceanospirillales' ~ 'Oceanospirillales-Other',
                           is.na(group) & order == 'Rhodobacterales' ~ 'Rhodobacterales-Other',
                           is.na(group) & order == 'Alteromonadales' ~ 'Alteromonadales-Other',
                           is.na(group) & order == 'Thiotrichales' ~ 'Thiotrichales-Other',
                           is.na(group) & order == 'Vibrionales' ~ 'Vibrionales-Other',
                           TRUE ~ 'Other-Other'),
         group = factor(group, levels = levels(colour_options$asv_clumping$group))) %>%
  # filter(!is.na(group)) %>%
  select(-Top_order, -Top_genus) %>%
  rename(the_colour = group) %>%
  mutate(across(health:healthXtimepoint, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepoint, ~. < alpha),
         health_d = health & up_disease,
         health_h = health & !up_disease) %>%
  rename('Health State' = health,
         'Health State (Diseased)' = health_d,
         'Health State (Healthy)' = health_h,
         'Sampling Time' = timepoint,
         'Health State x Sampling Time' = healthXtimepoint) %>%
  select(-up_disease, -`Health State`) %>%
  mutate(top_asv = if_else(asv_id %in% shap_importance$asv_id, 'yes', 'no')) %>%
  # mutate(top_asv = case_when(asv_id %in% c('ASV25', 'ASV8', 'ASV38') ~ 'Pathogen',
  #                            asv_id %in% c('ASV26', 'ASV30', 'ASV361', 'ASV51') ~ 'Opportunist',
  #                            TRUE ~ top_asv)) %>%
  left_join(ml_classifications,
            by = 'asv_id') %>%
  
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        base_annotations = list(
          'Intersection size' = intersection_size(
            counts = TRUE,
            mapping = aes(fill = type)
          ) + 
            # scale_fill_manual(values = c('yes' = 'red', 'no' = '#595959', 
            #                              'Pathogen' = 'orange', 'Opportunist' = 'green')) +
            theme(legend.position = 'left')
        ),
        annotations = list(
          'Order' = ggplot(mapping = aes(fill = the_colour)) +
            geom_bar(stat = 'count', position = 'fill', show.legend = FALSE) +
            scale_y_continuous(labels=scales::percent_format()) +
            scale_fill_manual(values = microbe_colors) +
            labs(y = 'Abundance')
        ),
        name = NULL
  )

plot_grid(base_plot, colour_options$legend,  rel_widths = c(1, .25))
ggsave('../Results/Fig5_traditional_complexUpset.png', height = 12, width = 10, bg = 'white')


select(fixed_models, -where(is.list)) %>%
  mutate(across(health:healthXtimepoint, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepoint, ~. < alpha)) %>%
  filter(if_any(contains('health'))) 
