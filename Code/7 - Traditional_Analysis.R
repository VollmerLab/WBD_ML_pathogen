#CHECK OTHER FILL
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



#### Full Taxonomy Genera ####
n_genera <- ungroup(field_data) %>%
  count(genus)

#### Make Plot ####
colour_options <- read_rds('../intermediate_files/asv_colors.rds')
microbe_colors <- set_names(colour_options$color_palette$hex,
          colour_options$color_palette$group)
levels(colour_options$asv_clumping$Top_order)

base_plot_data <- select(fixed_models, -where(is.list)) %>%
  left_join(mutate(colour_options$asv_clumping,
                   genus = str_remove_all(genus, '\\<i\\>|\\</i\\>')),
            by = c('order', 'genus')) %>%
  mutate(group = case_when(!is.na(group) ~ group,
                           is.na(group) & order == 'Oceanospirillales' ~ 'Oceanospirillales-Other',
                           is.na(group) & order == 'Vibrionales' ~ 'Vibrionales-Other',
                           is.na(group) & order == 'Alteromonadales' ~ 'Alteromonadales-Other',
                           is.na(group) & order == 'Verrucomicrobiales' ~ 'Verrucomicrobiales-Other',
                           is.na(group) & order == 'Campylobacterales' ~ 'Campylobacterales-Other',
                           TRUE ~ 'Other-Other'),
         group = factor(group, levels = levels(colour_options$asv_clumping$group))) %>%
  # filter(!is.na(group)) %>%
  select(-Top_order, -Top_genus) %>%
  rename(the_colour = group) %>%
  mutate(across(health:healthXtimepoint, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepoint, ~. < alpha),
         health_d = health & up_disease,
         health_h = health & !up_disease) %>%
  rename('Disease State' = health,
         'Disease Associated' = health_d,
         'Healthy Associated' = health_h,
         'Time' = timepoint,
         'Disease State x Time' = healthXtimepoint) %>%
  select(-up_disease, -`Disease State`) %>%
  mutate(top_asv = if_else(asv_id %in% shap_importance$asv_id, 'yes', 'no')) %>%
  # mutate(top_asv = case_when(asv_id %in% c('ASV25', 'ASV8', 'ASV38') ~ 'Pathogen',
  #                            asv_id %in% c('ASV26', 'ASV30', 'ASV361', 'ASV51') ~ 'Opportunist',
  #                            TRUE ~ top_asv)) %>%
  left_join(ml_classifications,
            by = 'asv_id') %>% 
  mutate(type = if_else(is.na(type), 'no', 'yes')) %>%
  # mutate(type = case_when(is.na(type) ~ 'no',
  #                         type == 'Pathogen' ~ 'Pathogen',
  #                         TRUE ~ 'yes')) %>%
  identity()

base_plot_data %>%
  filter(`Healthy Associated`) %>%
  count(genus) %>%
  filter(genus != 'NA')


base_plot <- base_plot_data %>%
  # filter(!if_all(where(is.logical), ~!.)) %>%
  mutate(Time = Time | `Disease State x Time`, .keep = 'unused') %>%
  # filter(`Disease Associated` | `Healthy Associated`) %>%
  rename(Disease = `Disease Associated`,
         Healthy = `Healthy Associated`) %>%
  # relocate(Time, .after = `Disease State x Time`) %>%
  # select(-Time) %>%
  
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        base_annotations = list(
          'Intersection size' = intersection_size(
            counts = TRUE,
            mapping = aes(fill = type)
          ) + 
            # scale_fill_manual(values = c('yes' = 'red', 'no' = '#595959', 
            #                              'Pathogen' = 'orange', 'Opportunist' = 'green')) +
            # scale_fill_manual(values = c('yes' = 'gray', 'no' = '#595959', 'Pathogen' = 'red')) +
            scale_fill_manual(values = c('no' = 'gray65', 'yes' = 'gray25')) +
            scale_y_continuous(labels=scales::comma_format(), 
                               expand = expansion(mult = c(0.01, 0.05))) +
            labs(tag = 'B') +
            theme_classic() +
            theme(legend.position = 'none',
                  panel.background = element_rect(colour = 'black'),
                  # axis.ticks.y = element_line(colour = 'black'),
                  # axis.minor.ticks.y.left = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title = element_text(colour = 'black', size = 14),
                  axis.text = element_text(colour = 'black', size = 10),
                  plot.tag = element_text(colour = 'black', size = 18, 
                                            face = 'bold', hjust = 0))
        ),
        annotations = list(
          'Order' = ggplot(mapping = aes(fill = the_colour)) +
            geom_bar(fill = 'white', colour = 'black', stat = 'count', 
                     position = 'fill', show.legend = FALSE) +
            geom_bar(stat = 'count', position = 'fill', show.legend = FALSE) +
            scale_y_continuous(labels=scales::percent_format(), 
                               expand = expansion(mult = c(0.01, 0.05))) +
            scale_fill_manual(values = microbe_colors) +
            labs(y = 'Relative Abundance',
                 tag = 'A') +
            theme_classic() +
            theme(legend.position = 'none',
                  panel.background = element_rect(colour = 'black'),
                  axis.ticks.y = element_line(colour = 'black'),
                  axis.minor.ticks.y.left = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title = element_text(colour = 'black', size = 14),
                  axis.text = element_text(colour = 'black', size = 10),
                  plot.tag = element_text(colour = 'black', size = 18, 
                                            face = 'bold', hjust = 0))
        ),
        
        set_sizes=(
          upset_set_size(geom = geom_bar(fill = c('gray65', '#F21A00', '#3B9AB2')),
                         position = 'right') + #'gray65'
            geom_text(aes(label = after_stat(count)), hjust = 1.1, stat = 'count', 
                      colour = 'white', fontface = 'bold', size = 4) +
          # + annotate(geom='text', label='@', x='Drama', y=850, color='white', size=3) +
            # expand_limits(y = 200) +
            labs(tag = 'D') +
            theme_classic() +
            theme(panel.background = element_rect(colour = 'black'),
                  axis.ticks.x = element_line(colour = 'black'),
                  axis.title = element_text(colour = 'black', size = 14),
                  axis.text.x = element_text(angle = 0,  colour = 'black', size = 10),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.tag = element_text(colour = 'black', size = 18, 
                                          face = 'bold', hjust = 0))
        ),
        
        matrix = (
          intersection_matrix() +
            labs(tag = 'C') +
            theme_classic() 
        ),
        # min_size = 10,
        sort_intersections=FALSE,
        intersections=list(c('Healthy', 'Time'),
                           'Healthy',
                           c('Disease', 'Time'),
                           'Disease',
                           'Time',
                           'Outside of known sets'),
        sort_sets = FALSE, 
        name = NULL,
        themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(color = c('black', '#F21A00', '#3B9AB2'),
                                                                                                    size = 14, face = 'bold'),
                                                                         panel.background = element_rect(colour = 'black'),
                                                                         panel.border = element_rect(colour = 'black', fill = 'transparent'),
                                                                         plot.tag = element_text(colour = 'black', size = 18, 
                                                                                                 face = 'bold', hjust = 0),
                                                                         axis.ticks = element_line(colour = 'black'),
                                                                         axis.line = element_line(colour = "black", 
                                                                                                  linewidth = rel(1))))),
        upset_stripes(colors = NA)
  )

base_plot

# plot_grid(base_plot, 
#           cowplot::plot_grid(NULL, colour_options$legend, NULL, ncol = 1, rel_heights = c(0.1, 1, 0.4)), 
#           rel_widths = c(1, .25))
# ggsave('../Results/Fig4_traditional_complexUpset.png', 
#        height = 12, width = 10, bg = 'white')


ggdraw(base_plot) +
  draw_plot(colour_options$legend, 
            x = 0.75, y = 0.25, width = 0.25, height = 0.75)
ggsave('../Results/Fig4_traditional_complexUpset.png', 
       height = 12, width = 10, bg = 'white')
ggsave('../Results/Fig4.svg', 
       height = 12, width = 10, bg = 'white')



select(fixed_models, -where(is.list)) %>%
  mutate(across(health:healthXtimepoint, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepoint, ~. < alpha)) %>%
  filter(if_any(contains('health'))) %>%
  count(up_disease)



select(fixed_models, -where(is.list)) %>%
  mutate(across(health:healthXtimepoint, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepoint, ~. < alpha)) %>%
  filter(if_any(contains('health'))) %>%
  count(up_disease, genus) %>%
  count(up_disease)


#### Composition Analysis ####
tst_data <- base_plot_data %>%
  
  filter(!if_all(where(is.logical), ~!.)) %>%
  mutate(Time = Time | `Disease State x Time`, .keep = 'unused') %>%
  filter(`Disease Associated` | `Healthy Associated`) %>%
  select(-Time) %>%

  group_by(across(c(genus, where(is.logical)))) %>%
  summarise(n = n_distinct(asv_id),
            .groups = 'drop') %>%
  rowwise %>%
  
  mutate(allFalse = all(!c_across(where(is.logical)))) %>%
  ungroup %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = where(is.logical),
               names_to = 'term',
               values_to = 'sig') %>%
  filter(sig) %>%
  group_by(row_id, genus, n) %>%
  summarise(intersection = str_c(term, collapse = '; '),
            .groups = 'drop') %>%
  select(-row_id) 

# tst_data <- tst_data %>%
#   mutate(new_intersection = case_when(intersection == 'allFalse' ~ 'allFalse',
#                                       intersection == 'Sampling Time' ~ 'Sampling Time',
#                                       str_detect(intersection, 'Health State \\(Healthy\\)') ~ 'healthy',
#                                       str_detect(intersection, 'Health State \\(Diseased\\)') ~ 'disease',
#                                       TRUE ~ intersection)) %>%
#   group_by(genus, new_intersection) %>%
#   summarise(n = sum(n), 
#             .groups = 'drop')

tst_data %>%
  filter(genus != 'NA') %>%
  # filter(sum(n) > 5, .by = 'genus') %>%
  # filter(sum(n) > 7, .by = 'intersection') %>%
  pivot_wider(names_from = 'intersection',
              values_from = 'n',
              values_fill = 0L) %>%
  column_to_rownames('genus') %>%
  chisq.test(simulate.p.value = FALSE)


# x <- 1; k <- 3; m <- 24; N <- 342
row_fisher <- function(x, k, m, N, direction = 'two.sided'){
  ##https://dputhier.github.io/ASG/practicals/go_statistics_td/go_statistics_td_2015.html
  data.frame(significant = c(x, m - x),
             not_significant = c(k - x, N - m - (k - x))) %>%
    fisher.test(alternative = direction) %>%
    tidy
}


tst_data %>%
  filter(genus != 'NA') %>%
  rename(n_asv = n) %>%
  # mutate(n_genus = sum(n_asv), .by = 'genus') %>%
  left_join(rename(n_genera, n_genus = n),
            by = 'genus') %>%
  mutate(n_intersect = sum(n_asv), .by = 'intersection') %>%
  mutate(total = sum(n_genera$n)) %>%
  
  # filter(sum(n_asv) > 5, .by = 'genus') %>%
  filter(n_genus > 5) %>%
  filter(sum(n_asv) > 7, .by = 'intersection') %>%
  rowwise() %>%
  mutate(row_fisher(n_asv, n_genus, n_intersect, total, direction = 'greater')) %>%
  ungroup %>%
  # group_by(new_intersection) %>%
  # mutate(p_adj = p.adjust(p.value, 'fdr')) %>%
  filter(p.value < 0.05) %>%
  select(genus, intersection, n_asv, n_genus, estimate, p.value)


######### EXPERIMENTAL ########
#### Overrepresentation of Orders ####
tst_data <- select(fixed_models, -where(is.list)) %>%
  mutate(across(health:healthXtimepoint, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepoint, ~. < alpha)) %>%
  filter(if_any(contains('health'))) %>%
  count(order, family, genus, up_disease) %>%
  mutate(up_disease = if_else(up_disease, 'disease', 'healthy')) %>%
  pivot_wider(names_from = up_disease,
              values_from = n,
              values_fill = 0) %>%
  full_join(count(fixed_models, order, family, genus) %>%
              rename(total = n),
            by = c('order', 'family', 'genus')) %>%
  mutate(across(c(healthy, disease), ~replace_na(., 0L))) %>%
  arrange(-total) 

tst_data %>%
  filter(total > 2) 

bind_rows(filter(tst_data, total > 3),
          summarise(tst_data, across(where(is.integer), sum), .by = c('order', 'family')) %>%
            filter(total > 5),
          summarise(tst_data, across(where(is.integer), sum), .by = c('order')) %>%
            filter(total > 5)) 


tst_data %>%
  filter(total > 2) %>%
  mutate(neither = total - healthy - disease) %>%
  select(-total) %>%
  filter(healthy > 0 | disease > 0) %>%
  rowwise %>%
  # mutate(binom.test(x = (healthy + disease), 
  #                   n = (healthy + disease + neither),
  #                   alternative = 'greater') %>% 
  #          broom::tidy() %>%
  #          select(p.value) %>%
  #          rename(p_associated = p.value)) %>%
  # filter(p_associated < 0.05) %>%
  mutate(binom.test(x = healthy, 
                    n = (healthy + disease + neither),
                    alternative = 'greater') %>% 
           broom::tidy() %>%
           select(p.value) %>%
           rename(p_healthy = p.value),
         
         binom.test(x = disease,
                    n = (healthy + disease + neither),
                    alternative = 'greater') %>% 
           broom::tidy() %>%
           select(p.value) %>%
           rename(p_disease = p.value)) %>%
  ungroup() %>%
  filter(p_healthy < 0.05 | 
           p_disease < 0.05)


matrix(
  c(72, 20, 28, 
    17, 0, 3, 
    1, 8, 0,
    8, 0, 1,
    6, 0, 0,
    5, 0, 0), 
  nrow = 6, 
  byrow = TRUE,
  dimnames = list(
    genus = c("Endozoicomonas", "Alteromonas", 'Thalassotalea', 
              'Pelagibacter', 'Halomonas', 'Sansalvadorimonas'),
    health = c("healthy", "disease", 'neither'))) %>%
  chisq.test(simulate.p.value = TRUE)




tst_data %>%
  mutate(neither = total - healthy - disease) %>%
  select(-total) %>%
  filter(healthy > 0 | disease > 0) %>%
  rowwise %>%
  mutate(binom.test(x = healthy, 
                    n = (healthy + disease + neither),
                    alternative = 'greater') %>% 
           broom::tidy() %>%
           select(p.value) %>%
           rename(p_healthy = p.value),
         
         binom.test(x = disease,
                    n = (healthy + disease + neither),
                    alternative = 'greater') %>% 
           broom::tidy() %>%
           select(p.value) %>%
           rename(p_disease = p.value)) %>%
  ungroup() %>%
  filter(p_healthy < 0.05 | 
           p_disease < 0.05)
  

3/10
(29+11+9+12+4+3) 

row_fisher <- function(x, k, m, N, direction = 'two.sided'){
  ##https://dputhier.github.io/ASG/practicals/go_statistics_td/go_statistics_td_2015.html
  data.frame(significant = c(x, m - x),
             not_significant = c(k - x, N - m - (k - x))) %>%
    fisher.test(alternative = direction) %>%
    tidy
}

select(tst_data, -disease) %>%
  column_to_rownames('order') %>%
  chisq.test()

select(tst_data, -disease) %>%
  mutate(total_healthy = sum(healthy),
         total_total = sum(total)) %>%
  filter(total > 1) %>%
  rowwise %>%
  mutate(row_fisher(healthy, total, total_healthy, total_total, direction = 'greater')) %>%
  filter(p.value < 0.05)

select(tst_data, -healthy) %>%
  mutate(total_disease = sum(disease),
         total_total = sum(total)) %>%
  filter(total > 1) %>%
  rowwise %>%
  mutate(row_fisher(disease, total, total_disease, total_total, direction = 'greater')) %>%
  filter(p.value < 0.05)



tst_data %>%
  mutate(neither = total - disease - healthy) %>%
  # filter(total > 5) %>%
  select(-total) %>%
  column_to_rownames('order') %>%
  chisq.test(B = 2000, simulate.p.value = TRUE)

library(chisq.posthoc.test)

tst_data %>%
  mutate(neither = total - disease - healthy) %>%
  # filter(total > 5) %>%
  select(-total, -neither) %>%
  filter(healthy > 0 | disease > 0) %>%
  column_to_rownames('order') %>%
  chisq.test(B = 2000, simulate.p.value = TRUE)

tst_data %>%
  mutate(neither = total - disease - healthy) %>%
  # filter(total > 5) %>%
  select(-total) %>%
  # filter(healthy > 0 | disease > 0) %>%
  column_to_rownames('order') %>%
  chisq.posthoc.test(method = 'fdr',
                     B = 2000, simulate.p.value = TRUE) %>%
  as_tibble() %>% 
  pivot_wider(names_from = Value,
              values_from = c('healthy', 'disease', 'neither')) %>%
  pivot_longer(cols = -Dimension,
               names_to = c('health_association', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  mutate(across(c(Residuals, `p values`), ~parse_number(.))) %>%
  filter(`p values` < 0.05)

#### Try 2 ####
tmp_data <- select(fixed_models, -where(is.list)) %>%
  mutate(across(health:healthXtimepoint, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepoint, ~. < alpha)) %>%
  count(order, health, timepoint, healthXtimepoint) %>%
  # filter(order == 'Alteromonadales') %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = where(is.logical)) %>%
  filter(value) %>%
  group_by(row_id, order, n) %>%
  summarise(combo = str_c(name, collapse = ' + '),
            .groups = 'drop') %>%
  pivot_wider(names_from = combo, values_from = n,
              values_fill = 0L) %>%
  select(-row_id) %>%
  group_by(order) %>%
  summarise(across(everything(), sum)) 

column_to_rownames(tmp_data, 'order') %>%
  chisq.test(B = 2000, simulate.p.value = TRUE)

tmp_data %>%
  pivot_longer(cols = -order) %>%
  rowwise %>%
  mutate(name = str_split(name, ' \\+ ')) %>%
  unnest(name) %>%
  group_by(order, name) %>%
  summarise(value = sum(value),
            .groups = 'drop') %>%
  pivot_wider() %>%
  full_join(count(fixed_models, order) %>%
              rename(total = n),
            by = 'order') %>%
  mutate(across(where(is.integer), ~replace_na(., 0))) %>%
  select(order, health, total) %>%
  column_to_rownames('order') %>%
  chisq.test(B = 2000, simulate.p.value = TRUE)

#### Try 3 ####
upset_data <- select(fixed_models, -where(is.list)) %>%
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
  group_by(order, across(where(is.logical))) %>%
  summarise(n = n(),
            .groups = 'drop') %>%
  # filter(order %in% levels(colour_options$asv_clumping$Top_order)) %>%
  mutate(row_id = row_number()) %>%
  # rowwise %>%
  mutate(allFALSE = !`Sampling Time` & !`Health State x Sampling Time` & 
           !`Health State (Diseased)` & !`Health State (Healthy)`) %>%
  pivot_longer(cols = where(is.logical)) %>%
  filter(value) %>%
  select(-value) %>%
  group_by(order, n, row_id) %>%
  summarise(name = str_c(name, collapse = '+')) %>%
  select(-row_id) %>%
  group_by(order, name) %>%
  summarise(n = sum(n),
            .groups = 'drop') %>%
  pivot_wider(names_from = name,
              values_from = n,
              values_fill = 0L) 

upset_data %>%
  column_to_rownames('order') %>%
  chisq.test(simulate.p.value = TRUE)


upset_data %>%
  column_to_rownames('order') %>%
  chisq.posthoc.test(method = 'none',
                     simulate.p.value = TRUE) %>%
  as_tibble() %>% 
  pivot_wider(names_from = Value,
              values_from = -c(Dimension, Value)) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = -Dimension,
               names_to = c('health_association', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  mutate(across(c(Residuals, `p values`), ~parse_number(.)),
         effect = case_when(`p values` < 0.05 & Residuals > 0 ~ '+',
                            `p values` < 0.05 & Residuals < 0 ~ '-',
                            TRUE ~ '0')) %>%
  # filter(`p values` < 0.05) %>%
  select(-`p values`, -Residuals) %>%
  pivot_wider(names_from = health_association,
              values_from = 'effect',
              values_fill = '0') %>%
  arrange(Dimension)
  


tmp <- upset_data %>%
  pivot_longer(cols = where(is.integer),
               names_to = 'category',
               values_to = 'n') %>%
  nest_by(category) %>%
  mutate(data = list(full_join(data,
                               count(fixed_models, order) %>%
                                 rename(total = n),
                               by = 'order') %>%
                       mutate(across(where(is.integer), ~replace_na(., 0L)))))

tmp %>%
  mutate(data = list(filter(data, n > 0)),
         n_order = nrow(data)) %>%
  filter(n_order > 1) %>%
  mutate(column_to_rownames(data, 'order') %>%
           chisq.test(simulate.p.value = TRUE) %>%
           broom::tidy()) %>%
  filter(p.value < 0.05) %>%
  mutate(column_to_rownames(data, 'order') %>%
           chisq.posthoc.test(method = 'fdr', 
                              simulate.p.value = TRUE) %>%
           as_tibble() %>%
           select(-total) %>%
           pivot_wider(names_from = 'Value',
                       values_from = 'n') %>%
           mutate(across(c(Residuals, `p values`), ~as.character(.)),
                  across(c(Residuals, `p values`), ~parse_number(.)),
                  effect = case_when(`p values` < 0.05 & Residuals > 0 ~ '+',
                                     `p values` < 0.05 & Residuals < 0 ~ '-',
                                     TRUE ~ '0')) %>%
           select(Dimension, effect) %>%
           pivot_wider(names_from = 'Dimension',
                       values_from = 'effect')) %>% View

