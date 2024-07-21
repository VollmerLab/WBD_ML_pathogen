#### Libraries ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(lubridate)
library(emmeans)
library(broom)
library(patchwork)
library(multidplyr)
library(ComplexUpset)
library(cowplot)

#### N ASVs ####
alpha <- 0.05

#### Functions ####
make_fc_table <- function(form, model, data){
  x_vars <- str_split(as.character(form)[2], ' [\\+\\*\\|] ') %>% unlist %>% unique %>%
    str_remove_all('\\(|\\)')
  
  out <- emmeans(model, form, data = data) %>%
    contrast(if_else(any(x_vars == 'exp_dis'), 'revpairwise', 'pairwise')) %>%
    as_tibble %>%
    select(any_of(x_vars), estimate)
  
  if(any(x_vars == 'exp_dis')){
    out <- out %>%
      filter(!(exp_dis == 1 & exp_hea == 1)) %>%
      mutate(treatment = case_when(exp_dis == 0 & exp_hea == 0 ~ 'pre',
                                   exp_dis == 1 ~ 'postDisease',
                                   exp_hea == 1 ~ 'postHealth')) %>%
      select(treatment, estimate)
    x_vars <- 'treatment'
  } 
  
  out %>%
    pivot_wider(names_from = any_of(x_vars),
                values_from = estimate,
                names_prefix = 'fcDH_')
}

test_within_timepoint <- function(model, data){
  tp_contrasts <- emmeans(model, ~health | timepoint, data = data) %>%
    contrast('pairwise') %>%
    as_tibble
  
  select(tp_contrasts, timepoint, p.value) %>%
    pivot_wider(names_from = timepoint,
                values_from = p.value,
                names_prefix = 'p.withinYear_') %>%
    mutate(tp_contrasts = list(tp_contrasts))
  
}

cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('lmerTest', 'emmeans', 'tibble', 'nlme',
                           'dplyr', 'tidyr', 'stringr'))
cluster_copy(cluster, c('make_fc_table', 'test_within_timepoint'))

#### Data ####
field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  mutate(timepoint = str_c(season, year, sep = '_'),
         across(domain:genus, str_replace_na)) %>%
  nest_by(across(c(domain:genus, asv_id)))


taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz',
                     show_col_types = FALSE) %>%
  mutate(across(where(is.character), 
                ~str_replace_na(., replacement = '')))

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
  
  mutate(make_fc_table(~health | timepoint, model, data),
         test_within_timepoint(model, data)) %>%
  
  collect %>%
  ungroup %>%
  rename_with(~str_c('p_', .), .cols = c('health', 'timepoint', 'healthXtimepoint')) %>%
  mutate(across(c(starts_with('p_'), starts_with('p.within')), 
                ~p.adjust(., 'fdr'), 
                .names = 'fdr_{.col}')) %>%
  rename_with(~str_remove(., '_p'))


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
  select(-Top_order, -Top_genus) %>%
  rename(the_colour = group) %>%
  mutate(across(fdr_health:fdr_healthXtimepoint, ~. < alpha),
         health_d = fdr_health & up_disease,
         health_h = fdr_health & !up_disease) %>%
  
  rowwise %>%
  mutate(consistent = all(c_across(starts_with('fdr.within')) < alpha) & 
           (all(c_across(starts_with('fcDH')) > 0) | 
              all(c_across(starts_with('fcDH')) < 0)),
         .keep = 'unused') %>%
  ungroup %>%
  
  
  rename('Disease State' = fdr_health,
         'Disease Associated' = health_d,
         'Healthy Associated' = health_h,
         'Time' = fdr_timepoint,
         'Disease State x Time' = fdr_healthXtimepoint) %>%
  select(-up_disease, -`Disease State`) %>%
  mutate(type = if_else(!consistent, 'no', 'yes')) %>%
  select(-consistent) %>%
  identity()

base_plot_data %>%
  filter(`Healthy Associated`) %>%
  count(genus) %>%
  filter(genus != 'NA')

base_plot_data %>%
  mutate(genus = if_else(genus == 'NA', NA_character_, genus)) %>%
  summarise(n_disease = sum(`Disease Associated`),
            n_healthy = sum(`Healthy Associated`),
            
            genus_disease = n_distinct(genus[`Disease Associated`], na.rm = TRUE),
            genus_healthy = n_distinct(genus[`Healthy Associated`], na.rm = TRUE))



base_plot_data %>%
  filter(type == 'yes') %>%
  mutate(genus = if_else(genus == 'NA', NA_character_, genus)) %>%
  summarise(n_disease = sum(`Disease Associated`),
            n_healthy = sum(`Healthy Associated`),
            
            genus_disease = n_distinct(genus[`Disease Associated`], na.rm = TRUE),
            genus_healthy = n_distinct(genus[`Healthy Associated`], na.rm = TRUE))


set_size_replacement <- base_plot_data %>%
  mutate(Time = Time | `Disease State x Time`, .keep = 'unused') %>%
  select(Time, `Disease Associated`, `Healthy Associated`, type) %>%
  group_by(type) %>%
  summarise(across(where(is.logical), sum)) %>%
  pivot_longer(cols = -type,
               names_to = 'group',
               values_to = 'n') %>%
  mutate(print_value = if_else(type == 'no', sum(n), n),
         .by = group) %>%
  mutate(group = factor(group, levels = rev(c('Healthy Associated', 'Disease Associated', 'Time'))),
         type = str_c(group, type, sep= '_')) %>%
  ggplot(aes(y = group, x = n, fill = type)) +
  geom_col() +
  geom_text(aes(x = if_else(print_value == 2, 2 + 10, print_value), label = print_value), hjust = 1.1, #stat = 'count', 
            colour = 'white', fontface = 'bold', size = 4) +
  scale_fill_manual(values = c('Healthy Associated_no' = '#3B9AB2', 'Healthy Associated_yes' = colorspace::darken('#3B9AB2', amount = 0.5), 
                               'Disease Associated_no' = '#F21A00', 'Disease Associated_yes' = colorspace::darken('#F21A00', amount = 0.5), 
                               'Time_no' = 'gray65', 'Time_yes' = 'gray25')) +
  labs(tag = 'D',
       x = 'Set Size') +
  theme_classic() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black'),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.x = element_text(angle = 0,  colour = 'black', size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.tag = element_text(colour = 'black', size = 18, 
                                face = 'bold', hjust = 0))


base_plot <- base_plot_data %>%
  mutate(Time = Time | `Disease State x Time`, .keep = 'unused') %>%
  rename(Diseased = `Disease Associated`,
         Healthy = `Healthy Associated`) %>%
  
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        
        base_annotations = list(
          'Intersection size' = intersection_size(
            counts = TRUE,
            mapping = aes(fill = type)
          ) + 
            geom_text(data = . %>% summarise(n = sum(type == 'yes'), .by = intersection) %>%
                        filter(n > 0),
                      aes(x = intersection, y = 5 + n, label = n), inherit.aes = FALSE) +
            scale_fill_manual(values = c('no' = 'gray65', 'yes' = 'gray25')) +
            scale_y_continuous(labels=scales::comma_format(), 
                               expand = expansion(mult = c(0.01, 0.05))) +
            labs(tag = 'B') +
            theme_classic() +
            theme(legend.position = 'none',
                  panel.background = element_rect(colour = 'black'),
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
                           c('Diseased', 'Time'),
                           'Diseased',
                           'Time',
                           'Outside of known sets'
                           ),
        
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
base_plot[[6]]
set_size_replacement
base_plot[[6]] <- set_size_replacement


ggdraw(base_plot) +
  draw_plot(colour_options$legend, 
            x = 0.75, y = 0.25, width = 0.25, height = 0.75)
ggsave('../Results/Fig3_traditional_complexUpset.png', 
       height = 12, width = 10, bg = 'white')
ggsave('../Results/Fig3.svg', 
       height = 12, width = 10, bg = 'white')

#### Composition Analysis ####
genus_associations <- base_plot_data %>%
  filter(type == 'yes',
         `Disease Associated`) %>%
  
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

genus_associations %>%
  filter(genus != 'NA') %>%
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

genus_associations %>%
  filter(genus != 'NA') %>%
  rename(n_asv = n) %>%
  left_join(rename(n_genera, n_genus = n),
            by = 'genus') %>%
  mutate(n_intersect = sum(n_asv), .by = 'intersection') %>%
  mutate(total = sum(n_genera$n)) %>%
  rowwise() %>%
  mutate(row_fisher(n_asv, n_genus, n_intersect, total, direction = 'greater')) %>%
  ungroup %>%
  filter(p.value < 0.05) %>%
  select(genus, intersection, n_asv, n_genus, estimate, p.value)

#### Output list of ASVs ####
base_plot_data %>%
  select(asv_id, `Disease Associated`, 
         `Healthy Associated`, type) %>%
  rename(consistent = type) %>%
  mutate(consistent = consistent == 'yes') %>%
  write_csv('../intermediate_files/DA_asv_groups.csv')


