library(tidyverse)
library(emmeans)
library(ggbeeswarm)
library(ggtext)
library(patchwork)

alpha <- 0.05

#### Functions ####
top_classification <- function(data, threshold = 80){
  rename_with(data, .cols = domain:species,
              ~str_c(., '_name')) %>%
    pivot_longer(cols = -asv_id,
                 names_to = c('taxon_level', '.value'),
                 names_pattern = '(.*)_(.*)') %>%
    mutate(taxon_level = factor(taxon_level, ordered = TRUE,
                                levels = c('domain', 'phylum', 'class', 'order', 
                                           'family', 'genus', 'species'))) %>%
    filter(confidence >= threshold) %>%
    group_by(asv_id) %>%
    filter(taxon_level == max(taxon_level)) %>%
    ungroup %>%
    left_join(select(data, -ends_with('confidence')), by = 'asv_id')
  
}

#### Data ####
taxonomy <- read_csv('../../intermediate_files/update_taxonomy.csv',
                     show_col_types = FALSE) %>%
  top_classification() %>%
  left_join(read_csv('../../intermediate_files/update_taxonomy.csv',
                     show_col_types = FALSE) %>%
              select(asv_id, ends_with('confidence')),
            by = 'asv_id') %>%
  bind_rows(read_csv('../../intermediate_files/taxonomy.csv.gz', show_col_types = FALSE) %>%
              filter(asv_id == 'ASV131') %>%
              mutate(name = genus,
                     taxon_level = 'genus'))
  

  
model_list <- read_csv('../../Results/equivilant_top_models.csv.gz',
                       show_col_types = FALSE) %>%
  filter(pract_equiv >= 0.8) %>%
  rename(wflow_id = model) %>%
  pull(wflow_id)

asv_rankings <- read_csv('../../Results/asv_importance.csv.gz',
                            show_col_types = FALSE) %>%
  filter(p_adjust < alpha) %>%
  filter(!is.na(response)) %>%
  select(asv_id, all_of(model_list), response, SE, starts_with('asymp'), p_adjust) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(asv_id = fct_reorder(asv_id, response, .desc = TRUE)) %>%
  mutate(name = if_else(taxon_level == 'genus', str_c(name, ' sp.'), name),
         name = if_else(taxon_level %in% c('genus', 'species'), str_c('<i>', name, '</i>'), name),
         # name = str_c(name, ' (', asv_id, ')'),
         name = str_c(asv_id, ' - ', name),
         taxon_name = case_when(taxon_level == 'family' ~ str_c(name, sep = '; '),
                                 taxon_level == 'order' ~ name,
                                 TRUE ~ str_c(name, family, sep = '; ')),
         taxon_name = fct_reorder(taxon_name, response, .desc = TRUE)) %>%
  filter(!is.na(taxon_name))

asv_shaps <- read_csv('../../intermediate_files/model_shaps.csv.gz', show_col_types = FALSE) %>%
  group_by(wflow_id) %>%
  # mutate(across(ends_with('value'), ~scale(.)[,1])) %>%
  # mutate(across(ends_with('shap'), ~scale(.)[,1])) %>%
  ungroup %>%
  filter(wflow_id %in% model_list) %>%
  pivot_longer(starts_with('ASV'),
               names_to = c('asv_id', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  filter(asv_id %in% as.character(asv_rankings$asv_id)) %>%
  mutate(asv_id = factor(asv_id, levels = levels(asv_rankings$asv_id))) %>%
  left_join(select(asv_rankings, asv_id, taxon_name),
            by = 'asv_id')

field_models <- read_rds('../../intermediate_files/field_asv_models.rds.gz') %>%
  select(asv_id, all_health_fit) %>%
  filter(asv_id %in% as.character(asv_rankings$asv_id)) %>%
  rowwise(asv_id) %>%
  reframe(contrast(all_health_fit, 'pairwise') %>% broom::tidy(conf.int = TRUE)) %>%
  
  mutate(asv_id = factor(asv_id, levels = levels(asv_rankings$asv_id))) %>%
  mutate(timepoint = str_c(year, season, sep = '_') %>% 
           str_replace_all(c('S' = 'Jul', 'W' = 'Jan')),
         fdr = p.adjust(p.value, 'fdr'),
         .by = c(year, season)) %>%
  left_join(select(asv_rankings, asv_id, taxon_name),
            by = 'asv_id')


tank_models <- read_rds('../../intermediate_files/tank_asv_models.rds.gz') %>% 
  select(asv_id, tank_posthoc) %>%
  filter(asv_id %in% as.character(asv_rankings$asv_id)) %>%
  unnest(tank_posthoc) %>%
  mutate(asv_id = factor(asv_id, levels = levels(asv_rankings$asv_id))) %>%
  # filter(contrast %in% c('DvH', 'DvN', 'DDvDH', 'DDvNH', 'DHvNH')) %>%
  mutate(contrast = case_when(contrast == 'DvH' ~ 'Outcome',
                              contrast == 'DvN' ~ 'Exposure',
                              contrast == 'DDvDH' ~ 'Post (Diseased / Healthy)',
                              contrast == 'DDvNH' ~ 'Post (Diseased / Control)',
                              contrast == 'PostvPreD' ~ 'Diseased (Post / Pre)',
                              TRUE ~ contrast)) %>%
  left_join(summarise(field_models, all_sig = all(fdr < alpha), .by = asv_id),
            by = 'asv_id') %>%
  mutate(across(c(estimate, p.value, conf.low, conf.high),
                ~if_else(!all_sig, NA_real_, .))) %>%
  mutate(fdr = p.adjust(p.value, 'fdr'),
         .by = c(contrast, all_sig)) %>% #
  left_join(select(asv_rankings, asv_id, taxon_name),
            by = 'asv_id')

#
#### Ranking Plot ####
count(asv_rankings, genus)

rank_plot <- asv_rankings %>%
  ggplot(aes(x = response, y = taxon_name, fill = p_adjust < alpha)) + #
  geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL), width = 0.1) + 
  geom_point(data = . %>%
               select(asv_id, taxon_name, starts_with('base')) %>%
               pivot_longer(cols = starts_with('base'),
                            names_to = 'model',
                            values_to = 'rank',
                            names_prefix = 'base_'),
             aes(x = rank, colour = model), fill = 'black',
             position = position_jitter(height = 0.2, width = 2.5),
             size = 1) +
  geom_point(shape = 'circle filled') +
  scale_fill_manual(values = c('TRUE' = 'black', 'FALSE' = 'white')) +
  scale_colour_discrete(labels = c('forest' = 'RF', 
                                 'lasso' = 'LASSO',
                                 'mlp' = 'MLP',
                                 'svmLinear' = 'SVM'),
                        breaks = c('mlp', 'svmLinear',
                                   'lasso', 'forest')) +
  scale_x_continuous(breaks = c(0, 25, 50, 75)) +
  guides(fill = 'none',
         colour = guide_legend(override.aes = list(size = 4))) +
  labs(x = 'ASV Rank',
       y = NULL,
       colour = 'ML Model') +
  theme_classic() +
  theme(legend.title = element_markdown(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        legend.key = element_blank(),
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.y = element_markdown())

#### Shap Plots ####
shap_associations <- asv_shaps %>%
  group_by(sample_id, asv_id, health, value, taxon_name) %>%
  summarise(shap = mean(shap),
            .groups = 'drop') %>%
  group_by(asv_id) %>%
  reframe(broom::tidy(lm(shap ~ value))) %>%
  filter(term != '(Intercept)') %>%
  arrange(desc(asv_id)) %>%
  mutate(disease_associated = estimate > 0) %>%
  select(asv_id, disease_associated)

count(shap_associations, disease_associated)

shap_associations %>%
  mutate(association = if_else(disease_associated, 'Disease', 'Healthy'), 
         .keep = 'unused') %>%
  left_join(taxonomy, by = 'asv_id') %>%
  select(-taxon_level:-confidence) %>%
  select(-contains('domain'), -contains('phylum'), -contains('class')) %>%
  rename_with(.cols = c(order:species), ~str_c(., '_name')) %>%
  pivot_longer(cols = -c(asv_id, association),
               names_to = c('taxon_level', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  mutate(name = str_c(name, ' (', str_replace_na(scales::percent(confidence, scale = 1)), ')'),
         .keep = 'unused') %>%
  pivot_wider(names_from = taxon_level,
              values_from = name) %>%
  rename(shap_assocation = association) %>%
  write_csv('../../Results/shap_assications.csv')
    

shap_plot <- asv_shaps %>%
  mutate(value = value - min(value), .by = c(asv_id, wflow_id)) %>%
  # filter(wflow_id == 'base_mlp') %>%
  mutate(health = case_when(health == 'D' ~ 'Diseased',
                            health == 'H' ~ 'Healthy')) %>%
  group_by(sample_id, asv_id, health, value, taxon_name) %>%
  summarise(shap = mean(shap),
            .groups = 'drop') %>%
  
  mutate(value = if_else(health == 'Diseased', value, -1 * value),
         value = case_when(value > 10 ~ 10,
                           value < -10 ~ -10,
                           TRUE ~ value)) %>%
  
  ggplot(aes(y = taxon_name, x = shap, #, shape = health
             fill = value)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_quasirandom(aes(group = 1), groupOnX = FALSE, dodge.width = 1, 
                   shape = 'circle filled') +
  scale_fill_gradient2(midpoint = 0, 
                       low = wesanderson::wes_palette("Zissou1", 2, type = "continuous")[1],
                       high = wesanderson::wes_palette("Zissou1", 2, type = "continuous")[2],
                       labels = c('Healthy', '5', '0', '5', 'Diseased')) +
  # guides(shape = guide_legend(override.aes = list(size = 4))) +
  labs(fill = 'Normalized<br>log<sub>2</sub>(CPM)',
       shape = 'Disease<br>State',
       y = NULL,
       x = 'SHAP') +
  theme_classic() +
  theme(legend.title = element_markdown(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        legend.key = element_blank(),
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.y = element_markdown())

#
#### Field Plot ####
field_plot <- field_models %>%
  mutate(timepoint = str_replace_all(timepoint, '_', ' '),
         timepoint_sig = if_else(fdr < alpha, timepoint, NA_character_),
         direction_sig = case_when(fdr > alpha ~ 'Neither',
                                   estimate > 0 ~ 'Diseased',
                                   estimate < 0 ~ 'Healthy'), 
         direction_sig = factor(direction_sig, levels = c('Healthy', 'Neither', 'Diseased')),
         timepoint = factor(timepoint, levels = c('2017 Jul', '2017 Jan', '2016 Jul', '2016 Jan'))) %>%
  ggplot(aes(y = taxon_name, x = estimate, #shape = timepoint,
             fill = direction_sig, group = timepoint)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.1, position = position_dodge(0.5),
                show.legend = FALSE) + 
  geom_point(position = position_dodge(0.5), shape = 'circle filled') +
  scale_fill_manual(values = set_names(c('white', wesanderson::wes_palette("Zissou1", 2, type = "continuous")),
                                       c('Neither', 'Healthy', 'Diseased')),
                    breaks = c('Diseased', 'Healthy')) +
  scale_x_continuous(limits = c(-7, 12), breaks = c(-5, 0, 5)) +
  # scale_fill_manual(values = c('TRUE' = 'black', 'FALSE' = 'white')) +
  # scale_shape_manual(values = c('2016 Jan' = 'square filled', '2016 Jul' = 'diamond filled', 
  #                      '2017 Jan' = 'triangle filled', '2017 Jul' = 'triangle down filled')) +
  guides(fill = guide_legend(override.aes = list(shape = 'circle filled', size = 4)),
         shape = guide_legend(override.aes = list(size = 4, fill = 'black'))) +
  labs(x = 'Field log<sub>2</sub>(D/H)',
       y = NULL,
       shape = 'Time',
       fill = 'Disease\nAssociation') +
    theme_classic() +
    theme(axis.title.x = element_markdown(),
          legend.title = element_text(colour = 'black', size = 14),
          legend.text = element_text(colour = 'black', size = 10),
          legend.key = element_blank(),
          panel.background = element_rect(colour = 'black'),
          axis.text = element_text(colour = 'black', size = 10),
          axis.title = element_text(colour = 'black', size = 14),
          axis.text.y = element_markdown()) 



#### Tank Plot ####
#Maybe have untested as grey - maybe make them blank
tank_plot <- tank_models %>%
  filter(!contrast %in% c('Outcome', 'PostvPreH', 'DHvNH', 
                          'Post (Diseased / Control)',
                          'Diseased (Post / Pre)')) %>%
  mutate(contrast = case_when(contrast == 'Post (Diseased / Healthy)' ~ 'Outcome',
                              contrast == 'Diseased (Post / Pre)' ~ 'Time',
                              TRUE ~ contrast),
         contrast = factor(contrast, levels = rev(c('Exposure', 'Outcome', 'Field'))),
         
         direction_sig = case_when(!all_sig ~ 'Untested',
                                   fdr > alpha ~ 'Neither',
                                   estimate > 0 ~ 'Diseased',
                                   estimate < 0 ~ 'Healthy'), 
         direction_sig = factor(direction_sig, levels = c('Healthy', 'Neither', 'Diseased', 
                                                          'Untested'))) %>%
  ggplot(aes(y = taxon_name, x = estimate, shape = contrast, #
             fill = direction_sig)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.1, position = position_dodge(0.5),
                show.legend = FALSE) + 
  geom_point(position = position_dodge(0.5), show.legend = TRUE) +
  geom_label(aes(x = 0, label = if_else(is.na(estimate), 'nt', NA_character_)),
             fill = alpha('white', 0.5), label.size = NA) +
  
  guides(fill = guide_legend(override.aes = list(shape = 'circle filled', size = 4)),
         shape = guide_legend(override.aes = list(size = 4, fill = 'black'))) +
  scale_x_continuous(limits = c(-7, 12), breaks = c(-5, 0, 5)) +
  scale_fill_manual(values = set_names(c('white', #'grey50',
                                         wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous")),
                                       c('Neither', 'Healthy', 'Diseased')),
                    breaks = c('Diseased', 'Healthy'), drop = FALSE) +
  scale_shape_manual(values = c('Field' = 'circle filled', 
                                'Exposure' = 'square filled', 
                                'Outcome' = 'diamond filled'), 
                     breaks = c('Exposure', 'Outcome', 'Field'),
                     drop = FALSE) +
  labs(x = 'Tank log<sub>2</sub>(D/H)',
       y = NULL,
       fill = 'Disease\nAssociation',
       shape = 'Effect') +
  theme_classic() +
  theme(axis.title.x = element_markdown(),
        legend.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        legend.key = element_blank(),
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.y = element_markdown())


#### Combined Plots ####
(rank_plot + shap_plot + field_plot + tank_plot) +
  plot_layout(nrow = 1, guides = 'collect', axes = 'collect_y') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag.position = 'topleft',
        plot.tag.location = 'panel',
        plot.tag = element_text(vjust = 5, size = 16),
        plot.margin = margin(t = 10))
ggsave('../../Results/Fig5_overview_results.png', height = 7, width = 12)



