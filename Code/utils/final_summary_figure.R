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
                                levels = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'))) %>%
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
  select(asv_id, all_of(model_list), response, SE, starts_with('asymp')) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(asv_id = fct_reorder(asv_id, response, .desc = TRUE)) %>%
  mutate(name = if_else(taxon_level == 'genus', str_c(name, ' sp.'), name),
         name = if_else(taxon_level %in% c('genus', 'species'), str_c('<i>', name, '</i>'), name),
         name = str_c(name, ' (', asv_id, ')'),
         taxon_name = case_when(taxon_level == 'family' ~ str_c(order, name, sep = '; '),
                                 taxon_level == 'order' ~ name,
                                 TRUE ~ str_c(order, family, name, sep = '; ')),
         taxon_name = fct_reorder(taxon_name, response, .desc = TRUE))

asv_shaps <- read_csv('../../intermediate_files/model_shaps.csv.gz', show_col_types = FALSE) %>%
  group_by(wflow_id) %>%
  # mutate(across(ends_with('value'), ~scale(.)[,1])) %>%
  mutate(across(ends_with('shap'), ~scale(.)[,1])) %>%
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
  rowwise(asv_id) %>%
  reframe(update(tank_posthoc, side = '=') %>%
            broom::tidy(conf.int = TRUE) %>%
            select(-p.value),
          broom::tidy(tank_posthoc, conf.int = TRUE) %>%
            select(p.value)) %>%
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
         .by = contrast) %>%
  left_join(select(asv_rankings, asv_id, taxon_name),
            by = 'asv_id')

#### Ranking Plot ####
rank_plot <- asv_rankings %>%
  ggplot(aes(x = response, y = taxon_name)) +
  geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL), width = 0.1) + 
  geom_point() +
  labs(x = 'Feature Rank',
       y = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.y = element_markdown())

#### Shap Plots ####
shap_plot <- asv_shaps %>%
  # filter(wflow_id == 'base_mlp') %>%
  mutate(health = case_when(health == 'D' ~ 'Diseased',
                            health == 'H' ~ 'Healthy')) %>%
  group_by(sample_id, asv_id, health, value, taxon_name) %>%
  summarise(shap = mean(shap),
            .groups = 'drop') %>%
  ggplot(aes(y = taxon_name, x = shap, #, shape = health
             colour = value)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_quasirandom(aes(group = 1), groupOnX = FALSE, dodge.width = 1) +
  scale_color_gradient2(midpoint = log(10, base = 2), 
                        breaks = seq(8, 20, by = 4),
                        trans = scales::log_trans(base = 2)) +
  # guides(shape = guide_legend(override.aes = list(size = 4))) +
  labs(colour = 'Normalized<br>log<sub>2</sub>(CPM)',
       shape = 'Disease<br>State',
       y = NULL,
       x = 'SHAP') +
  theme_classic() +
  theme(legend.title = element_markdown(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.y = element_markdown())

#### Field Plot ####
field_plot <- field_models %>%
  mutate(timepoint = str_replace_all(timepoint, '_', ' '),
         timepoint_sig = if_else(fdr < alpha, timepoint, NA_character_)) %>%
  ggplot(aes(y = taxon_name, x = estimate, shape = timepoint,
             fill = fdr < alpha)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.1, position = position_dodge(0.5),
                show.legend = FALSE) + 
  geom_point(position = position_dodge(0.5)) +
  scale_fill_manual(values = c('TRUE' = 'black', 'FALSE' = 'white')) +
  scale_shape_manual(values = c('2016 Jan' = 'square filled', '2016 Jul' = 'diamond filled', 
                       '2017 Jan' = 'triangle filled', '2017 Jul' = 'triangle down filled')) +
  guides(fill = guide_legend(override.aes = list(shape = 'circle filled', size = 4)),
         shape = guide_legend(override.aes = list(size = 4, fill = 'black'))) +
  labs(x = 'log<sub>2</sub>(D/H)',
       y = NULL,
       shape = 'Time',
       fill = 'Significance') +
    theme_classic() +
    theme(axis.title.x = element_markdown(),
          legend.title = element_markdown(colour = 'black', size = 14),
          legend.text = element_text(colour = 'black', size = 10),
          panel.background = element_rect(colour = 'black'),
          axis.text = element_text(colour = 'black', size = 10),
          axis.title = element_text(colour = 'black', size = 14),
          axis.text.y = element_markdown()) 

#### Tank Plot ####
tank_plot <- tank_models %>%
  filter(!contrast %in% c('Outcome', 'PostvPreH', 'DHvNH', 'Post (Diseased / Control)')) %>%
  mutate(contrast = case_when(contrast == 'Post (Diseased / Healthy)' ~ 'Outcome',
                              contrast == 'Diseased (Post / Pre)' ~ 'Time',
                              TRUE ~ contrast),
         contrast = factor(contrast, levels = c('Time', 'Exposure', 'Outcome')),
         contrast_sig = if_else(fdr < alpha, contrast, NA_character_),
         contrast_sig = factor(contrast_sig, levels = levels(contrast))) %>%
  ggplot(aes(y = taxon_name, x = estimate, colour = contrast,
             fill = contrast_sig)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.1, position = position_dodge(0.5),
                show.legend = FALSE) + 
  geom_point(position = position_dodge(0.5), shape = 'circle filled') +
  guides(fill = 'none',
         colour = guide_legend(override.aes = list(shape = 'circle', size = 4))) +
  scale_fill_discrete(na.value = 'white') +
  labs(x = 'log<sub>2</sub>(Contrast)',
       y = NULL,
       colour = 'Contrast',
       alpha = 'Significance') +
  theme_classic() +
  theme(axis.title.x = element_markdown(),
        legend.title = element_markdown(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.y = element_markdown())


#### Combined Plots ####
(rank_plot + shap_plot + field_plot + tank_plot) +
  plot_layout(nrow = 1, guides = 'collect', axes = 'collect_y') +
  plot_annotation(tag_levels = 'A')
ggsave('../../Results/overview_results.png', height = 7, width = 12)
