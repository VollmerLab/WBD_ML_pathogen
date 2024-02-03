library(tidyverse)
library(magrittr)
library(patchwork)
library(vegan)
library(phyloseq)

#### Data ####
field_data_cpm <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  select(asv_id, sample_id, log2_cpm_norm) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm)


if(!file.exists('../intermediate_files/taxonomy.csv.gz')){
  taxonomy <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
           show_col_types = FALSE) %>%
    select(asv_id, domain:species) %>%
    distinct %>%
    mutate(across(everything(), str_replace_na)) %>%
    write_csv('../intermediate_files/taxonomy.csv.gz')
} else {
  taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz', 
                       show_col_types = FALSE) %>%
    mutate(across(everything(), str_replace_na))
}

field_metadata <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                           show_col_types = FALSE) %>%
  select(sample_id, year, season, site, health) %>%
  distinct %>%
  mutate(timepoint = str_c(year, season, sep = '_'),
         across(where(is.character), as.factor)) %>%
  select(sample_id, site, health, timepoint)

#### Alpha Diversity Metrics ####
microbiome_data <- read_rds("../intermediate_files/prepped_microbiome.rds.gz")


alpha_metric_analysis <- microbiome_data %>%
  subset_taxa(rownames(tax_table(microbiome_data)) %in% colnames(field_data_cpm)) %>%
  microbiome::alpha() %>%
  as_tibble(rownames = 'sample_id') %>%
  right_join(field_metadata,
             by = 'sample_id') %>%
  select(sample_id, site, health, timepoint, observed, diversity_shannon, diversity_inverse_simpson) %>%
  rename(richness = observed) %>%
  pivot_longer(cols = c(richness, diversity_shannon, diversity_inverse_simpson),
               names_to = 'metric', 
               values_to = 'value') %>%
  # mutate(value = case_when(metric == 'diversity_shannon' ~ log(value),
  #                          metric == 'diversity_inverse_simpson' ~ log(value),
  #                          TRUE ~ value)) %>%
  mutate(log_value = log(value),
         sqrt_value = sqrt(value),
         inv_value = 1/value) %>%
  nest_by(metric)

library(emmeans)
source('~/R/R Functions/diagPlots.R')
## Richness ##
richness_aov <- aov(value ~ health * timepoint * site, data = alpha_metric_analysis$data[[3]])
diag.plots(richness_aov, col.nos = c(2:4), data = alpha_metric_analysis$data[[3]])

#Fix site/timepoint variance
richness_gls <- gls(value ~ health * site * timepoint, data = alpha_metric_analysis$data[[3]],
                    weights = varComb(varIdent(form = ~1|health), 
                                      varIdent(form = ~1|timepoint), 
                                      varIdent(form = ~1|site)))

diag.plots(richness_gls, col.nos = c(2:4), data = alpha_metric_analysis$data[[3]])

lmerTest::lmer(value ~ health + (1 | site) + (1 | timepoint),
     data = alpha_metric_analysis$data[[3]]) %>% anova(ddf = 'Kenward-Roger')

## Shannon ##
shannon_aov <- aov(value ~ health * site * timepoint, data = alpha_metric_analysis$data[[2]])
diag.plots(shannon_aov, col.nos = c(2:4), data = alpha_metric_analysis$data[[2]])

shannon_gls <- gls(value ~ health * site * timepoint, data = alpha_metric_analysis$data[[2]],
                    weights = varComb(varIdent(form = ~1|health), 
                                      varIdent(form = ~1|timepoint), 
                                      varIdent(form = ~1|site)))

diag.plots(shannon_gls, col.nos = c(2:4), data = alpha_metric_analysis$data[[2]])

lmerTest::lmer(value ~ health + (1 | site) + (1 | timepoint),
               data = alpha_metric_analysis$data[[2]]) %>% anova(ddf = 'Kenward-Roger')

#Best is no transformation


## Simpson ##
simpson_aov <- aov(value ~ health * timepoint * site, data = alpha_metric_analysis$data[[1]])
diag.plots(simpson_aov, col.nos = c(2:4), data = alpha_metric_analysis$data[[1]])

simpson_gls <- gls(value ~ health * site * timepoint, data = alpha_metric_analysis$data[[1]],
                   weights = varComb(varIdent(form = ~1|health), 
                                     varIdent(form = ~1|timepoint), 
                                     varIdent(form = ~1|site)))

diag.plots(simpson_gls, col.nos = c(2:4), data = alpha_metric_analysis$data[[1]])


simpson_gls_log <- gls(log_value ~ health * timepoint * site, data = alpha_metric_analysis$data[[1]],
                   weights = varComb(varIdent(form = ~1|health), 
                                     varIdent(form = ~1|timepoint), 
                                     varIdent(form = ~1|site)))
diag.plots(simpson_gls_log, col.nos = c(2:4), data = alpha_metric_analysis$data[[1]])

lmerTest::lmer(value ~ health * timepoint + (1 | site),
               data = alpha_metric_analysis$data[[1]]) %>% anova(ddf = 'Kenward-Roger')

#### Plot/Outputs ####
bind_rows(
  richness = as_tibble(anova(richness_gls, type = 'marginal'), rownames = 'term'),
  shannon = as_tibble(anova(shannon_gls, type = 'marginal'), rownames = 'term'),
  simpson = as_tibble(anova(simpson_gls, type = 'marginal'), rownames = 'term'),
  .id = 'metric'
) %>%
  mutate(`F-value` = sprintf(`F-value`, fmt = '%#.3f'),
         `p-value` = scales::pvalue(`p-value`)) %>%
  write_csv('../Results/table2_alpha_diversity_anova.csv')

significant_metrics <- bind_rows(
  Richness = emmeans(richness_gls, ~health | timepoint * site) %>%
    contrast('pairwise') %>%
    as_tibble(),
  
  Shannon = emmeans(shannon_gls, ~health | timepoint * site) %>%
    contrast('pairwise') %>%
    as_tibble(),
  
  `Inverse Simpson` = emmeans(simpson_gls, ~health | timepoint * site) %>%
    contrast('pairwise') %>%
    as_tibble(),
  
  .id = 'metric'
) %>%
  filter(metric != 'Richness',
         p.value < 0.05) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '\n')) %>%
  select(metric, timepoint, site)

bind_rows(
  Richness = emmeans(richness_gls, ~health | timepoint * site) %>%
    as_tibble(),
  
  `Shannon` = emmeans(shannon_gls, ~health | timepoint * site) %>%
    as_tibble(),
  
  `Inverse Simpson` = emmeans(simpson_gls, ~health | timepoint * site) %>%
    as_tibble(),
  
  .id = 'metric'
) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '\n'),
         health = str_replace_all(health, c('D' = 'Diseased', 'H' = 'Healthy'))) %>%
  ggplot(aes(x = site, y = emmean, ymin = lower.CL, ymax = upper.CL,
             colour = health)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_text(data = significant_metrics,
            aes(x = site, y = Inf, label = '*'), 
            size = 10, inherit.aes = FALSE,
            vjust = 1) +
  facet_grid(metric ~ timepoint, scales = 'free_y', 
             switch = 'y') +
  labs(x = 'Site',
       y = NULL,
       colour = 'Health\nState') +
  theme_classic() +
  theme(strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 14),
        panel.background = element_rect(colour = 'black'),
        # legend.position = c(0.9, 0.9),
        legend.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 12))
ggsave('../Results/Fig3_alpha_diversity.png', height = 10, width = 10)


#### PERMANOVA ####
if(file.exists('../intermediate_files/field_other_adonis.rds.gz')){
  the_adonis <- read_rds('../intermediate_files/field_adonis.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  
  adon_y_data <- column_to_rownames(field_data_cpm, 'sample_id')
  
  the_adonis <- adonis2(vegdist(adon_y_data, binary = FALSE, method = 'bray') ~ 
                          health * site * timepoint, 
                        permutations = 9999, by = 'terms',
                        data = field_metadata, parallel = clus)
  
  stopCluster(cl = clus)
  write_rds(the_adonis, '../intermediate_files/field_adonis.rds.gz')
}


the_adonis

#### VEGAN ####
if(file.exists('../intermediate_files/field_tank_nmds.rds.gz')){
  the_nmds <- read_rds('../intermediate_files/field_tank_nmds.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  the_nmds <- column_to_rownames(field_data_cpm, 'sample_id') %>%
    metaMDS(distance = 'bray', binary = FALSE, trymax = 1000, 
            parallel = clus)
  stopCluster(cl = clus)
  write_rds(the_nmds, '../intermediate_files/field_tank_nmds.rds.gz')
}



if(file.exists('../intermediate_files/field_tank_asvArrows.rds.gz')){
  asv_fit <- read_rds('../intermediate_files/field_tank_asvArrows.rds.gz')
} else {
  asv_fit <- field_data %>%
    select(-sample_id) %>%
    envfit(the_nmds, env = .,
           permutations = 9999)
  write_rds(asv_fit, '../intermediate_files/field_tank_asvArrows.rds.gz')
}


#### Plot ####
colony_points_nmds <- scores(the_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(field_data %>%
              select(sample_id, health),
            by = 'sample_id') %>%
  mutate(tank_field = if_else(str_detect(sample_id, 'Bin'),
                              'tank', 'field')) %>%
  # left_join(tank_metadata, by = 'sample_id') %>%
  left_join(field_metadata, by = 'sample_id')

filter(colony_points_nmds, NMDS1 > 0.8)

health_plot <- colony_points_nmds %>% 
  mutate(health = if_else(health == 'D', 'Diseased', 'Healthy')) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = health)) +
  labs(colour = 'Health\nState') +
  theme_classic() +
  theme(legend.position = 'right',
        panel.background = element_rect(colour = 'black'))
# ggsave('../Results/nmds_field_tank.png', height = 7, width = 7)


#### Figure out why ####

year_plot <- colony_points_nmds %>% 
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = as.character(year))) +
  labs(colour = 'Year') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'))

season_plot <- colony_points_nmds %>% 
  mutate(season = if_else(season == 'W', 'January', 'July')) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = season)) +
  labs(colour = 'Season') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'))

site_plot <- colony_points_nmds %>% 
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = site)) +
  labs(colour = 'Site') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'))

timepoint_plot <- colony_points_nmds %>% 
  mutate(timepoint = str_c(season, year, sep = ' '),
         timepoint = str_replace_all(timepoint, c('S' = 'Jul',
                                                  'W' = 'Jan')),
         timepoint = factor(timepoint, levels = c('Jan 2016',
                                                  'Jul 2016',
                                                  'Jan 2017',
                                                  'Jul 2017'))) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = timepoint)) +
  labs(colour = 'Timepoint') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'))

health_plot / 
  (site_plot | timepoint_plot) & 
  plot_annotation(tag_levels = 'A')
ggsave('../Results/Fig4_NMDS_field.png', height = 7, width = 7)

