library(multcomp)
library(multcompView)
library(tidyverse)
library(magrittr)
library(patchwork)
library(vegan)
library(phyloseq)
library(lubridate)

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
  select(sample_id, year, season, site, health, lib.size) %>%
  distinct %>%
  mutate(timepoint = str_c(year, season, sep = '_'),
         across(where(is.character), as.factor)) %>%
  select(sample_id, site, health, timepoint, lib.size)

#### Alpha Diversity Metrics ####
microbiome_data <- read_rds("../intermediate_files/prepped_microbiome.rds.gz")

# raw_counts <- otu_table(microbiome_data) %>%
#   as.data.frame %>%
#   as_tibble(rownames = 'sample_id') %>%
#   mutate(the_asvs = cbind(across(starts_with('ASV'))),
#          .keep = 'unused')

alpha_metrics <- microbiome_data %>%
  # subset_samples(is.na(tank)) %>%
  # filter_taxa(function (x) {sum(x > 0) > 1}, prune = TRUE) %>%
  # subset_taxa(rownames(tax_table(microbiome_data)) %in% colnames(field_data_cpm)) %>%
  microbiome::alpha(zeroes = FALSE) %>%
  as_tibble(rownames = 'sample_id') 
write_csv(alpha_metrics, '../intermediate_files/alpha_metrics.csv')

alpha_metric_analysis <- alpha_metrics %>%
  right_join(field_metadata,
             by = 'sample_id') %>%
  select(sample_id, site, health, timepoint, lib.size, observed, 
         diversity_shannon, diversity_inverse_simpson) %>%
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

library(nlme)
library(emmeans)
source('~/R/R Functions/diagPlots.R')
## Richness ##
richness_aov <- lme(value ~ health * timepoint * lib.size, 
                    random = ~1 | site,
                    data = alpha_metric_analysis$data[[3]])
diag.plots(richness_aov, col.nos = c(2:4), data = alpha_metric_analysis$data[[3]])


#Fix site/timepoint variance
richness_gls <- lme(value ~ health * timepoint, 
                    random = ~1 | site, 
                    data = alpha_metric_analysis$data[[3]],
                    weights = varComb(varIdent(form = ~1|health), 
                                      varIdent(form = ~1|timepoint), 
                                      varIdent(form = ~1|site)))

diag.plots(richness_gls, col.nos = c(2:4), data = alpha_metric_analysis$data[[3]])
car::Anova(richness_gls, type = '2')

write_csv(alpha_metric_analysis$data[[3]], '~/../Desktop/test_data.csv')

#Try log 
richness_gls_sqrt <- lme(sqrt_value ~ health * timepoint, 
                    random = ~1 | site, 
                    data = alpha_metric_analysis$data[[3]],
                    weights = varComb(varIdent(form = ~1|health), 
                                      varIdent(form = ~1|timepoint), 
                                      varIdent(form = ~1|site)))
diag.plots(richness_gls_sqrt, col.nos = c(2:4), data = alpha_metric_analysis$data[[3]])
anova(richness_gls_sqrt)
car::Anova(richness_gls_sqrt, type = '2')


## Shannon ##
shannon_aov <- lme(value ~ health * timepoint, 
                   random = ~1 | site, 
                   data = alpha_metric_analysis$data[[2]])
diag.plots(shannon_aov, col.nos = c(2:4), data = alpha_metric_analysis$data[[2]])

shannon_gls <- lme(value ~ health * timepoint, 
                   random = ~1 | site, 
                   data = alpha_metric_analysis$data[[2]],
                    weights = varComb(varIdent(form = ~1|health), 
                                      varIdent(form = ~1|timepoint), 
                                      varIdent(form = ~1|site)))

diag.plots(shannon_gls, col.nos = c(2:4), data = alpha_metric_analysis$data[[2]])


#Best is no transformation


## Simpson ##
simpson_aov <- lme(value ~ health * timepoint, 
                   random = ~1 | site, 
                   data = alpha_metric_analysis$data[[1]])
diag.plots(simpson_aov, col.nos = c(2:4), data = alpha_metric_analysis$data[[1]])


simpson_gls_log <- lme(log_value ~ health * timepoint, 
                       random = ~1 | site,
                       data = alpha_metric_analysis$data[[1]],
                   weights = varComb(varIdent(form = ~1|health), 
                                     varIdent(form = ~1|timepoint), 
                                     varIdent(form = ~1|site)))
diag.plots(simpson_gls_log, col.nos = c(2:4), data = alpha_metric_analysis$data[[1]])

lmerTest::lmer(value ~ health * timepoint + (1 | site),
               data = alpha_metric_analysis$data[[1]]) %>% anova(ddf = 'Kenward-Roger')

#### Plot/Outputs ####
as_tibble(car::Anova(simpson_gls_log, type = '2'), rownames = 'term')

bind_rows(
  richness = as_tibble(car::Anova(richness_gls_sqrt, type = '2'), rownames = 'term'),
  shannon = as_tibble(car::Anova(shannon_gls, type = '2'), rownames = 'term'),
  simpson = as_tibble(car::Anova(simpson_gls_log, type = '2'), rownames = 'term'),
  .id = 'metric'
) %>%
  mutate(`Chisq` = sprintf(Chisq, fmt = '%#.3f'),
         `p-value` = scales::pvalue(`Pr(>Chisq)`),
         .keep = 'unused') %>%
  filter(term != '(Intercept)') %>%
  write_csv('../Results/table2_alpha_diversity_anova.csv')


significant_metrics <- bind_rows(
  Richness = emmeans(richness_gls_sqrt, ~health | timepoint) %>%
    contrast('pairwise') %>%
    as_tibble(),
  
  Shannon = emmeans(shannon_gls, ~health | timepoint) %>%
    contrast('pairwise') %>%
    as_tibble(),
  
  `Inverse Simpson` = emmeans(simpson_gls_log, ~health | timepoint) %>%
    contrast('pairwise') %>%
    as_tibble(),
  
  .id = 'metric'
) %>%
  filter(p.value < 0.05) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '-'),
         timepoint = ymd(str_c(timepoint, '1', sep = '-'))) %>%
  select(metric, timepoint)

significant_metrics_time <- bind_rows(
  Richness = emmeans(richness_gls_sqrt, ~timepoint) %>%
    cld(Letters = LETTERS) %>%
    as_tibble(),
  
  Shannon = emmeans(shannon_gls, ~timepoint) %>%
    cld(Letters = LETTERS) %>%
    as_tibble(),
  
  `Inverse Simpson` = emmeans(simpson_gls_log, ~timepoint) %>%
    cld(Letters = LETTERS) %>%
    as_tibble(),
  
  .id = 'metric'
) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '-'),
         timepoint = ymd(str_c(timepoint, '1', sep = '-')),
         .group = str_trim(.group)) %>%
  select(metric, timepoint, .group)



ref_grid(richness_gls_sqrt) %>%
  update(tran = poisson(link = 'sqrt')) %>%
  regrid(transform = 'response') %>%
  emmeans(~health) %>%
  contrast('pairwise')

bind_rows(
  Richness = ref_grid(richness_gls_sqrt) %>%
    update(tran = poisson(link = 'sqrt')) %>%
    emmeans(~health | timepoint, type = 'response') %>%
    as_tibble() %>%
    rename(emmean = response),
  
  `Shannon` = emmeans(shannon_gls, ~health | timepoint) %>%
    as_tibble(),
  
  `Inverse Simpson` = ref_grid(simpson_gls_log) %>%
    update(tran = poisson(link = 'log')) %>%
    emmeans(~health | timepoint, type = 'response') %>%
    as_tibble() %>%
    rename(emmean = response),
  
  .id = 'metric'
) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '-'),
         timepoint = ymd(str_c(timepoint, '1', sep = '-')),
         health = str_replace_all(health, c('D' = 'Diseased', 'H' = 'Healthy'))) %>%
  ggplot(aes(x = timepoint, y = emmean, ymin = lower.CL, ymax = upper.CL,
             colour = health)) +
  geom_pointrange(position = position_dodge2(50)) +
  geom_text(data = significant_metrics,
            aes(x = timepoint, y = Inf, label = '*'), 
            size = 10, inherit.aes = FALSE,
            vjust = 1) +
  scale_colour_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, type = "continuous"),
                                         c('Healthy', 'Diseased'))) +
  # facet_wrap(~metric, scales = 'free_y')
  scale_x_date(breaks = ymd(c('2016-01-01', '2016-07-01', 
                              '2017-01-01', '2017-07-01')), 
               date_labels = '%Y\n%b') +
  facet_grid(metric ~ ., scales = 'free_y', 
             switch = 'y') +
  labs(x = NULL,
       y = NULL,
       colour = 'Disease\nState') +
  theme_classic() +
  theme(strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 14),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank(),
        # legend.position = c(0.9, 0.9),
        legend.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 12))
ggsave('../Results/alpha_diversity.png', height = 7, width = 7)


bind_rows(
  Richness = ref_grid(richness_gls_sqrt) %>%
    update(tran = poisson(link = 'sqrt')) %>%
    emmeans(~timepoint, type = 'response') %>%
    as_tibble() %>%
    rename(emmean = response),
  
  `Shannon` = emmeans(shannon_gls, ~timepoint) %>%
    as_tibble(),
  
  `Inverse Simpson` = ref_grid(simpson_gls_log) %>%
    update(tran = poisson(link = 'log')) %>%
    emmeans(~timepoint, type = 'response') %>%
    as_tibble() %>%
    rename(emmean = response),
  
  .id = 'metric'
) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '-'),
         timepoint = ymd(str_c(timepoint, '1', sep = '-'))) %>% filter(timepoint == ymd('2017-07-01'))
  left_join(significant_metrics_time, by = c('timepoint', 'metric')) %>%
  ggplot(aes(x = timepoint, y = emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange(position = position_dodge2(50)) +
  geom_text(aes(x = timepoint, y = Inf, label = .group), 
            size = 10, inherit.aes = FALSE,
            vjust = 1) +
  scale_colour_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, type = "continuous"),
                                         c('Healthy', 'Diseased'))) +
  # facet_wrap(~metric, scales = 'free_y')
  scale_x_date(breaks = ymd(c('2016-01-01', '2016-07-01', 
                              '2017-01-01', '2017-07-01')), 
               date_labels = '%Y\n%b') +
  facet_grid(metric ~ ., scales = 'free_y', 
             switch = 'y') +
  labs(x = NULL,
       y = NULL,
       colour = 'Disease\nState') +
  theme_classic() +
  theme(strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 14),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank(),
        # legend.position = c(0.9, 0.9),
        legend.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 12))


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
  asv_fit <- field_data_cpm %>%
    select(-sample_id) %>%
    envfit(the_nmds, env = .,
           permutations = 9999)
  write_rds(asv_fit, '../intermediate_files/field_tank_asvArrows.rds.gz')
}


#### Plot ####
colony_points_nmds <- scores(the_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(field_metadata, by = 'sample_id')

filter(colony_points_nmds, NMDS1 > 0.8)

health_plot <- colony_points_nmds %>% 
  mutate(health = if_else(health == 'D', 'Diseased', 'Healthy')) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = health, shape = health)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(colour = 'Disease\nState',
       shape = 'Disease\nState') +
  theme_classic() +
  theme(legend.position = 'right',
        panel.background = element_rect(colour = 'black'),
        strip.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10))
# ggsave('../Results/nmds_field_tank.png', height = 7, width = 7)

timepoint_plot <- colony_points_nmds %>% 
  mutate(timepoint = str_replace_all(timepoint, c('W' = 'Jan', 'S' = 'Jul'))) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(data=. %>% select(-timepoint), 
             colour="grey60", size = 0.75) +
  geom_point(aes(colour = site, shape = health), size = 1.5) +
  facet_wrap(~timepoint, labeller = labeller(timepoint = ~str_replace(., '_', ' - '))) +
  guides(shape = 'none') +
  labs(colour = 'Site') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10))

health_plot / 
  (timepoint_plot) & 
  plot_annotation(tag_levels = 'A')
# ggsave('../Results/Fig4_NMDS_field.png', height = 7, width = 7)


colony_points_nmds %>% 
  mutate(timepoint = str_replace_all(timepoint, c('W' = 'Jan', 'S' = 'Jul')),
         health = if_else(health == 'D', 'Diseased', 'Healthy')) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(shape = site, colour = health, fill = health), size = 1.5) +
  facet_wrap(~timepoint, labeller = labeller(timepoint = ~str_replace(., '_', ' - '))) +
  scale_shape_manual(values = c('CK14' = 'circle filled', 'CK4' = 'square filled', 'HS' = 'diamond filled', 
                                'SB' = 'triangle filled', 'TS' = 'triangle down filled')) +
  scale_colour_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, type = "continuous"),
                                         c('Healthy', 'Diseased'))) +
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, type = "continuous"),
                                         c('Healthy', 'Diseased'))) +
  guides(fill = 'none',
         shape = guide_legend(override.aes = list(fill = 'black', size = 3)),
         colour = guide_legend(override.aes = list(size = 3))) +
  labs(colour = 'Disease\nState',
       shape = 'Site') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10))
ggsave('../Results/Fig3_NMDS_field.png', height = 7, width = 7)

colony_points_nmds %>% 
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(data=. %>% select(-site), 
             colour="grey") +
  geom_point(aes(colour = timepoint)) +
  facet_wrap(~site) +
  labs(colour = 'Site') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'))
