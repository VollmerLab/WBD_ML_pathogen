library(tidyverse)
library(magrittr)
library(patchwork)
library(vegan)

#### Data ####
taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz', show_col_types = FALSE) %>%
  mutate(across(everything(), str_replace_na))

field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  select(asv_id, sample_id, n_reads, health) %>%
  pivot_wider(names_from = asv_id, values_from = n_reads) %>%
  mutate(health = factor(health))

field_metadata <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                           show_col_types = FALSE) %>%
  select(sample_id, year, season, site) %>%
  distinct


tank_raw <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                     show_col_types = FALSE) %>%
  select(-site, -dataset, -year, -season, -anti, -health) %>%
  select(-domain:-species) %>%
  rename(time_treat = time) %>%
  mutate(time = str_extract(time_treat, '[0-9]') %>% str_c('T', .)) %>%
  filter(geno != 'GE') %>%
  separate(anti_health, into = c('with_anti', 'disease_exposure')) %>%
  rename(health = resist) %>%
  mutate(disease_exposure = if_else(time == 'T0', 'N', disease_exposure),
         health = str_replace_all(health, c('R' = 'H', 'S' = 'D')) %>% factor()) %>%
  filter(with_anti == 'N', (disease_exposure == 'D' | time == 'T0')) %>%
  mutate(sample_id = str_c(sample_id, time, sep = '.'), .keep = 'unused') %>%
  # filter(str_detect(sample_id, 'T8')) %>%
  identity()

tank_data <- tank_raw %>%
  select(asv_id, health, sample_id, n_reads) %>%
  pivot_wider(names_from = asv_id, values_from = n_reads) 

tank_metadata <- tank_raw %>%
  select(sample_id, geno, tank) %>%
  distinct %>%
  mutate(time = str_extract(sample_id, 'T[0-9]+$'))

#### PERMANOVA ####
adon_y_data <- bind_rows(field_data,
          tank_data) %>%
  filter(!str_detect(sample_id, 'Bin')) %>%
  mutate(across(where(is.numeric), ~as.integer(. > 0))) %>%
  select(-health) %>%
  column_to_rownames('sample_id')

adon_x_data <- full_join(field_metadata, 
          select(field_data, sample_id, health),
          by = 'sample_id') %>%
  mutate(year = as.character(year),
         across(where(is.character), as.factor))

if(file.exists('../intermediate_files/field_other_adonis.rds.gz')){
  the_health_adonis <- read_rds('../intermediate_files/field_health_adonis.rds.gz')
  the_other_adonis <- read_rds('../intermediate_files/field_other_adonis.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  
  the_health_adonis <- adonis2(vegdist(adon_y_data, binary = TRUE, method = 'raup') ~ 
                          health, 
                        permutations = 9999, by = 'terms',
                        data = adon_x_data, parallel = clus)
  
  the_other_adonis <- adonis2(vegdist(adon_y_data, binary = TRUE, method = 'raup') ~ 
                                 season + site + year, 
                               permutations = 9999, by = 'margin',
                               data = adon_x_data, parallel = clus)
  stopCluster(cl = clus)
  write_rds(the_health_adonis, '../intermediate_files/field_health_adonis.rds.gz')
  write_rds(the_other_adonis, '../intermediate_files/field_other_adonis.rds.gz')
}


the_health_adonis
the_other_adonis


#### VEGAN ####
if(file.exists('../intermediate_files/field_tank_nmds.rds.gz')){
  the_nmds <- read_rds('../intermediate_files/field_tank_nmds.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  the_nmds <- bind_rows(field_data,
                        tank_data) %>%
    filter(!str_detect(sample_id, 'Bin')) %>%
    mutate(across(where(is.numeric), ~as.integer(. > 0))) %>%
    select(-health) %>%
    column_to_rownames('sample_id') %>%
    #remove identical points
    metaMDS(distance = 'raup', binary = TRUE, trymax = 1000, 
            parallel = clus)
  stopCluster(cl = clus)
  write_rds(the_nmds, '../intermediate_files/field_tank_nmds.rds.gz')
}



if(file.exists('../intermediate_files/field_tank_asvArrows.rds.gz')){
  asv_fit <- read_rds('../intermediate_files/field_tank_asvArrows.rds.gz')
} else {
  asv_fit <- bind_rows(field_data,
                       tank_data) %>%
    filter(!str_detect(sample_id, 'Bin')) %>%
    select(-sample_id) %>%
    envfit(the_nmds, env = .,
           permutations = 9999)
  write_rds(asv_fit, '../intermediate_files/field_tank_asvArrows.rds.gz')
}


#### Plot ####
colony_points_nmds <- scores(the_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(bind_rows(field_data,
                      tank_data) %>%
              select(sample_id, health),
            by = 'sample_id') %>%
  mutate(tank_field = if_else(str_detect(sample_id, 'Bin'),
                              'tank', 'field')) %>%
  left_join(tank_metadata, by = 'sample_id') %>%
  left_join(field_metadata, by = 'sample_id')

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

health_plot / 
  (site_plot | year_plot | season_plot) & 
  plot_annotation(tag_levels = 'A')
ggsave('../Results/Fig2_NMDS_field.png', height = 7, width = 15)

