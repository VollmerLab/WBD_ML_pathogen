#### Libraries ####
library(phyloseq)
library(tidyverse)
library(vegan)

#### N ASVs ####
alpha <- 0.05
rerun_models <- FALSE

#### Data ####
the_nmds <- read_rds('../../intermediate_files/field_tank_nmds.rds.gz')
asv_fit <- read_rds('../../intermediate_files/field_tank_asvArrows.rds.gz')

taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz',
                     show_col_types = FALSE) %>%
  mutate(across(where(is.character), 
                ~str_replace_na(., replacement = '')))

shap_importance <- read_csv('../Results/asv_importance.csv.gz',
                            show_col_types = FALSE) %>%
  filter(p_adjust < alpha) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(asv_id = fct_reorder(asv_id, median_rank)) %>%
  select(asv_id, domain:genus) %>%
  distinct %>%
  arrange(asv_id)


field_data <- read_csv('../../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE)

#### NMDS ####
species_fits <- scores(the_nmds)$species %>%
  as_tibble(rownames = 'asv_id') %>%
  left_join(tibble(asv_id = names(asv_fit$vectors$r),
                   r2 = asv_fit$vectors$r,
                   p = asv_fit$vectors$pvals),
            by = 'asv_id') %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(useful = str_c(family, genus, asv_id, sep = '_') %>%
           str_replace_all('__', '_')) %>% 
  filter(family %in% unique(shap_importance$family)) %>%
  mutate(important_asv = asv_id %in% unique(shap_importance$asv_id)) %>%
  identity()


species_fits %>%
  arrange(important_asv) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_segment(aes(colour = important_asv),
               xend = 0, yend = 0, show.legend = FALSE) +
  geom_point(aes(colour = important_asv, size = important_asv)) +
  geom_point(data = scores(the_nmds)$sites %>%
               as_tibble(rownames = 'sample_id') %>%
               left_join(select(field_data, sample_id, health) %>%
                           distinct(),
                         by = 'sample_id'),
             aes(fill = health),
             size = 2, shape = 21) +
  geom_text(data = . %>% filter(important_asv),
            aes(label = useful), show.legend = FALSE) +
  facet_wrap(~family) +
  scale_colour_manual(values = c('TRUE' = 'black', 'FALSE' = 'grey50')) +
  scale_size_manual(values = c('TRUE' = 5, 'FALSE' = 1)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = 'black')),
         size = 'none') +
  labs(colour = 'ASV Importance',
       fill = 'Health State') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'))
ggsave('../../Results/nmds_important_asvs.png', height = 15, width = 15)
