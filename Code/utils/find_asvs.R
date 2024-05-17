#### Libraries ####
library(phyloseq)
library(tidyverse)

genera_to_find <- c('MD3-55', 'Sphingobium')

#### Data ####
# microbiome_raw <- read_rds("../Data/field_tank_newPS_deciphersilva.rds") 

microbiome_analysis <- read_rds('../../intermediate_files/prepped_microbiome.rds.gz')

field_data <- read_csv('../../intermediate_files/normalized_field_asv_counts.csv',
                       show_col_types = FALSE)
tank_data <- read_csv('../../intermediate_files/normalized_tank_asv_counts.csv',
                      show_col_types = FALSE)

#### Find ASVs ####
asvs_interest <- tax_table(microbiome_analysis) %>%
  as.data.frame %>%
  as_tibble(rownames = 'asv_id') %>%
  filter(genus %in% genera_to_find)

asvs_interest %>%
  count(genus, species)
 

filter(field_data, asv_id %in% asvs_interest$asv_id) 
filter(tank_data, asv_id %in% asvs_interest$asv_id) 


all_asvs_search <- otu_table(microbiome_analysis) %>%
  as.data.frame %>%
  as_tibble(rownames = 'sample_id') %>%
  select(sample_id, all_of(asvs_interest$asv_id))




asv_presence <- all_asvs_search %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = 'asv_id',
               values_to = 'n_reads') %>%
  left_join(asvs_interest, by = 'asv_id') %>%
  left_join(sample_data(microbiome_analysis) %>%
              as_tibble(rownames = 'sample_id') %>%
              select(sample_id, health, dataset),
            by = 'sample_id') %>%
  group_by(genus, asv_id) %>%
  summarise(n_sample = sum(n_reads > 0),
            n_field = sum(n_reads > 0 & dataset == 'field'),
            n_tank = sum(n_reads > 0 & dataset == 'tank'),
            n_healthy = sum(n_reads > 0 & health == 'H'),
            n_disease = sum(n_reads > 0 & health == 'D'),
            total_sample = n(),
            .groups = 'drop') 

write_csv(asv_presence, '../../Results/aquarick_shig_counts.csv')

