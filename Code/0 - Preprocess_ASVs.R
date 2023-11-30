#### Libraries ####

library(tidyverse)
library(magrittr)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(edgeR)
library(ggvenn)
library(ggupset)


#### Functions ####
# data <- otu_tmm
filter_missingness <- function(data, prop_missing){
  #data = a DGElist object, model_sample = samples to use for the calculation of % missingness
  #prop_missing = maximum percentage of samples which can be 0 and still keep ASV in dataset
  keep <- rowMeans(data$counts == 0) <= prop_missing
  
  message('\n')
  message('ASVs Removed for not being expressed in enough samples: ', scales::comma(sum(!keep)))
  message('ASVs Kept by Filter: ', scales::comma(sum(keep)))
  message('\n')
  data[keep, keep.lib.sizes = FALSE]
}

#### Data ####
microbiome_raw <- read_rds("../Data/field_tank_newPS_deciphersilva.rds") 

microbiome_data <- microbiome_raw %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  subset_taxa(domain == "Bacteria" &
                phylum != "Cyanobacteria" &
                !is.na(phylum) &
                family != "Mitochondria" &
                class != "Chloroplast" & 
                order != "Chloroplast") %>%
  prune_samples(sample_sums(.) > 1000, .) %>%
  subset_samples(site != 'MG') #no disease samples

metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  mutate(health = if_else(health == 'N', 'H', health))

count(metadata, dataset, year, season)
count(metadata, dataset, health)
count(metadata, dataset, site)

count(metadata, dataset)

filter(metadata, dataset == 'tank',
       anti != 'A') %>% 
  count(tank, anti_health)

#### Sort out Experimental Design ####
sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  filter(dataset == 'tank') %>%
  select(time, tank, geno, anti, health, anti_health, resist, sample_id) %>%
  count(anti, health, anti_health, time)

#### Split Field & Tank Metadata ####
field_meta <- metadata %>% 
  filter(dataset == 'field')

tank_meta <- metadata %>%
  filter(dataset == 'tank')

#### Check overlap btw tank & field across years ####
sample_point_asvs <- phyloseq_filter_prevalence(microbiome_data, prev.trh = 0.1) %>%
  otu_table() %>%
  as.data.frame() %>%
  as_tibble(rownames = 'sample_id') %>%
  mutate(sample_id = str_remove(sample_id, '^X')) %>%
  left_join(metadata,
            by = 'sample_id') %>%
  mutate(across(starts_with('ASV'), ~. > 0)) %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = 'asv_id',
               values_to = 'pa') %>%
  filter(pa) %>%
  select(-pa) %>%
  select(sample_id, dataset, year, season, asv_id) %>%
  group_by(dataset, year, season, asv_id) %>%
  summarise(n_samples = n_distinct(sample_id),
            .groups = 'drop') %>%
  mutate(data_year_season = str_c(dataset, year, season, sep = '_'), 
         .keep = 'unused') %>%
  group_by(asv_id) %>%
  summarise(n_sample_points = n_distinct(data_year_season),
            sample_points = list(data_year_season)) 

sample_point_asvs %>%
  ggplot(aes(x = sample_points)) +
  geom_bar() +
  scale_y_continuous() +
  scale_x_upset() +
  theme_classic() +
  theme_combmatrix(combmatrix.label.make_space = TRUE)

asvs_to_keep <- sample_point_asvs %>%
  filter(n_sample_points == 5) %>%
  pull(asv_id)
length(asvs_to_keep)

#### Check preliminary clustering ####
# library(vegan)
# prelim_nmds <- microbiome_data %>%
#   subset_taxa(rownames(tax_table(microbiome_data)) %in% asvs_to_keep) %>%
#   subset_samples(dataset == 'field') %>%
#   otu_table() %>%
#   as.data.frame() %>%
#   metaMDS()
# 
# scores(prelim_nmds, 'sites') %>%
#   as_tibble(rownames = 'sample_id') %>%
#   left_join(metadata, by = 'sample_id') %>%
#   
#   ggplot(aes(x = NMDS1, y = NMDS2, colour = interaction(year, season), shape = health)) +
#   geom_point() +
#   facet_wrap(~site)
#Remove MG because it has no disease samples

#### Normalize ASV counts ####
otu_tmm <- microbiome_data %>%
  subset_taxa(rownames(tax_table(microbiome_data)) %in% asvs_to_keep) %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  
  #Add any other filtering here
  filter_missingness(0.9) %>%
  
  edgeR::calcNormFactors(method = 'TMMwsp')

#### Output CPM ####
full_data <- full_join(cpm(otu_tmm, log = TRUE, prior.count = 0.5,
              normalized.lib.sizes = TRUE) %>%
            as_tibble(rownames = 'asv_id') %>%
            pivot_longer(cols = -asv_id,
                         names_to = 'sample_id',
                         values_to = 'log2_cpm_norm'),
          
          cpm(otu_tmm, log = FALSE, 
              normalized.lib.sizes = TRUE) %>%
            as_tibble(rownames = 'asv_id') %>%
            pivot_longer(cols = -asv_id,
                         names_to = 'sample_id',
                         values_to = 'cpm_norm'),
          by = c('asv_id', 'sample_id')) %>%
  
  full_join(cpm(otu_tmm, log = FALSE,
                normalized.lib.sizes = FALSE) %>%
              as_tibble(rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'cpm'),
            by = c('asv_id', 'sample_id')) %>%
  full_join(as_tibble(otu_tmm$counts, rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'n_reads'),
            by = c('asv_id', 'sample_id')) %>%
  
  left_join(metadata, 
            by = 'sample_id') %>%
  left_join(as_tibble(otu_tmm$samples, rownames = 'sample_id') %>%
              select(-group),
            by = 'sample_id') %>%
  left_join(tax_table(microbiome_data) %>%
              as.data.frame() %>%
              as_tibble(rownames = 'asv_id'),
            by = c('asv_id'))


filter(full_data, dataset == 'field') %>%
  janitor::remove_empty(which = 'cols') %>%
  write_csv('../intermediate_files/normalized_field_asv_counts.csv')

filter(full_data, dataset == 'field') %>%
  select(asv_id, sample_id, log2_cpm_norm, health) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm)


  write_csv('../intermediate_files/normalized_field_asv_counts.csv')

filter(full_data, dataset == 'tank') %>%
  write_csv('../intermediate_files/normalized_tank_asv_counts.csv')

