#### Libraries ####

library(tidyverse)
library(magrittr)
library(phyloseq)
library(microbiome)
library(taxize)
library(metagMisc)
library(edgeR)
library(ggvenn)
library(ggupset)

#### Tank samples to remove for double sequencing?? ####
#Definitely remove yellow - maybe remove green/orange
#Yellow has 2 samples and should only have 1, green/orange should both have 1 and do but were diseased prior to first sampling
# remove_samples <- c('Bin1_N_D_YE', 'Bin2_N_D_YE', 'Bin5_N_D_OR', 'Bin5_N_D_GR')
remove_samples <- 'none'

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

top_classification <- function(data, threshold = 80){
  rename_with(data, .cols = domain:species,
              ~str_c(., '_name')) %>%
    pivot_longer(cols = -asv_id,
                 names_to = c('taxon_level', '.value'),
                 names_pattern = '(.*)_(.*)') %>%
    mutate(taxon_level = factor(taxon_level, ordered = TRUE,
                                levels = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'))) %>%
    filter(confidence >= threshold) %>%
    select(-confidence) %>%
    pivot_wider(names_from = taxon_level, 
                values_from = 'name')
}

synonomize_taxonomy <- function(taxonomy, lowest_level = 'genus'){
  all_levels <- c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'asv_id')
  remove_levels <- all_levels[(which(all_levels == lowest_level) + 1):length(all_levels)]
  
  duplicated_higher_taxonomy <- select(taxonomy, -asv_id, -all_of(remove_levels)) %>%
    distinct %>%
    filter(!is.na(!!sym(lowest_level))) %>%
    group_by(!!sym(lowest_level)) %>%
    filter(n() > 1) %>%
    arrange(!!sym(lowest_level)) %>%
    pull(!!sym(lowest_level)) %>%
    unique
  
  if(lowest_level == 'species'){
    genus_species_mismatch <- select(taxonomy, -asv_id, -all_of(remove_levels)) %>%
      distinct %>%
      filter(!is.na(!!sym(lowest_level))) %>%
      filter(!str_detect(species, genus)) %>%
      pull(species) %>%
      unique %>%
      sort
    
    duplicated_higher_taxonomy <- unique(c(duplicated_higher_taxonomy, genus_species_mismatch))
  }
  
  if(length(duplicated_higher_taxonomy) > 0){
    update_classification <- classification(duplicated_higher_taxonomy, db = 'ncbi')
    
    refresh_taxa <- tibble(!!sym(lowest_level) := names(update_classification)) %>%
      mutate(row = row_number()) %>%
      rowwise(!!sym(lowest_level)) %>%
      reframe(update_classification[[row]]) %>%
      filter(rank %in% all_levels,
             rank != lowest_level) %>%
      select(-id) %>%
      pivot_wider(names_from = rank, 
                  values_from = name) %>%
      mutate(domain = 'Bacteria')
    
    
    
    out <- filter(taxonomy, !(!!sym(lowest_level) %in% refresh_taxa[[lowest_level]])) %>%
      bind_rows(filter(taxonomy, !!sym(lowest_level) %in% refresh_taxa[[lowest_level]]) %>%
                  select(asv_id, !!sym(lowest_level), all_of(remove_levels)) %>%
                  left_join(refresh_taxa, by = lowest_level))
  } else {
    out <- taxonomy
  }
  out
}

phyloseq_filter_prevalence_fix <- function (physeq, prev.trh = 0.05, 
                                            abund.trh = NULL, threshold_condition = "OR", 
                                            abund.type = "total") {
  if (prev.trh > 1 | prev.trh < 0) {
    stop("Prevalence threshold should be non-negative value in the range of [0, 1].\n")
  }
  if (!is.null(abund.trh)) {
    if (abund.trh <= 0) {
      stop("Abundance threshold should be non-negative value larger 0.\n")
    }
  }
  prevalenceThreshold <- prev.trh * phyloseq::nsamples(physeq)
  prevdf <- prevalence(physeq, package = 'base')
  if (abund.type == "total") {
    prevdf$AbundFilt <- prevdf$TotalAbundance
  }
  if (abund.type == "mean") {
    prevdf$AbundFilt <- prevdf$MeanAbundance
  }
  if (abund.type == "median") {
    prevdf$AbundFilt <- prevdf$MedianAbundance
  }
  if (is.null(abund.trh)) {
    tt <- prevdf$Prevalence >= prevalenceThreshold
  }
  if (!is.null(abund.trh)) {
    if (threshold_condition == "OR") {
      tt <- (prevdf$Prevalence >= prevalenceThreshold | 
               prevdf$AbundFilt >= abund.trh)
    }
    if (threshold_condition == "AND") {
      tt <- (prevdf$Prevalence >= prevalenceThreshold & 
               prevdf$AbundFilt >= abund.trh)
    }
  }
  keepTaxa <- prevdf$Taxa[tt]
  res <- phyloseq::prune_taxa(taxa = keepTaxa, x = physeq)
  return(res)
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
  prune_samples(sample_sums(.) > 1000, .) %>% #for some reason must do twice
  subset_samples(site != 'MG') %>% #no disease samples
  prune_samples(!str_detect(rownames(sample_data(.)), 
                            str_c(remove_samples, collapse = '|')), .)

## Reclassify with Decipher 
if(!file.exists('../Data/update_taxonomy.csv')){
  Biostrings::writeXStringSet(file="../intermediate_files/all_asvs.fasta", refseq(microbiome_data))
} else {
  new_taxonomy <- read_csv('../Data/update_taxonomy.csv', 
                           show_col_types = FALSE) %>%
    top_classification()
}


old_taxonomy <- tax_table(microbiome_data) %>%
  as.data.frame %>%
  as_tibble(rownames = 'asv_id')

updated_taxonomy <- bind_rows(new_taxonomy,
                              filter(old_taxonomy, !asv_id %in% new_taxonomy$asv_id))

reconciled_higher_taxonomy <- synonomize_taxonomy(updated_taxonomy, 'species') %>%
  synonomize_taxonomy('genus') %>%
  synonomize_taxonomy('family') %>%
  synonomize_taxonomy('order') %>%
  synonomize_taxonomy('class') %>%
  synonomize_taxonomy('phylum') 

tax_table(microbiome_data) <- column_to_rownames(reconciled_higher_taxonomy, 'asv_id') %>% 
  as.matrix

metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  mutate(health = if_else(health == 'N', 'H', health),
         geno = if_else(geno == 'GE', 'GR', geno),
         sample_id = str_replace(sample_id, 'GE', 'GR'))


sequenced_tank_samples <- filter(metadata, dataset == 'tank') %>%
  mutate(fragment_id = str_remove(sample_id, 'P[0-9]_')) %>%
  pull(fragment_id) %>%
  unique

#### Read in original Tank Data ####
original_tank <- read_csv('../Data/su17_tank_surv_data.csv', 
         show_col_types = FALSE) %>%
  mutate(fragment_id = str_c(str_c('Bin', Block), if_else(antibiotics == 'no_antibiotics', 
                                                          'N', 'A'),
                             if_else(exposure == 'healthy', 'N', 'D'), 
                             str_sub(genotype, 1, 2) %>% str_to_upper(), sep = '_'),
         .before = everything()) %>%
  filter(fragment_id %in% sequenced_tank_samples) %>%
  filter(!genotype %in% c('Grey', 'Blue')) %>% #Did not sequence these genotypes for 16s
  rename_with(~str_replace_all(., c('21' = 'jul14am',	'30' = 'jul14pm',	
                                    '45' = 'jul15am',	'54' = 'jul15pm',	
                                    '69' = 'jul16am',	'78' = 'jul16pm',	
                                    '93' = 'jul17am',	'102' = 'jul17pm',
                                    '117' = 'jul18am', '126' = 'jul18pm',	
                                    '141' = 'jul19am',	'150' = 'jul19pm',	
                                    '165' = 'jul20am',	'174' = 'jul20pm'))) %>%
  
  pivot_longer(cols = starts_with('jul'),
               names_to = 'date',
               values_to = 'state') %>%
  
  mutate(tank = str_c('Bin', Block),
         health = if_else(state == 0, 'H', 'D'),
         geno = str_sub(genotype, 1, 2) %>% str_to_upper(),
         # date = fct_inorder(date),
         #no time 0_anti in this dataset - must keep the metadata from before - all are healthy
         time = case_when(date %in% c('jul15pm') ~ '2_exp',
                          date %in% c('jul20pm') ~ '8_exp',
                          TRUE ~ 'other'),
         .keep = 'unused') %>%
  filter(time != 'other') %>%
  select(-Condition) %>%
  
  mutate(plate = case_when(time == '2_exp' ~ 'P5',
                           time == '8_exp' ~ 'P6',
                           TRUE ~ NA_character_),
         sample_id = str_c(plate, fragment_id, sep = '_'),
         anti = if_else(antibiotics == 'no_antibiotics', 'N', 'A'),
         exposure = if_else(exposure == 'healthy', 'N', 'D'),
         anti_exposure = str_c(anti, exposure, sep = '_'),
         site = 'tank', season = 'S', dataset = 'tank',
         year = '2017') %>%
  select(sample_id, health, time, tank, anti_exposure, anti,
         exposure, geno, plate, year, season, site, dataset)


#### Swap tank data to original tank data rather than post-processed data ####
pre_exposure <- filter(metadata, dataset == 'tank') %>%
  filter(time == '0_anti') %>%
  select(-resist, -anti_health) %>%
  mutate(exposure = 'pre')

updated_tank_data <- filter(metadata, dataset == 'tank', time != '0_anti') %>%
  select(-health, -resist, -anti_health) %>%
  left_join(select(original_tank, sample_id, health, exposure),
            by = 'sample_id') %>%
  bind_rows(pre_exposure) %>%
  mutate(tank = str_c(tank, anti, exposure, sep = '_'),
         tank = str_remove(tank, '_pre'))

metadata <- select(metadata, -resist, -anti_health) %>%
  filter(dataset != 'tank') %>%
  bind_rows(updated_tank_data)

#### Summary Stats ####
count(metadata, dataset, year, season)
count(metadata, dataset, health)
count(metadata, dataset, health, year, season) %>%
  pivot_wider(names_from = health,
              values_from = n)

count(metadata, dataset)

filter(metadata, dataset == 'tank',
       anti != 'A') %>% 
  count(tank, anti, exposure)

#### Experimental Design ####
metadata %>%
  filter(dataset == 'tank') %>%
  select(time, tank, geno, anti, exposure, health, sample_id) %>%
  count(anti, exposure, health, time)

#### Split Field & Tank Metadata ####
field_meta <- metadata %>% 
  filter(dataset == 'field')

tank_meta <- metadata %>%
  filter(dataset == 'tank')

str_subset(tank_meta$sample_id, 'GE')

#### Check overlap btw tank & field across years ####
sample_point_asvs <- phyloseq_filter_prevalence_fix(microbiome_data, prev.trh = 0.1) %>%
  otu_table() %>%
  as.data.frame() %>%
  as_tibble(rownames = 'sample_id') %>%
  mutate(sample_id = str_remove(sample_id, '^X'),
         sample_id = str_replace(sample_id, 'GE', 'GR')) %>%
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

#### Normalize ASV counts ####
otu_tmm <- microbiome_data %>%
  subset_taxa(rownames(tax_table(microbiome_data)) %in% asvs_to_keep) %>%
  phyloseq_filter_prevalence_fix(prev.trh = 0.1) %>% 
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

#### Output Files for Analysis ####
filter(full_data, dataset == 'field') %>%
  janitor::remove_empty(which = 'cols') %>%
  write_csv('../intermediate_files/normalized_field_asv_counts.csv')

filter(full_data, dataset == 'tank') %>%
  write_csv('../intermediate_files/normalized_tank_asv_counts.csv')


sample_data(microbiome_data) <- column_to_rownames(metadata, 'sample_id')
prepped_microbiome_data <- subset_samples(microbiome_data, is.na(anti) | anti == 'N')

write_rds(prepped_microbiome_data, '../intermediate_files/prepped_microbiome.rds.gz')
