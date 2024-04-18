## 1 - find all cysteniphlum in our data
## 2 - find all cysteniphilum in rosales

#### Libraries ####
library(Biostrings)
library(msa)
library(phyloseq)
library(tidyverse)


#### Functions ####
read_blca <- function(blca_out){
  read_delim(blca_out, delim = '\t', 
             col_names = c('asv_id', 'taxonomy'), show_col_types = FALSE) %>%
    mutate(superkingdom = str_extract(taxonomy, 'superkingdom:.*;[0-9\\.]+;phylum') %>% str_remove_all('superkingdom:|;phylum'),
           phylum = str_extract(taxonomy, 'phylum:.*;[0-9\\.]+;class') %>% str_remove_all('phylum:|;class'),
           class = str_extract(taxonomy, 'class:.*;[0-9\\.]+;order') %>% str_remove_all('class:|;order'),
           order = str_extract(taxonomy, 'order:.*;[0-9\\.]+;family') %>% str_remove_all('order:|;family'),
           family = str_extract(taxonomy, 'family:.*;[0-9\\.]+;genus') %>% str_remove_all('family:|;genus'),
           genus = str_extract(taxonomy, 'genus:.*;[0-9\\.]+;species') %>% str_remove_all('genus:|;species'),
           species = str_extract(taxonomy, 'species:.*;[0-9\\.]+') %>% str_remove_all('species:'),
           .keep = 'unused') %>%
    rename(domain = superkingdom) %>%
    # select(asv_id, superkingdom, phylum) %>%
    
    mutate(across(-asv_id, ~str_extract(., ';.*$') %>% str_remove(';') %>% 
                    parse_number(), .names = '{.col}_confidence'),
           across(where(is.character), ~str_remove(., ';.*$')))
}

#### Data ####
panama_taxonomy <- read_blca('../../intermediate_files/all_asvs.fasta.blca.out')

rosales_taxonomy <- read_blca('../../intermediate_files/rosales_2019.fasta.blca.out')

taxonomy <- bind_rows(ours = panama_taxonomy,
                      rosales = rosales_taxonomy,
                      .id = 'study')

panama_sequences <- read_rds('../../intermediate_files/prepped_microbiome.rds.gz') %>%
  refseq

rosales_sequences <- readDNAStringSet('../../intermediate_files/rosales_2019.fasta')

all_sequences <- c(panama_sequences, rosales_sequences)

#### Check for cysteniphilum ####
the_targets <- filter(taxonomy, genus == 'Cysteiniphilum') %>%
  select(-domain:-family, -domain_confidence:-family_confidence) %>%
  mutate(species = str_remove(species, genus) %>% str_trim,
         general_species = str_c(genus, species, sep = ' '),
         species = str_c(genus, ' (', scales::percent(genus_confidence, scale = 1, accuracy = 1, suffix = ''), ') ',
                       species, ' (', scales::percent(species_confidence, scale = 1, accuracy = 1, suffix = ''), ')'),
         .keep = 'unused') 

the_targets %>%
  group_by(study, general_species) %>%
  summarise(n = n(),
            .groups = 'drop')


#
aligned_subset <- msa(all_sequences[the_targets$asv_id], verbose = TRUE)

msa::msaConvert(aligned_subset, type = 'ape::DNAbin') %>%
  ape::dist.dna(model = 'N') %>%
  as.matrix %>%
  as_tibble(rownames = 'asv_id1') %>%
  pivot_longer(cols = -asv_id1,
               names_to = 'asv_id2',
               values_to = 'n_differences') %>%
  mutate(study1 = if_else(str_detect(asv_id1, 'ASV'), 'ours', 'rosales'),
         study2 = if_else(str_detect(asv_id2, 'ASV'), 'ours', 'rosales')) %>%
  left_join(select(the_targets, asv_id, species) %>%
              dplyr::rename(species1 = species),
            by = c('asv_id1' = 'asv_id')) %>%
  left_join(select(the_targets, asv_id, species) %>%
              dplyr::rename(species2 = species),
            by = c('asv_id2' = 'asv_id')) %>% 
  
  # filter(asv_id1 == 'ASV25' | asv_id2 == 'ASV25',
  #        study1 != study2)
  
  filter(study1 != study2,
         n_differences == 0)




library(bioseq)
dist(aligned_subset)
aligned_subset %>%
  as.character() %>%
  as_dna() %>%
  
  
  enframe() %>%
  mutate(name = if_else(str_detect(name, 'PA_'), 'sequence', 'sequence_compare')) %>%
  pivot_wider()