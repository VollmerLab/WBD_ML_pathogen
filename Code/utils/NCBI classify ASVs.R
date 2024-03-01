library(tidyverse)

alpha <- 0.05

taxonomy <- read_csv('../../intermediate_files/taxonomy.csv.gz',
                     show_col_types = FALSE) %>%
  mutate(across(where(is.character), 
                ~str_replace_na(., replacement = '')))

taxonomy_confidence <- read_csv('../../intermediate_files/update_taxonomy.csv',
                                show_col_types = FALSE)

shap_importance <- read_csv('../../Results/asv_importance.csv.gz',
                            show_col_types = FALSE) %>%
  filter(p_adjust < alpha) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(asv_id = fct_reorder(asv_id, median_rank)) %>%
  select(asv_id, domain:genus) %>%
  distinct %>%
  arrange(asv_id)

#### Get sequences of ASVs and BLAST ####
if(!file.exists('../../intermediate_files/interesting_asv_16s.fasta')){
  read_rds('../../intermediate_files/prepped_microbiome.rds.gz') %>%
    prune_taxa(taxa_names(.) %in% as.character(shap_importance$asv_id), .) %>%
    refseq %>%
    Biostrings::writeXStringSet(filepath = '../../intermediate_files/interesting_asv_16s.fasta')
}

## BLAST on NCBI ##

if(file.exists('../../intermediate_files/interesting_asv_blast_classifications.csv.gz')){
  the_taxa_complete <- read_csv('../../intermediate_files/interesting_asv_blast_classifications.csv.gz', show_col_types = FALSE)
} else {
  # library(multidplyr)
  library(rentrez)
  library(taxize)
  
  find_taxonomy <- function(sub_id){
    nuccore_search <- entrez_search(db = 'nuccore', term = sub_id)
    
    taxon_link <- entrez_link(dbfrom = 'nuccore', id = nuccore_search$ids, db = 'taxonomy')
    
    taxon_result <- entrez_fetch(db = 'taxonomy', id = taxon_link$links$nuccore_taxonomy, rettype = 'txt')
    taxon_result
  }
  
  get_higher_taxa <- function(species_name){
    initial_taxa <- classification(species_name, db = 'ncbi') 
    
    initial_taxa[[1]] %>%
      as_tibble() %>%
      filter(rank %in% c('domain', 'kingdom', 'phylum', 
                         'class', 'order', 'family', 
                         'genus', 'species')) %>%
      select(rank, name) %>%
      pivot_wider(names_from = rank, values_from = name)
  }
  
  # the_cluster <- new_cluster(parallel::detectCores() - 1)
  # cluster_library(the_cluster, c('dplyr', 'rentrez', 'taxize',
  #                                'tibble', 'tidyr', 'purrr'))
  # cluster_copy(the_cluster, c('find_taxonomy', 'get_higher_taxa'))
  
  the_taxa <- read_csv('../../intermediate_files/UP5UBP83013-Alignment-HitTable.csv',
                       col_names = c("asv_id", "subject.id", "identity", "alignment.length", 
                                     "mismatches", "gap.opens", "q.start", "q.end", "s.start", 
                                     "s.end", "evalue", "bit.score"),
                       show_col_types = FALSE) %>%
    nest_by(subject.id) %>%
    # partition(the_cluster) %>%
    mutate(ncbi_species_raw = possibly(find_taxonomy, otherwise = NA_character_)(subject.id)) %>%
    # collect %>%
    ungroup %>%
    mutate(ncbi_species = str_remove(ncbi_species_raw, '^.*Error.*\\\n') %>%
             str_extract(' .*\\\n') %>% 
             str_remove('^ ') %>% 
             str_remove('\\\n$'))
  
  successful_finds <- filter(the_taxa, !is.na(ncbi_species))
  failed_finds <- filter(the_taxa, is.na(ncbi_species))
  
  z <- 0
  if(nrow(failed_finds) > 0 & z <= 5){
    z <- z + 1
    new_finds <- failed_finds %>%
      select(-starts_with('ncbi')) %>%
      rowwise() %>%
      mutate(ncbi_species_raw = possibly(find_taxonomy, otherwise = NA_character_)(subject.id)) %>%
      # collect %>%
      ungroup %>%
      mutate(ncbi_species = str_remove(ncbi_species_raw, '^.*Error.*\\\n') %>%
               str_extract(' .*\\\n') %>% 
               str_remove('^ ') %>% 
               str_remove('\\\n$'))
    
    new_success <- filter(new_finds, !is.na(ncbi_species))
    successful_finds <- bind_rows(successful_finds, new_success)
    
    failed_finds <- filter(new_finds, is.na(ncbi_species))
  }
  the_taxa <- bind_rows(successful_finds, failed_finds)
  
  the_taxa_complete <- the_taxa %>%
    nest_by(ncbi_species) %>%
    # partition(the_cluster) %>%
    mutate(possibly(get_higher_taxa, 
                    otherwise = tibble('domain' = NA_character_, 'kingdom' = NA_character_, 'phylum' = NA_character_, 
                                       'class' = NA_character_, 'order' = NA_character_, 'family' = NA_character_, 
                                       'genus' = NA_character_, 'species' = NA_character_))(ncbi_species)) %>%
    # collect %>%
    ungroup %>%
    unnest(data, keep_empty = TRUE) %>%
    unnest(data, keep_empty = TRUE) %>%
    janitor::remove_empty(which = 'cols')
  
  
  write_csv(the_taxa_complete, '../../intermediate_files/interesting_asv_blast_classifications.csv.gz')
}


filter(taxonomy_confidence, asv_id %in% str_c('ASV', c(30, 26))) %>% 
  View

the_taxa_complete %>%
  filter(asv_id %in% c('ASV8', 'ASV38')) %>%
  group_by(asv_id) %>%
  filter(evalue == min(evalue)) %>%
  ungroup %>%
  filter(!str_detect(ncbi_species, 'sp\\.'),
         !str_detect(ncbi_species, 'bacterium')) %>%
  group_by(ncbi_species, asv_id) %>%
  filter(identity == max(identity),
         bit.score == max(bit.score)) %>%
  select(ncbi_species, asv_id) %>%
  distinct


the_taxa_complete %>%
  filter(asv_id %in% str_c('ASV', c(30, 26))) %>%
  group_by(asv_id) %>%
  filter(evalue == min(evalue)) %>%
  ungroup %>%
  filter(!str_detect(ncbi_species, 'sp\\.'),
         !str_detect(ncbi_species, 'bacterium')) %>%
  group_by(ncbi_species, asv_id) %>%
  filter(identity == max(identity),
         bit.score == max(bit.score)) %>%
  select(family, ncbi_species, asv_id) %>%
  distinct 


the_taxa_complete %>%
  filter(asv_id %in% str_c('ASV', c(361, 51))) %>%
  group_by(asv_id) %>%
  filter(evalue == min(evalue)) %>%
  ungroup %>%
  filter(!str_detect(ncbi_species, 'sp\\.'),
         !str_detect(ncbi_species, 'bacterium')) %>%
  group_by(ncbi_species, asv_id) %>%
  filter(identity == max(identity),
         bit.score == max(bit.score)) %>%
  select(family, ncbi_species, asv_id) %>%
  distinct 

the_taxa_complete %>%
  filter(asv_id %in% c('ASV40')) %>%
  group_by(asv_id) %>%
  filter(evalue == min(evalue)) %>%
  ungroup %>%
  filter(!str_detect(ncbi_species, 'sp\\.'),
         !str_detect(ncbi_species, 'bacterium')) %>%
  group_by(ncbi_species, asv_id) %>%
  filter(identity == max(identity),
         bit.score == max(bit.score)) %>%
  select(family, ncbi_species, asv_id) %>%
  distinct



#### GENERAL FILTERING ####
#https://chanzuckerberg.zendesk.com/hc/en-us/articles/360050963352-A-Guide-to-BLAST
#Most hits to a single species regardless of score
the_taxa_complete %>%
  filter(!is.na(genus)) %>%
  group_by(asv_id) %>%
  count(phylum, class, order, family, genus, species) %>%
  filter(n == max(n))


#Max bitscore then "vote" by most hits to a single species
the_taxa_complete %>%
  filter(!is.na(genus)) %>%
  group_by(asv_id) %>%
  filter(bit.score == max(bit.score)) %>%
  count(phylum, class, order, family, genus, species) %>%
  filter(n == max(n)) 


#Min e-value then max alignment length then "vote"
the_taxa_complete %>%
  filter(!is.na(genus)) %>%
  group_by(asv_id) %>%
  filter(evalue == min(evalue)) %>%
  filter(alignment.length == max(alignment.length)) %>%
  count(phylum, class, order, family, genus, species) %>%
  filter(n == max(n)) %>% View

the_taxa_complete %>%
  filter(asv_id == 'ASV8') %>% 
  arrange(evalue, alignment.length) %>% filter(species == 'Vibrio coralliilyticus')

#Chat GPT suggestion
holistic_metric_function <- function(e_value, percent_identity, bit_score, alignment_length, w = c(-3, 2, 1, 1)){
  w[1] * log(e_value + 1e-300) + w[2] * percent_identity + w[3] * bit_score + w[4] * alignment_length
}


the_taxa_complete %>%
  # filter(!is.na(genus)) %>%
  mutate(holistic = holistic_metric_function(evalue, identity, bit.score, alignment.length, c(-5, 3, 4, 2))) %>%
  group_by(asv_id) %>%
  filter(holistic == max(holistic)) %>%
  count(phylum, class, order, family, genus, species) %>%
  filter(n == max(n)) %>% View


recursive_vote <- function(data, taxonomic_levels) {
  result <- data
  for (level in taxonomic_levels) {
    if(level == 'species'){
      if(sum(str_detect(result$species, 'sp\\.')) < nrow(result)){
        result <- filter(result, !str_detect(species, 'sp\\.'))
      }
    }
    result <- result %>%
      count(!!sym(level)) %>%
      filter(n == max(n)) %>%
      select(all_of(level)) %>%
      left_join(result, by = level)
  }
  
  select(result, all_of(taxonomic_levels)) %>%
    distinct
}

the_taxa_complete %>%
  nest_by(asv_id) %>%
  reframe(recursive_vote(data, c('phylum', 'class', 'order', 'family', 'genus', 'species'))) 

the_taxa_complete %>%
  filter(!is.na(genus)) %>%
  mutate(holistic = holistic_metric_function(evalue, identity, bit.score, alignment.length, c(-5, 3, 4, 2))) %>%
  group_by(asv_id) %>%
  filter(holistic == max(holistic)) %>%
  ungroup %>%
  nest_by(asv_id) %>%
  reframe(recursive_vote(data, c('phylum', 'class', 'order', 'family', 'genus', 'species'))) %>% 
  mutate(asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>% View

filter(the_taxa_complete, asv_id == 'ASV25') %>%
  mutate(holistic = holistic_metric_function(evalue, identity, bit.score, alignment.length, c(-5, 3, 4, 2))) %>%
  arrange(-holistic) %>% View
