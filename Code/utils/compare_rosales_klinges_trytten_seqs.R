library(Biostrings)
library(msa)
library(bioseq)
library(phyloseq)
library(tidyverse)
library(magrittr)

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

#### Prep ####
if(!file.exists('../../intermediate_files/klinges_2022.fasta')){
  read_tsv('../../intermediate_files/klinges_taxonomy_unrare_fullseq.txt', 
           col_names = c('sequence', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
           show_col_types = FALSE, skip = 1) %>%
    mutate(asv_id = str_c('klinges_ASV', row_number())) %$%
    set_names(sequence, asv_id) %>%
    DNAStringSet() %>%
    writeXStringSet('../../intermediate_files/klinges_2022.fasta')
}

if(!file.exists('../../intermediate_files/trytten_unpub.fasta')){
  florida_sequences <- read_rds('../../../Emily_Microbiome/intermediate_files/preprocess_microbiome.rds') %>%
    refseq
  names(florida_sequences) <- str_c('FL_ASV', 1:length(names(florida_sequences)), sep = '_')
  writeXStringSet(florida_sequences, '../../intermediate_files/trytten_unpub.fasta')
}


#### Data ####
top_asvs <- read_csv('../../Results/Table45_asv_table.csv', 
                     show_col_types = FALSE) %>%
  filter(likely_type != '-') %>%
  select(ID, likely_type) %>%
  filter(likely_type == 'Pathogen') %>%
  select(ID) %>%
  dplyr::rename(asv_id = ID) %>%
  mutate(asv_id = str_c('PA_', asv_id))

#### ASV Differential Expression ####
#https://github.com/srosales712/coral_disease_resistance_microbiome/blob/master/Figures_DR_MB1.ipynb
#each "chunk" below from different block - stated after block
rosales_asvs <- c("e2829560ad1e796bebb56d909c7fb12d",
                  "5c518d771ed47c631781c4e204c796b0",
                  "fa87a3ed888b44bab9d77e04cbb99665",
                  ##select ASVs that were significant for outcome from ANCOM results - palmata
                  
                  "7eb68c2ff12bb8a0a46d036c37f8f26e",
                  "fa87a3ed888b44bab9d77e04cbb99665",
                  ##select ASVs that were significant for outcome from ANCOM results - cervcornis
                  
                  "3f8018e35043db82841f695939b95ece",
                  "87624cf95e61dd0ce58bab24686f1eb4",
                  "c6af4b65b2f25b905461f80f13b89edd",
                  "7ebeae7195a9f372afdc9d94744cc349",
                  "5b79cf6d5a5a9bf0bb866aed449eff44",
                  ##select ASVs that were significant for genotype from ANCOM results - cervicornis
                  
                  "8415141e03996c7309d7d1fd8aee6967",
                  "7253c0f3989180f5b68f83c04ecfbda2",
                  "7ebeae7195a9f372afdc9d94744cc349"
                  ##select ASVs that were significant for genotype from ANCOM results - palmata
                  ) %>%
  unique


klinges_asvs <- str_c('klinges_ASV',
                      c(1, 3, 4, 67,
                        #From supplementary table S7 - differential abundance or variance by nutrient treatment
                        
                        2, 6, 3
                        #Stated in results as being "dominant"
                        )) %>%
  unique


trytten_asvs <- read_csv('../../../Emily_Microbiome/intermediate_files/classified_significant_asvs.csv.gz',
                         show_col_types = FALSE) %>%
  filter(DiseaseOutcome) %>%
  pull(asv_id) %>%
  str_c('FL_', .)

sweet_asvs <- str_c('sweet_',
                    c('KC737024', 'KC737026', 'KC737032'))

#### Taxonomy Data ####
panama_taxonomy <- read_blca('../../intermediate_files/all_asvs.fasta.blca.out') %>%
  mutate(asv_id = str_c('PA_', asv_id))

rosales_taxonomy <- read_blca('../../intermediate_files/rosales_2019.fasta.blca.out')

klinges_taxonomy <- read_blca('../../intermediate_files/klinges_2022.fasta.blca.out')

trytten_taxonomy <- read_blca('../../intermediate_files/trytten_unpub.fasta.blca.out')

sweet_taxonomy <- read_blca('../../intermediate_files/sweet_2014.fasta.blca.out') %>%
  mutate(asv_id = str_c('sweet_', asv_id) %>% str_remove('\\.1$'))

#### Sequence Data ####
panama_sequences <- read_rds('../../intermediate_files/prepped_microbiome.rds.gz') %>%
  refseq
names(panama_sequences) <- str_c('PA_', names(panama_sequences))

rosales_sequences <- readDNAStringSet('../../intermediate_files/rosales_2019.fasta')

#https://github.com/graceklinges/disease-susceptible-micro/blob/main/outputs/taxonomy.txt
klinges_sequences <- readDNAStringSet('../../intermediate_files/klinges_2022.fasta')

trytten_sequences <- readDNAStringSet('../../intermediate_files/trytten_unpub.fasta')

sweet_sequences <- readDNAStringSet('../../intermediate_files/sweet_2014.fasta')
names(sweet_sequences) <- str_c('sweet_', str_extract(names(sweet_sequences), 'K[A-Za-z0-9]+'))

#### ASV Matching ####
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
    ungroup
}

data <- matched_sequences$compare_taxa[[4]]; level = matched_sequences$taxon_level[[4]]; id = matched_sequences$name[[4]]
filter_taxa <- function(data, level, id){
  filter(data, !!sym(level) == id)
}

filter_taxa(matched_sequences$compare_taxa[[12]], level = matched_sequences$taxon_level[[12]], id = matched_sequences$name[[12]]) %>%
  top_classification()

matched_sequences <- inner_join(panama_taxonomy, 
           top_asvs,
           by = 'asv_id') %>%
  top_classification %>%
  mutate(taxon_level = as.character(taxon_level)) %>%
  select(-confidence) %>%
  mutate(sequence = as_dna(as.character(panama_sequences[asv_id]))) %>%
  expand_grid(datasouce = c('rosales', 'klinges', 'trytten', 'sweet')) %>%
  mutate(compare_taxa = case_when(datasouce == 'rosales' ~ list(rosales_taxonomy),
                                  datasouce == 'klinges' ~ list(klinges_taxonomy),
                                  datasouce == 'trytten' ~ list(trytten_taxonomy),
                                  datasouce == 'sweet' ~ list(sweet_taxonomy),
                                  TRUE ~ list(NULL)),
         compare_seqs = case_when(datasouce == 'rosales' ~ list(rosales_sequences),
                                  datasouce == 'klinges' ~ list(klinges_sequences),
                                  datasouce == 'trytten' ~ list(trytten_sequences),
                                  datasouce == 'sweet' ~ list(sweet_sequences),
                                  TRUE ~ list(NULL))) %>%
  filter(lengths(compare_taxa) > 0) %>%
  rowwise() %>%
  mutate(compare_taxa = list(filter_taxa(compare_taxa, level = taxon_level, id = name) %>%
                               top_classification())) %>%
  filter(nrow(compare_taxa) > 0) %>%
  mutate(compare_seqs = list(compare_seqs[compare_taxa$asv_id])) %>%
  
  
  mutate(compare = list(as_dna(as.character(compare_seqs)) %>%
           tibble(asv_id = names(.),
                  sequence = .) %>%
           left_join(compare_taxa,
                     by = 'asv_id') %>%
             rename_with(~str_c(., '_compare'))), 
         .keep = 'unused') %>%
  unnest(compare)

#### Align Sequences ####
aligned_sequences <- matched_sequences %>%
  mutate(seqLength = str_length(sequence),
         seqLength_compare = str_length(sequence_compare)) %>%
  rowwise %>%
  mutate(msa(c(DNAStringSet(sequence), DNAStringSet(sequence_compare)), verbose = FALSE) %>%
           as.character() %>%
           as_dna() %>%
           enframe() %>%
           mutate(name = if_else(str_detect(name, 'PA_'), 'sequence', 'sequence_compare')) %>%
           pivot_wider(),
         .keep = 'unused')

#### Trim to Core ####
# seq1 <- aligned_sequences$sequence[[1]]
# seq2 <- aligned_sequences$sequence_compare[[1]]
trim_to_core <- function(seq1, seq2){
  leading_missing <- c(seq1, seq2) %>%
    str_extract('^-*') %>% 
    str_length()
  
  tail_missing <- c(seq1, seq2) %>%
    str_extract('-*$') %>% 
    str_length()
  
  the_core <- c(seq1, seq2) %>%
    seq_remove_position(1, max(leading_missing)) %>%
    seq_remove_position(str_length(.) - max(tail_missing), str_length(.)) %>%
    enframe() %>%
    mutate(name = if_else(str_detect(name, '1'), 'sequence', 'sequence_compare')) %>%
    pivot_wider()
   
  the_core
}

trimmed_sequences <- mutate(aligned_sequences, 
                            trim_to_core(sequence, sequence_compare),
                            .keep = 'unused')

#### Count Sequence Differences ####
make_base_position <- function(the_seq){
  as.character(the_seq) %>%
    str_split('') %>%
    unlist %>%
    tibble(bp = .) %>%
    mutate(position = row_number()) %>%
    select(position, bp)
}

mismatched_sequences <- trimmed_sequences %>%
  mutate(across(c(sequence, sequence_compare), ~list(make_base_position(.))),
         sequence_compare = list(rename_with(sequence_compare, ~str_c(., '_compare')))) %>%
  mutate(both_seqs = list(full_join(sequence, sequence_compare, by = c('position' = 'position_compare'))),
         .keep = 'unused') %>%
  unnest(both_seqs) %>%
  group_by(across(asv_id:seqLength_compare)) %>%
  summarise(n_base = n(),
            n_original = sum(bp != '-'),
            n_compare = sum(bp_compare != '-'),
            n_match = sum(bp == bp_compare),
            n_mismatch = sum(bp != bp_compare),
            .groups = 'drop')
  

mismatched_sequences %>%
  filter((asv_id_compare %in% rosales_asvs & datasouce == 'rosales') |
           (asv_id_compare %in% klinges_asvs & datasouce == 'klinges') | #(datasouce == 'klinges') | #Klinges significance is associated with nutrients - all in "disease susceptible" corals
           (asv_id_compare %in% trytten_asvs & datasouce == 'trytten') |
           (asv_id_compare %in% sweet_asvs & datasouce == 'sweet')) %>%
  mutate(pct_mismatch = n_mismatch / n_compare) %>%
  select(asv_id, name, datasouce, asv_id_compare, name_compare, n_mismatch, n_compare, pct_mismatch) %>%
  select(-name_compare)


rosales_original_taxonomy <- read_tsv('../../Data/rosales/taxonomy.tsv', show_col_types = FALSE) %>%
  filter(`Feature ID` %in% c('5b79cf6d5a5a9bf0bb866aed449eff44', '7eb68c2ff12bb8a0a46d036c37f8f26e'))

#Rosales: 5b79cf6d5a5a9bf0bb866aed449eff44 i.e. our ASV25
# Sup. Fig. 3
# Present/abundant in genotypes: C20,21,22,24,28,29,30
