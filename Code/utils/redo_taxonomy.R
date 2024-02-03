library(tidyverse)

#### 
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
  subset_samples(site != 'MG') %>% #no disease samples
  prune_samples(!str_detect(rownames(sample_data(.)), 
                            str_c(remove_samples, collapse = '|')), .)


## Reclassify with Decipher 
Biostrings::writeXStringSet(file="../intermediate_files/all_asvs.fasta", refseq(microbiome_data))

#### Redo with Decipher ####
#http://www2.decipher.codes/ClassifyOrganisms.html


load('~/../Downloads/decipherOut_50pct.RData')

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
new_taxonomy <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(new_taxonomy) <- ranks
as_tibble(new_taxonomy, rownames = 'asv_id') %>%
  write_csv('../intermediate_files/update_taxonomy.csv')

#### Redo with BCLA ####

new_taxonomy <- read_delim('../../intermediate_files/all_asvs.fasta.blca.out', delim = '\t', 
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

write_csv(new_taxonomy, '../../intermediate_files/update_taxonomy.csv')
