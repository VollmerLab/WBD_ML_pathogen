library(tidyverse)
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