library(tidyverse)
library(magrittr)
library(broom)

alpha <- 0.05
asv_ranks <- read_csv('../../Results/asv_importance.csv.gz', show_col_types = FALSE)
taxonomy <- read_csv('../../intermediate_files/taxonomy.csv.gz', show_col_types = FALSE)


final_asvs <- read_csv('../../Results/Table45_asv_table.csv')


top_taxonomy <- asv_ranks %>%
  filter(p_adjust < 0.05) %>%
  select(asv_id) %>%
  left_join(taxonomy, by = 'asv_id')

count(top_taxonomy, family) %>% arrange(-n)
count(top_taxonomy, family, genus) %>% arrange(-n)



filter(final_asvs, Field != 'n/c') %>%
  select(ID) %>%
  filter(ID != 'ASV40') %>%
  left_join(taxonomy,
            by = c('ID' = 'asv_id')) %>%
  count(genus) %>% 
  arrange(-n)

filter(taxonomy, asv_id == 'ASV40')
read_csv('../../intermediate_files/update_taxonomy.csv', show_col_types = FALSE) %>%
  filter(asv_id == 'ASV40')

read_csv('../../Results/ASV_oraGSEA.csv.gz', show_col_types = FALSE) %>%
  filter(ora_p.value < 0.05) %>%
  pull(family)


#### Overrepresentation of Families ####
row_fisher <- function(x, k, m, N, direction = 'two.sided'){
  ##https://dputhier.github.io/ASG/practicals/go_statistics_td/go_statistics_td_2015.html
  data.frame(significant = c(x, m - x),
             not_significant = c(k - x, N - m - (k - x))) %>%
    fisher.test(alternative = direction) %>%
    tidy
}

ora_asv <- function(asv_ranks, min_size = 5, lowest_level = 'family'){
  asv_ranks %>%
    mutate(p_adjust = if_else(is.na(p_adjust), 1, p_adjust)) %>%
    select(asv_id, p_adjust) %>%
    left_join(taxonomy, by = 'asv_id') %>%
    group_by(across(domain:!!sym(lowest_level))) %>%
    summarise(n_sig_fam = sum(p_adjust < alpha),
              total_fam = n(),
              .groups = 'rowwise') %>%
    mutate(total_sig = sum(.$n_sig_fam),
           total_asv = sum(.$total_fam)) %>%
    filter(total_fam >= min_size) %>%
    summarise(row_fisher(n_sig_fam, total_fam, total_sig, total_asv, direction = 'greater'),
              .groups = 'drop') %>%
    rename_with(.cols = -c(domain:!!sym(lowest_level)),
                ~str_c('ora', ., sep = '_'))
}

overrepresented_families <- ora_asv(asv_ranks, min_size = 5, lowest_level = 'family')
filter(overrepresented_families, ora_p.value < 0.05)
