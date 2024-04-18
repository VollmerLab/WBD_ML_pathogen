#### Libraries ####
library(Biostrings)
library(phyloseq)
library(tidyverse)
library(lmerTest)
library(emmeans)

#### Functions ####
top_classification <- function(data, threshold = 80){
  rename_with(data, .cols = domain:species,
              ~str_c(., '_name')) %>%
    pivot_longer(cols = -asv_id,
                 names_to = c('taxon_level', '.value'),
                 names_pattern = '(.*)_(.*)') %>%
    mutate(taxon_level = factor(taxon_level, ordered = TRUE,
                                levels = c('domain', 'phylum', 'class', 'order', 
                                           'family', 'genus', 'species'))) %>%
    filter(confidence >= threshold) %>%
    group_by(asv_id) %>%
    filter(taxon_level == max(taxon_level)) %>%
    ungroup %>%
    left_join(select(data, -ends_with('confidence')), by = 'asv_id')
  
}

#### Data ####
phyloseq_object <- read_rds('../../intermediate_files/prepped_microbiome.rds.gz') 



taxonomy <- read_csv('../../intermediate_files/update_taxonomy.csv',
                     show_col_types = FALSE) %>%
  mutate(across(where(is.character), str_replace_na)) %>%
  top_classification() %>%
  left_join(read_csv('../../intermediate_files/update_taxonomy.csv',
                     show_col_types = FALSE) %>%
              select(asv_id, ends_with('confidence')),
            by = 'asv_id') %>%
  bind_rows(read_csv('../../intermediate_files/taxonomy.csv.gz', show_col_types = FALSE) %>%
              filter(asv_id == 'ASV131') %>%
              mutate(name = genus,
                     taxon_level = 'genus')) %>%
  mutate(across(contains('confidence'), ~scales::percent(., scale = 1, accuracy = 0.1))) %>%
  mutate(top_classification = str_c(name, ' (', str_replace_na(confidence), ')'),
         .keep = 'unused') %>%
  select(-taxon_level) %>%
  rename_with(.cols = domain:species, ~str_c(., '_name')) %>%
  pivot_longer(cols = c(ends_with('name'), ends_with('confidence')),
               names_pattern = '(.*)_(.*)',
               names_to = c('taxon_level', '.value')) %>%
  
  mutate(name = str_c(taxon_level, '__', name, ' (', str_replace_na(confidence), ')'),
         .keep = 'unused') %>%
  group_by(asv_id, top_classification) %>%
  summarise(taxonomy = str_c(str_replace_na(name), collapse = '; '),
            .groups = 'drop') %>%
  bind_rows(tax_table(phyloseq_object) %>%
              as.data.frame() %>%
              as_tibble(rownames = 'asv_id') %>%
              mutate(across(where(is.character), str_replace_na)) %>%
              filter(!asv_id %in% taxonomy$asv_id) %>%
              pivot_longer(cols = -asv_id) %>%
              mutate(name = str_c(name, value, sep = '__')) %>%
              group_by(asv_id) %>%
              summarise(taxonomy = str_c(name, collapse = '; ')))


ml_results <- read_csv('../../Results/asv_importance.csv.gz',
         show_col_types = FALSE) %>%
  select(-df) %>%
  rename_with(~str_replace(., 'base', 'rank')) %>%
  rename(rank_median = median_rank,
         rank_model = response,
         rank_model_se = SE,
         rank_model_lower95 = asymp.LCL,
         rank_model_upper95 = asymp.UCL,
         rank_fdr = p_adjust,
         rank_alternative = alternative)

ml_kept_asvs <- ml_results %>%
  filter(rank_fdr < 0.05) %>%
  pull(asv_id)

field_results <- read_rds('../../intermediate_files/field_asv_models.rds.gz') %>%
  mutate(ml_retained = asv_id %in% ml_kept_asvs) %>%
  
  rowwise(asv_id, ml_retained) %>%
  reframe(fieldModel_healthCoefficient = fixef(model) %>%
            magrittr::extract('healthH') %>%
            unname,
          
          as_tibble(summary_table, rownames = 'term') %>%
            mutate(term = str_replace_all(term, ':', 'X')) %>%
            rename(SS = `Sum Sq`,
                   MS = `Mean Sq`,
                   nDF = NumDF,
                   dDF = DenDF,
                   fStat = `F value`,
                   pValue = `Pr(>F)`) %>%
            pivot_wider(names_from = term,
                        values_from = c(everything(), -term),
                        names_vary = 'slowest',
                        names_glue = 'fieldModel_{.value}_{term}'),
          
          contrast(all_health_fit, 'pairwise') %>% 
            broom::tidy(conf.int = TRUE) %>%
            as_tibble %>%
            select(-term, -contrast, -null.value) %>%
            mutate(timepoint = str_c(year, season, sep = ''),
                   .keep = 'unused') %>%
            rename(mean = estimate,
                   se = std.error,
                   lower95 = conf.low,
                   upper95 = conf.high,
                   tStat = statistic,
                   pValue = p.value) %>%
            pivot_wider(names_from = timepoint, 
                        values_from = c(everything(), -timepoint),
                        names_vary = 'slowest',
                        names_glue = 'fieldLogFC_{.value}_{timepoint}')) %>%
  mutate(across(contains('fieldLogFC'), 
                ~if_else(!ml_retained, NA_real_, .)),
         
         across(contains('pValue'), 
                ~p.adjust(., 'fdr', n = sum(!is.na(.))), 
                .names = '{.col}_fdr')) %>%
  rename_with(~str_replace(., 'pValue', 'fdr') %>% 
                str_remove('_fdr$'), 
              .cols = ends_with('_fdr'))

field_kept_asvs <- field_results %>%
  filter(ml_retained) %>%
  select(asv_id, contains('fieldLogFC') & (contains('fdr') | contains('meanDH'))) %>%
  pivot_longer(cols = -asv_id,
               names_prefix = 'fieldLogFC_',
               names_pattern = '(.*)_(.*)',
               names_to = c('.value', 'season')) %>%
  group_by(asv_id) %>%
  summarise(field_retained = all(fdr < 0.05) & (all(meanDH > 0) | all(meanDH < 0))) %>%
  filter(field_retained) %>%
  pull(asv_id)

tank_results <- read_rds('../../intermediate_files/tank_asv_models.rds.gz') %>%
  filter(asv_id %in% field_kept_asvs) %>%
  select(asv_id, summary_table, tank_posthoc) %>%
  rowwise(asv_id) %>%
  reframe(as_tibble(summary_table) %>%
            rename(SS = `Sum Sq`,
                   MS = `Mean Sq`,
                   nDF = NumDF,
                   dDF = DenDF,
                   fStat = `F value`,
                   pValue = `Pr(>F)`) %>%
            rename_with(~str_c('tankModel', ., 
                               'timeXexposureXoutcome', 
                               sep = '_')),
          
          rename(tank_posthoc,
                 mean = estimate,
                 se = std.error,
                 lower95 = conf.low,
                 upper95 = conf.high,
                 tStat = statistic,
                 pValue = p.value) %>%
            pivot_wider(names_from = contrast, 
                        values_from = c(everything(), -contrast),
                        names_vary = 'slowest',
                        names_glue = 'tankLogFC_{.value}_{contrast}')) %>%
  
  mutate(across(contains('pValue'), 
                ~p.adjust(., 'fdr', n = sum(!is.na(.))), 
                .names = '{.col}_fdr')) %>%
  rename_with(~str_replace(., 'pValue', 'fdr') %>% 
                str_remove('_fdr$'), 
              .cols = ends_with('_fdr'))

field_data <- read_csv('../../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  select(asv_id:norm.factors)


whacky <- c('Bin1_N_D_YE', 'Bin2_N_D_YE', 'Bin5_N_D_OR', 'Bin5_N_D_GR')
tank_data <-  read_csv('../../intermediate_files/normalized_tank_asv_counts.csv', 
                       show_col_types = FALSE) %>% 
  # filter(asv_id == 'ASV25') %>%
  filter(anti == 'N',
         sample_id != 'P6_Bin2_N_N_GR') %>% 
  mutate(time = str_extract(time, '[0-9]') %>% as.integer(),
         fragment_id = str_c(geno, exposure, str_extract(tank, 'Bin[0-9]+'), sep = '_'),
         tank_id = str_c(exposure, str_extract(tank, 'Bin[0-9]+'), sep = '_'),
         treatment = if_else(exposure == 'pre', 'pre',
                             str_c(exposure, health, sep = '_')),
         time_treat = str_c(time, exposure, health, sep = '_')) %>%
  select(-dataset, -anti, -plate, -year:-site, -cpm_norm:-cpm, 
         -domain:-species) %>%
  
  #Get rid of samples which shouldnt exist
  filter(!(str_detect(sample_id, str_c(whacky, collapse = '|')) & time == 8))



#### Supplementary Files ####
#Field Metadata
field_metadata <- select(field_data, -asv_id, -log2_cpm_norm:-n_reads, -dataset) %>%
  distinct
write_csv(field_metadata, '../../Results/SupplementaryFile1.csv.gz')

#Tank Metadata
tank_metadata <- select(tank_data, -asv_id, -log2_cpm_norm:-n_reads) %>%
  distinct %>%
  select(sample_id, geno, fragment_id, tank_id, time_treat, lib.size, norm.factors)
write_csv(tank_metadata, '../../Results/SupplementaryFile2.csv.gz')

# ASV, sequence, taxonomy - fasta
the_sequence <- refseq(phyloseq_object)
names(the_sequence) <- names(the_sequence) %>%
  tibble(asv_id = .) %>%
  left_join(select(taxonomy, asv_id, taxonomy),
            by = 'asv_id') %>%
  mutate(asv_id = str_c(asv_id, taxonomy, sep = '; ')) %>%
  pull(asv_id)
writeXStringSet(the_sequence, filepath = '../../Results/SupplementaryFile3.fasta')

#Field & Tank ASV Count matrix
asv_counts <- bind_rows(select(field_data, asv_id, sample_id, n_reads),
                        select(tank_data, asv_id, sample_id, n_reads)) %>%
  pivot_wider(names_from = 'asv_id',
              values_from = 'n_reads')
write_csv(asv_counts, '../../Results/SupplementaryFile4.csv.gz')

#Field & Tank ASV log2CPM matrix
asv_log2cpm <- bind_rows(select(field_data, asv_id, sample_id, log2_cpm_norm),
                         select(tank_data, asv_id, sample_id, log2_cpm_norm)) %>%
  pivot_wider(names_from = 'asv_id',
              values_from = 'log2_cpm_norm')
write_csv(asv_log2cpm, '../../Results/SupplementaryFile5.csv.gz')

# ASV, taxonomy, Full ML rankings, field results, tank results - csv

complete_asv_results <- full_join(taxonomy, 
                                  ml_results,
                                  by = 'asv_id') %>%
  full_join(field_results,
            by = 'asv_id') %>%
  full_join(tank_results,
            by = 'asv_id') %>%
  mutate(passedFilter = !is.na(ml_retained), 
         .after = taxonomy) %>%
  mutate(diffAbundance_healthAssociation = case_when(fieldModel_fdr_health < 0.05 & fieldModel_healthCoefficient > 0 ~ 'H',
                                                     fieldModel_fdr_health < 0.05 & fieldModel_healthCoefficient < 0 ~ 'D',
                                                     is.na(fieldModel_healthCoefficient) ~ NA_character_,
                                                     TRUE ~ 'N'),
         .after = fieldModel_fdr_healthXyear) %>%
  select(-fieldModel_healthCoefficient) %>%
  mutate(field_consistent = case_when(!is.na(fieldLogFC_mean_2016S) ~ asv_id %in% field_kept_asvs,
                                      TRUE ~ NA),
         .after = fieldLogFC_fdr_2017W) %>%
  mutate(pathogen_classification = case_when(is.na(field_consistent) | !field_consistent ~ NA,
                                             asv_id %in% str_c('ASV', c(25, 8, 38)) ~ 'pathogen', 
                                             asv_id %in% str_c('ASV', c(26, 30, 361, 51)) ~ 'opportunist',
                                             TRUE ~ 'commensal'))


write_csv(complete_asv_results, '../../Results/SupplementaryFile6.csv.gz')


write_rds(phyloseq_object, '../../Results/SupplementaryFile7.rds.gz', compress = 'gz')
