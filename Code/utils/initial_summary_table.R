library(tidyverse)


#### Data ####
field_data <- read_csv('../../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  select(sample_id, health, year, season, site, lib.size) %>%
  distinct 

whacky <- c('Bin1_N_D_YE', 'Bin2_N_D_YE', 'Bin5_N_D_OR', 'Bin5_N_D_GR')
tank_data <- read_csv('../../intermediate_files/normalized_tank_asv_counts.csv', 
                          show_col_types = FALSE) %>% 
  select(sample_id, health, time, tank, anti, geno, plate, 
         year, season, site, dataset, exposure, lib.size) %>%
  distinct %>%
  # filter(asv_id == 'ASV25') %>%
  filter(anti == 'N',
         sample_id != 'P6_Bin2_N_N_GR') %>% 
  mutate(time = str_extract(time, '[0-9]') %>% as.integer(),
         fragment_id = str_c(geno, exposure, str_extract(tank, 'Bin[0-9]+'), sep = '_'),
         tank_id = str_c(exposure, str_extract(tank, 'Bin[0-9]+'), sep = '_'),
         treatment = if_else(exposure == 'pre', 'pre',
                             str_c(exposure, health, sep = '_')),
         time_treat = str_c(time, exposure, health, sep = '_')) %>%
  select(sample_id, health, time, tank, geno, exposure, 
         fragment_id, tank_id, treatment, lib.size) %>%
  
  #Get rid of samples which shouldnt exist
  filter(!(str_detect(sample_id, str_c(whacky, collapse = '|')) & time == 8))


alpha_metrics <- read_csv('../../intermediate_files/alpha_metrics.csv', show_col_types = FALSE)

#### Model richness & reads ####
library(MKinfer)
library(magrittr)
field_groups <- c('year', 'season')

test_differences_boot <- bind_rows('Field' = full_join(left_join(field_data, alpha_metrics, 
                                                            by = 'sample_id') %>%
                                                    group_by(across(all_of(field_groups))) %>%
                                                    filter(sum(health == 'D') > 1 & sum(health == 'H') > 1) %>%
                                                    summarise('reads.p' = boot.t.test(lib.size ~ health)$boot.p.value,
                                                              .groups = 'drop'),
                                                  
                                                  left_join(field_data, alpha_metrics, 
                                                            by = 'sample_id') %>%
                                                    group_by(across(all_of(field_groups))) %>%
                                                    filter(sum(health == 'D') > 1 & sum(health == 'H') > 1) %>%
                                                    summarise('richness.p' = boot.t.test(observed ~ health)$boot.p.value,
                                                              .groups = 'drop')) %>%
                                mutate(time = str_c(year, season, sep = ' - ') %>%
                                         str_replace_all(c('S$' = 'Jul', 'W$' = 'Jan')),
                                       exposure = '-',
                                       .keep = 'unused'),
                              
                              
                              'Tank' = bind_cols(left_join(tank_data, alpha_metrics, 
                                                           by = 'sample_id') %>%
                                                   filter(treatment %in% c('D_D', 'D_H')) %>%
                                                   summarise('reads.p' = boot.t.test(lib.size ~ health)$boot.p.value),
                                                 
                                                 left_join(tank_data, alpha_metrics, 
                                                           by = 'sample_id') %>%
                                                   filter(treatment %in% c('D_D', 'D_H')) %>%
                                                   summarise('richness.p' = boot.t.test(observed ~ health)$boot.p.value)) %>%
                                mutate(time = 'PostExposure', exposure = 'Disease'), 
                              .id = 'Source')


test_differences <- bind_rows('Field' = full_join(left_join(field_data, alpha_metrics, 
                                                            by = 'sample_id') %>%
                                                    group_by(across(all_of(field_groups))) %>%
                                                    filter(sum(health == 'D') > 1 & sum(health == 'H') > 1) %>%
                                                    summarise(broom::tidy(t.test(lib.size ~ health)),
                                                              .groups = 'drop') %>%
                                                    select(any_of(c('year', 'season', 'site')), statistic, parameter, p.value) %>%
                                                    rename(reads.t = statistic,
                                                           reads.df = parameter,
                                                           'reads.p' = p.value),
                                                  
                                                  left_join(field_data, alpha_metrics, 
                                                            by = 'sample_id') %>%
                                                    group_by(across(all_of(field_groups))) %>%
                                                    filter(sum(health == 'D') > 1 & sum(health == 'H') > 1) %>%
                                                    summarise(broom::tidy(t.test(observed ~ health)),
                                                              .groups = 'drop') %>%
                                                    select(any_of(c('year', 'season', 'site')), statistic, parameter, p.value) %>%
                                                    rename(richness.t = statistic,
                                                           richness.df = parameter,
                                                           'richness.p' = p.value)) %>%
                                mutate(time = str_c(year, season, sep = ' - ') %>%
                                         str_replace_all(c('S$' = 'Jul', 'W$' = 'Jan')),
                                       exposure = '-',
                                       .keep = 'unused'),
                              
                              
                              'Tank' = bind_cols(left_join(tank_data, alpha_metrics, 
                                                           by = 'sample_id') %>%
                                                   filter(treatment %in% c('D_D', 'D_H')) %>%
                                                   summarise(broom::tidy(t.test(lib.size ~ treatment))) %>%
                                                   select(statistic, parameter, p.value) %>%
                                                   rename(reads.t = statistic,
                                                          reads.df = parameter,
                                                          'reads.p' = p.value),
                                                 
                                                 left_join(tank_data, alpha_metrics, 
                                                           by = 'sample_id') %>%
                                                   filter(treatment %in% c('D_D', 'D_H')) %>%
                                                   summarise(broom::tidy(t.test(observed ~ treatment)))%>%
                                                   select(statistic, parameter, p.value) %>%
                                                   rename(richness.t = statistic,
                                                          richness.df = parameter,
                                                          'richness.p' = p.value)) %>%
                                mutate(time = 'PostExposure', exposure = 'Disease'), 
                              .id = 'Source')

#### Join together ####
field_prep <- field_data %>%
  left_join(alpha_metrics, 
            by = 'sample_id') %>%
  # filter(site == 'CK14', year == '2016', season == 'S') %>% select(sample_id:observed) %>% View
  group_by(across(c(all_of(field_groups), health))) %>%
  summarise(#n_healthy = sum(health == 'H'),
            #n_disease = sum(health == 'D'),
            n_fragments = n(),
            mean_reads = mean(lib.size),
            se_reads = sd(lib.size) / sqrt(n_fragments),
            mean_richness = mean(observed),
            se_richness = sd(observed) / sqrt(n_fragments),
            .groups = 'drop') %>%
  mutate(reads = case_when(is.na(se_reads) ~ scales::comma(mean_reads, accuracy = 1),
                           TRUE ~ str_c(scales::comma(mean_reads, accuracy = 1), scales::comma(se_reads, accuracy = 1), sep = ' ± ')),
         
         richness = case_when(is.na(se_richness) ~ scales::comma(mean_richness, accuracy = 1),
                              TRUE ~ str_c(scales::comma(mean_richness, accuracy = 1), scales::comma(se_richness, accuracy = 1), sep = ' ± ')),
         .keep = 'unused') %>%
  
  pivot_wider(names_from = 'health',
              values_from = c(n_fragments, reads, richness), 
              names_glue = "{health}_{.value}") %>%
  rename(n_healthy = D_n_fragments,
         n_disease = H_n_fragments) %>%
  rename_with(~str_replace_all(., c('D' = 'disease', 'H' = 'healthy')))

tank_prep <- tank_data %>%
  left_join(alpha_metrics, 
            by = 'sample_id') %>%
  # group_by(fragment_id) %>%
  # mutate(health = ifelse(any(health == 'D'), 'D', health)) %>%
  mutate(time = if_else(time == 0, 'pre', 'post')) %>%
  group_by(time, treatment) %>%
  summarise(n_genotypes = n_distinct(geno),
            n_healthy = sum(health == 'H'),
            n_disease = sum(health == 'D'),
            n_fragments = n_distinct(fragment_id),
            mean_reads = mean(lib.size),
            se_reads = sd(lib.size) / sqrt(n_fragments),
            mean_richness = mean(observed),
            se_richness = sd(observed) / sqrt(n_fragments),
            .groups = 'drop') %>%
  select(-n_fragments) %>%
  
  mutate(reads = case_when(is.na(se_reads) ~ scales::comma(mean_reads, accuracy = 1),
                           TRUE ~ str_c(scales::comma(mean_reads, accuracy = 1), scales::comma(se_reads, accuracy = 1), sep = ' ± ')),
         
         richness = case_when(is.na(se_richness) ~ scales::comma(mean_richness, accuracy = 1),
                              TRUE ~ str_c(scales::comma(mean_richness, accuracy = 1), scales::comma(se_richness, accuracy = 1), sep = ' ± ')),
         .keep = 'unused') %>%
  mutate(exposure = str_remove(treatment, '_[DH]') %>%
           str_replace_all(c('D' = 'Disease', 'N' = 'Control', 'pre' = '-')),
         .keep = 'unused', .after = 'time') %>%
  mutate(time = str_to_sentence(time) %>% str_c('Exposure'),
         time = fct_relevel(time, 'PreExposure')) %>%
  arrange(time, exposure)

tank_prep

tank_data %>%
  filter(treatment != 'pre',
         treatment != 'N_H') %>%
  count(treatment, geno) %>%
  pivot_wider(names_from = geno,
              values_from = n)

#### Merge ####
output_table <- bind_rows('Field' = field_prep,
          'Tank' = tank_prep,
          .id = 'Source') %>%
  mutate(time = coalesce(time, str_c(year, season, sep = ' - ')) %>%
           str_replace_all(c('S$' = 'Jul', 'W$' = 'Jan')),
         .keep = 'unused') %>%
  mutate(disease_reads = case_when(Source == 'Tank' & n_disease > 0 ~ reads,
                                   Source == 'Tank' & n_disease == 0 ~ '-',
                                   TRUE ~ disease_reads),
         healthy_reads = case_when(Source == 'Tank' & n_healthy > 0 ~ reads,
                                   Source == 'Tank' & n_healthy == 0 ~ '-',
                                   TRUE ~ healthy_reads),
         
         disease_richness = case_when(Source == 'Tank' & n_disease > 0 ~ richness,
                                   Source == 'Tank' & n_disease == 0 ~ '-',
                                   TRUE ~ disease_richness),
         healthy_richness = case_when(Source == 'Tank' & n_healthy > 0 ~ richness,
                                   Source == 'Tank' & n_healthy == 0 ~ '-',
                                   TRUE ~ healthy_richness),
         exposure = replace_na(exposure, '-')) %>%
  rowwise %>%
  mutate(n_genotypes = replace_na(n_genotypes, sum(n_healthy + n_disease))) %>%
  ungroup %>%
  select(-reads, -richness) %>%
  relocate(time, exposure, any_of(field_groups), n_genotypes, .after = Source) %>%
  group_by(across(c(Source, time, exposure, any_of(field_groups)))) %>%
  summarise(n_genotypes = max(n_genotypes),
            across(c(where(is.numeric), -n_genotypes), sum),
            across(where(is.character), ~ifelse(length(.) > 1, .[.!='-'], .)),
            .groups = 'drop') %>%
  arrange(Source, exposure, time) %>%
  left_join(select(test_differences_boot, time, exposure, any_of(field_groups), ends_with('p'))) %>%
  mutate(reads.p = scales::pvalue(reads.p, accuracy = 1/10000),
         reads.p = replace_na(reads.p, '-'),
         richness.p = scales::pvalue(richness.p, accuracy = 1/10000),
         richness.p = replace_na(richness.p, '-')) %>%
  relocate(reads.p, .after = healthy_reads) 
write_csv(output_table, '../../Results/Table 1.csv')


output_table_long <- output_table %>%
  rename(healthy_n = n_healthy,
         disease_n = n_disease) %>%
  pivot_longer(cols = c(starts_with('healthy'), starts_with('disease')),
               names_to = c('state', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  mutate(state = str_to_sentence(state)) %>%
  relocate(reads.p, .after = reads) %>%
  relocate(richness.p, .after = richness)

write_csv(output_table_long, '../../Results/Table 1 long.csv')
