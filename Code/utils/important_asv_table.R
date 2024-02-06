
library(tidyverse)
library(emmeans)
alpha <- 0.05

taxonomy <- read_csv('../../intermediate_files/taxonomy.csv.gz', 
                     show_col_types = FALSE) %>%
  mutate(across(everything(), ~str_replace_na(., '')))

taxonom_confidence <- read_csv('../../intermediate_files/update_taxonomy.csv',
         show_col_types = FALSE)


ml_model_out <- read_csv('../../Results/asv_importance.csv.gz', show_col_types = FALSE) %>%
  left_join(taxonomy, ., by = 'asv_id') %>%
  rename(FDR = p_adjust) %>%
  # group_by(family, genus) %>%
  # filter(any(FDR < 0.05)) %>%
  # ungroup %>%
  filter(FDR < 0.05) %>%
  arrange(domain, phylum, class, order, family, genus, median_rank) %>%
  select(-contains('method'), -contains('alternative'), -contains('domain')) %>%
  relocate(asv_id, .after = 'genus') %>%
  rename_with(~str_remove(., 'base_') %>% 
                str_replace_all(c('forest' = 'Random Forest',
                                  'knn' = 'KNN',
                                  'lasso' = 'LASSO',
                                  'mlp' = 'MLP',
                                  'null' = 'Null',
                                  'pls' = 'PLS',
                                  'svmLinear' = 'SVM'))) %>%
  mutate(#p.value = scales::pvalue(p.value, 0.0001),
         FDR = scales::pvalue(FDR, 0.001)) %>%
  arrange(median_rank)



field_asv_models <- read_rds('../../intermediate_files/field_asv_models.rds.gz') %>% 
  filter(asv_id %in% ml_model_out$asv_id) %>%
  rowwise %>%
  mutate(pct_site = summary(model)$varcor %>%
           as_tibble() %>%
           mutate(pct_var = vcov / sum(vcov)) %>%
           filter(grp == 'site') %>%
           pull(pct_var)) %>%
  ungroup



field_out <- select(field_asv_models, asv_id, tp_contrasts) %>%
  unnest(tp_contrasts) %>%
  group_by(year, season) %>%
  mutate(fdr = p.adjust(p.value, 'fdr')) %>% 
  ungroup %>%
  mutate(estimate = str_c(round(estimate, 1), round(SE, 2), sep = '+-'),
         estimate = if_else(fdr < 0.05, str_c(estimate, ' *'), estimate),
         .keep = 'unused') %>%
  select(-df, -t.ratio, -contains('p.value')) %>%
  pivot_wider(names_from = c(year, season),
              values_from = estimate) %>%
  relocate(`2017_S`, .after = `2016_W`) %>%
  select(-contrast)



tank_out <- read_rds('../../intermediate_files/tank_asv_models.rds.gz') %>% 
  filter(asv_id %in% ml_model_out$asv_id) %>%
  inner_join(mutate(field_asv_models, across(contains('p.within'), ~p.adjust(., method = 'fdr'))) %>%
               filter_at(vars(contains('p.within')), all_vars(. < alpha)) %>%
               select(asv_id),
             by = 'asv_id') %>%
  select(asv_id, tank_posthoc) %>%
  unnest(tank_posthoc) %>%
  group_by(contrast) %>%
  mutate(fdr = p.adjust(p.value, 'fdr')) %>% 
  ungroup %>%
  mutate(estimate = str_c(round(estimate, 1), round(std.error, 2), sep = '+-'),
         estimate = if_else(fdr < 0.05, str_c(estimate, ' *'), estimate),
         .keep = 'unused') %>%
  select(asv_id, contrast, estimate) %>%
  pivot_wider(names_from = c(contrast),
              values_from = estimate) %>%
  mutate(likely_type = case_when(str_detect(DDvDH, '\\*') & 
                                   str_detect(DvH, '\\*') & 
                                   str_detect(PostvPreD, '\\*') ~ 'Pathogen',
                                 
                                 str_detect(DvN, '\\*') & 
                                   str_detect(PostvPreD, '\\*') ~ 'Opportunist',
                        
                        TRUE ~ 'Commensalist'))


asv_table <- select(ml_model_out, phylum:species, 
                    asv_id, median_rank, FDR) %>%
  relocate(`asv_id`, .after = `species`) %>%
  left_join(select(taxonom_confidence, asv_id, ends_with('confidence')) %>%
              select(-domain_confidence), 
            by = 'asv_id') %>%
  
  rename_with(~str_c(., '_name'), .cols = phylum:species) %>%
  pivot_longer(cols = c(ends_with('name'), ends_with('confidence')),
               names_to = c('taxon_level', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  mutate(taxon_level = str_to_sentence(taxon_level),
         name = str_c(name, ' (', scales::percent(confidence, scale = 1), ')'),
         .keep = 'unused') %>%
  pivot_wider(names_from = taxon_level, values_from = name) %>%
  relocate(asv_id:FDR, .after = Species) %>%
  
  full_join(field_out,
            by = 'asv_id') %>%
  relocate(`2016_W`, .before = `2016_S`) %>%
  relocate(`2017_W`, .before = `2017_S`) %>%
  left_join(tank_out,
            by = 'asv_id') %>%
  relocate(`2016_W`, .before = `2016_S`) %>%
  relocate(`2017_W`, .before = `2017_S`) %>%
  rename(ID = asv_id) %>%
  mutate(across(where(is.character), ~str_replace_all(., '\\+-', ' ± '))) %>%
  rename('Median Rank' = median_rank) %>%
  relocate(`PostvPreD`, .after = `2017_S`) %>%
  relocate(`PostvPreH`, .after = `2017_S`) %>%
  relocate(`DvN`, .after = `PostvPreD`) 

write_csv(asv_table, '../../Results/Table45_asv_table.csv')  



#### Random Numbers ####
filter(field_asv_models, asv_id %in% c('ASV15', 'ASV49', 'ASV700')) %>%
  select(asv_id, tp_contrasts) %>%
  unnest(tp_contrasts) %>%
  filter(p.value > 0.05)

filter(field_asv_models, asv_id == 'ASV40') %>%
  select(model) %>%
  pull(model) %>%
  pluck(1) %>%
  emmeans(~health) %>%
  contrast('pairwise')

filter(tank_mid, asv_id == 'ASV40') %>%
  select(time:p.value)

filter(ml_model_out, !asv_id %in% c('ASV15', 'ASV49', 'ASV700', 'ASV40')) %>%
  count(family, genus)


filter(tank_mid, asv_id %in% str_c('ASV', c(25, 322, 361, 26, 8))) %>%
  select(asv_id, time:p.value) %>%
  filter(p.value < 0.05)

filter(taxonomy, asv_id == 'ASV8') %>%
  select(family, genus)

filter(taxonomy, asv_id %in% c('ASV26', 'ASV361', 'ASV25')) %>%
  select(family, genus, asv_id)
