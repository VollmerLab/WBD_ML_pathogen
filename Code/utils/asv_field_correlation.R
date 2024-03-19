library(tidyverse)

#### N ASVs ####
alpha <- 0.05

#### Data ####
taxonomy <- read_csv('../../intermediate_files/taxonomy.csv.gz',
                     show_col_types = FALSE) %>%
  mutate(across(where(is.character), 
                ~str_replace_na(., replacement = '')))

field_data <- read_csv('../../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE)

shap_importance <- read_csv('../../Results/asv_importance.csv.gz',
                            show_col_types = FALSE) %>%
  filter(p_adjust < alpha) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(asv_id = fct_reorder(asv_id, median_rank)) %>%
  select(asv_id, domain:genus) %>%
  distinct %>%
  arrange(asv_id)

#### Correlation Among Important ASVs ####
sig_correlated_asvs <- field_data %>%
  filter(asv_id %in% shap_importance$asv_id,
         health == 'D') %>%
  select(asv_id, sample_id, log2_cpm_norm) %>%
  full_join(., ., by = 'sample_id',
            relationship = "many-to-many") %>%
  filter(asv_id.x > asv_id.y) %>%
  group_by(asv_id.x, asv_id.y) %>%
  summarise(broom::tidy(cor.test(log2_cpm_norm.x, log2_cpm_norm.y)),
            .groups = 'drop') %>%
  mutate(p_adj = p.adjust(p.value, method = 'fdr')) 


sig_correlated_asvs %>%
  filter(asv_id.x == 'ASV25' | asv_id.y == 'ASV25') %>%
  arrange(-estimate)


sig_correlated_asvs %>%
  filter(asv_id.x == 'ASV25' | asv_id.y == 'ASV25')%>%
  filter(asv_id.x == 'ASV8' | asv_id.y == 'ASV8')


field_data %>%
  filter(asv_id %in% shap_importance$asv_id,
         health == 'D') %>%
  select(asv_id, sample_id, log2_cpm_norm) %>%
  full_join(., ., by = 'sample_id',
            relationship = "many-to-many") %>%
  filter(asv_id.x > asv_id.y) %>%
  filter(asv_id.x == 'ASV25' | asv_id.y == 'ASV25') %>%
  filter(asv_id.x == 'ASV8' | asv_id.y == 'ASV8') %>%
  filter(log2_cpm_norm.y > min(log2_cpm_norm.y),
         log2_cpm_norm.x > min(log2_cpm_norm.x)) %>%
  summarise(broom::tidy(cor.test(log2_cpm_norm.x, log2_cpm_norm.y, method = 'spearman')),
            .groups = 'drop') 
# 
# ggplot(aes(x = log2_cpm_norm.x, y = log2_cpm_norm.y)) +
# geom_point()


field_data %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  select(asv_id, sample_id, log2_cpm_norm, health) %>%
  full_join(., ., by = c('sample_id', 'health'),
            relationship = "many-to-many") %>%
  filter(asv_id.x > asv_id.y) %>%
  filter(asv_id.x == 'ASV25' | asv_id.y == 'ASV25') %>%
  filter(asv_id.x == 'ASV8' | asv_id.y == 'ASV8') %>%
  mutate(across(starts_with('log2'), ~. > min(.))) %>%
  count(health, log2_cpm_norm.x, log2_cpm_norm.y) %>%
  rename(ASV8 = log2_cpm_norm.x,
         ASV25 = log2_cpm_norm.y) %>%
  mutate(same_different = case_when(ASV8 & ASV25 ~ 'both',
                                    !ASV8 & !ASV25 ~ 'neither',
                                    ASV25 ~ 'only_ASV25',
                                    ASV8 ~ 'only_ASV8')) %>%
  group_by(health, same_different) %>%
  summarise(n = sum(n), .groups = 'drop') %>%
  pivot_wider(names_from = same_different,
              values_from = 'n') %>%
  relocate(neither, .after = only_ASV8) %>%
  column_to_rownames('health') %>%
  select(-both, -neither) %T>%
  print %>%
  chisq.test()
#When only one or the other is present ASV25 is more associated with field samples with disease
