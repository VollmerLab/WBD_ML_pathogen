# cpm ~ (control-healthyEnd, disease-healthyEnd, disease-diseaseEnd) * dosed/not 

#### Libraries ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
# library(rstanarm)
library(lubridate)
library(emmeans)
library(broom)
library(patchwork)
library(multidplyr)

#### N ASVs ####
alpha <- 0.05
rerun_models <- FALSE

#### Data ####
the_nmds <- read_rds('../intermediate_files/field_tank_nmds.rds.gz')
asv_fit <- read_rds('../intermediate_files/field_tank_asvArrows.rds.gz')

model_list <- read_csv('../Results/equivilant_top_models.csv.gz',
                       show_col_types = FALSE) %>%
  filter(pract_equiv >= 0.95) %>%
  rename(wflow_id = model) %>%
  pull(wflow_id)

taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz',
                     show_col_types = FALSE) %>%
  mutate(across(where(is.character), 
                ~str_replace_na(., replacement = '')))

shap_importance <- read_csv('../Results/asv_importance.csv.gz',
                            show_col_types = FALSE) %>%
  filter(p_adjust < alpha) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(asv_id = fct_reorder(asv_id, median_rank)) %>%
  select(asv_id, domain:genus) %>%
  distinct %>%
  arrange(asv_id)

field_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE)

whacky <- c('Bin1_N_D_YE', 'Bin2_N_D_YE', 'Bin5_N_D_OR', 'Bin5_N_D_GR')
tank_asv_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
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
         -norm.factors, -domain:-species) %>%
  
  #Get rid of samples which shouldnt exist
  filter(!(str_detect(sample_id, str_c(whacky, collapse = '|')) & time == 8))


tank_asv_data %>%
  count(tank_id, time) %>%
  pivot_wider(names_from = time, values_from = n)

tank_asv_data %>%
  count(geno, time) %>%
  pivot_wider(names_from = time, values_from = n)

#### NMDS ####
species_fits <- scores(the_nmds)$species %>%
  as_tibble(rownames = 'asv_id') %>%
  left_join(tibble(asv_id = names(asv_fit$vectors$r),
                   r2 = asv_fit$vectors$r,
                   p = asv_fit$vectors$pvals),
            by = 'asv_id') %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(useful = str_c(family, genus, asv_id, sep = '_') %>%
           str_replace_all('__', '_')) %>% 
  filter(family %in% unique(shap_importance$family)) %>%
  mutate(important_asv = asv_id %in% unique(shap_importance$asv_id)) %>%
  identity()
  

species_fits %>%
  arrange(important_asv) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_segment(aes(colour = important_asv),
               xend = 0, yend = 0, show.legend = FALSE) +
  geom_point(aes(colour = important_asv, size = important_asv)) +
  geom_point(data = scores(the_nmds)$sites %>%
               as_tibble(rownames = 'sample_id') %>%
               left_join(select(field_data, sample_id, health) %>%
                           distinct(),
                         by = 'sample_id'),
             aes(fill = health),
             size = 2, shape = 21) +
  geom_text(data = . %>% filter(important_asv),
            aes(label = useful), show.legend = FALSE) +
  facet_wrap(~family) +
  scale_colour_manual(values = c('TRUE' = 'black', 'FALSE' = 'grey50')) +
  scale_size_manual(values = c('TRUE' = 5, 'FALSE' = 1)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = 'black')),
         size = 'none') +
  labs(colour = 'ASV Importance',
       fill = 'Health State') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'))
ggsave('../Results/nmds_important_asvs.png', height = 15, width = 15)

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

%>%
  filter(p_adj < alpha)

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

96/sum(c(96,16,28,3))

#### Get sequences of ASVs and BLAST ####
if(!file.exists('../intermediate_files/interesting_asv_16s.fasta')){
  read_rds('../intermediate_files/prepped_microbiome.rds.gz') %>%
    prune_taxa(taxa_names(.) %in% as.character(shap_importance$asv_id), .) %>%
    refseq %>%
    Biostrings::writeXStringSet(filepath = '../intermediate_files/interesting_asv_16s.fasta')
}

## BLAST on NCBI ##

if(file.exists('../intermediate_files/interesting_asv_blast_classifications.csv.gz')){
  the_taxa_complete <- read_csv('../intermediate_files/interesting_asv_blast_classifications.csv.gz', show_col_types = FALSE)
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
  
  the_taxa <- read_csv('../intermediate_files/UP5UBP83013-Alignment-HitTable.csv',
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
  
  
  write_csv(the_taxa_complete, '../intermediate_files/interesting_asv_blast_classifications.csv.gz')
}

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

#### Model ASVs in Field ####
# model <- field_asv_models$model[[1]]
# data <- field_asv_models$data[[1]]

plot_field <- function(model, data, caption = FALSE, alpha = 0.05, alpha_filter = 1){
  if(all(class(model) == 'lm')){
    plot_data <- emmeans(model, ~health * year * season * site,
                         data = data) %>%
      as_tibble()
  }
  
  if(all(class(model) == 'lmerModLmerTest') | 'stanreg' %in% class(model)){
    plot_data <- emmeans(model, ~health * year * season,
                         data = data) %>%
      as_tibble()
  }
  
  the_caption <- NULL
  if(caption){
    if('stanreg' %in% class(model)){
      the_caption <- NULL
    } else {
      if(all(class(model) == 'lm')){
        the_aov <- anova(model)
        
        the_aov <- as_tibble(the_aov, rownames = 'effect') %>%
          filter(effect != 'Residuals') %>%
          rename(NumDF = Df) %>%
          mutate(DenDF = the_aov$Df[length(the_aov$Df)])
        
      } else {
        the_aov <- anova(model, ddf = 'Kenward-Roger') %>%
          as_tibble(rownames = 'effect') 
      }
      the_caption <- the_aov %>%
        filter(`Pr(>F)` < alpha_filter) %>%
        mutate(effect = str_replace_all(effect, c('_c' = '', '_' = ' ', ':' = ' x ')),
               effect = str_to_title(effect) %>% str_replace('X', 'x'),
               across(c(where(is.numeric), -`F value`, -`Pr(>F)`), ~round(., 1)),
               `F value` = round(`F value`, 3),
               output = str_c(effect, ': F(', NumDF, ', ', DenDF, ') = ', 
                              `F value`, '; p = ', scales::pvalue(`Pr(>F)`))) %>%
        pull(output) %>%
        str_c(collapse = '\n')
    }
  }
  
  out <- plot_data %>%
    rename_with(~str_replace(., 'HPD', 'CL')) %>%
    mutate(grouping = row_number()) %>%
    mutate(timepoint = if_else(season == 'S', 
                               ymd(str_c(year, 'July', 1, sep = '-')),
                               ymd(str_c(year, 'January', 1, sep = '-')))) %>%
    ggplot(aes(x = timepoint, y = emmean, colour = health,
               ymin = lower.CL, ymax = upper.CL, group = grouping)) +
    geom_pointrange(position = position_dodge(50)) +
    labs(caption = the_caption)
  
  if(all(class(model) == 'lmerModLmerTest')){
    timepoint_diff <- emmeans(model, ~health | year + season) %>%
      contrast('pairwise') %>%
      as_tibble() %>%
      mutate(grouping = row_number(),
             timepoint = if_else(season == 'S', 
                                 ymd(str_c(year, 'July', 1, sep = '-')),
                                 ymd(str_c(year, 'January', 1, sep = '-')))) 
    
    out <- out +
      geom_text(data = timepoint_diff, inherit.aes = FALSE,
                aes(x = timepoint, y = Inf, 
                    label = if_else(p.value < alpha, '*', '')),
                vjust = 1, size = 8)
  }
  
  out
}

extract_significance <- function(model, anova_summary = 'summary'){
  if('lmerModLmerTest' %in% class(model)){
    if(anova_summary == 'summary'){
      anova_table <- summary(model, ddf = 'Kenward-Roger')$coefficients
    } else if(anova_summary == 'anova'){
      anova_table <- anova(model, ddf = 'Kenward-Roger')
    }
    
    out <- as_tibble(anova_table, rownames = 'param') %>%
      select(param, starts_with('Pr')) %>%
      pivot_wider(names_from = param, 
                  values_from = starts_with('Pr'),
                  names_prefix = 'p.') %>%
      rename_with(~str_replace_all(., c(':' = 'X', '_' = '')) %>%
                    str_remove_all('\\(|\\)'))
    
    out$summary_table <- list(anova_table)
  }
  
  if('stanreg' %in% class(model)){
    anova_table <- summary(model, 
                           pars = names(fixef(model)),
                           probs = c(0.1, 0.9)) %>%
      as_tibble(rownames = 'param') %>%
      mutate(sig = (`10%` < 0 & `90%` < 0) | (`10%` > 0 & `90%` > 0))
    
    out <- select(anova_table, param, sig) %>%
      pivot_wider(names_from = param, 
                  values_from = sig,
                  names_prefix = 'p.') %>%
      rename_with(~str_replace_all(., c(':' = 'X', '_' = '')) %>%
                    str_remove_all('\\(|\\)'))
    
    out$summary_table <- list(anova_table)
  }
  
  if(all(class(model) == 'lm')){
    if(anova_summary == 'summary'){
      anova_table <- summary(model)$coefficients
    } else if(anova_summary == 'anova'){
      anova_table <- anova(model)
    }
    
    out <- as_tibble(anova_table, rownames = 'param') %>%
      select(param, starts_with('Pr')) %>%
      pivot_wider(names_from = param, 
                  values_from = starts_with('Pr'),
                  names_prefix = 'p.') %>%
      rename_with(~str_replace_all(., c(':' = 'X', '_' = '')) %>%
                    str_remove_all('\\(|\\)'))
    
    out$summary_table <- list(anova_table)
  }
  
  out
}


# model <- field_asv_models$model[[1]]; data <- field_asv_models$data[[1]]; form <- ~health | season * year
# model <- tank_asv_models$model[[1]]; data <- tank_asv_models$data[[1]]; form <- ~health | exp_dis + exp_hea

make_coef_table <- function(form, model, data){
  
  x_vars <- str_split(as.character(form)[2], ' [\\+\\*] ') %>% unlist %>% unique %>%
    str_remove_all('\\(|\\)')
  
  out <- emmeans(model, form, data = data) %>%
    as_tibble %>%
    select(all_of(x_vars), emmean)
  
  if(any(x_vars == 'exp_dis')){
    out <- out %>%
      filter(!(exp_dis == 1 & exp_hea == 1)) %>%
      mutate(treatment = case_when(exp_dis == 0 & exp_hea == 0 ~ str_c('pre', health, sep = '_'),
                                   exp_dis == 1 ~ str_c('postDisease', health, sep = '_'),
                                   exp_hea == 1 ~ str_c('postHealth', health, sep = '_'),)) %>%
      select(treatment, emmean)
    x_vars <- 'treatment'
  } 
  
  out %>%
    pivot_wider(names_from = all_of(x_vars),
                values_from = emmean,
                names_prefix = 'mu_')
}

make_fc_table <- function(form, model, data){
  x_vars <- str_split(as.character(form)[2], ' [\\+\\*\\|] ') %>% unlist %>% unique %>%
    str_remove_all('\\(|\\)')
  
  out <- emmeans(model, form, data = data) %>%
    contrast(if_else(any(x_vars == 'exp_dis'), 'revpairwise', 'pairwise')) %>%
    as_tibble %>%
    select(any_of(x_vars), estimate)
  
  if(any(x_vars == 'exp_dis')){
    out <- out %>%
      filter(!(exp_dis == 1 & exp_hea == 1)) %>%
      mutate(treatment = case_when(exp_dis == 0 & exp_hea == 0 ~ 'pre',
                                   exp_dis == 1 ~ 'postDisease',
                                   exp_hea == 1 ~ 'postHealth')) %>%
      select(treatment, estimate)
    x_vars <- 'treatment'
  } 
  
  out %>%
    pivot_wider(names_from = any_of(x_vars),
                values_from = estimate,
                names_prefix = 'fcDH_')
}

test_within_timepoint <- function(model, data){
  tp_contrasts <- emmeans(model, ~health | year * season, data = data) %>%
    contrast('pairwise') %>%
    as_tibble
  
  select(tp_contrasts, year, season, p.value) %>%
    pivot_wider(names_from = c(year, season),
                values_from = p.value,
                names_prefix = 'p.withinYear_') %>%
    mutate(tp_contrasts = list(tp_contrasts))
  
}

cluster <- new_cluster(parallel::detectCores() / 2)
cluster_library(cluster, c('rstanarm', 'lmerTest', 'emmeans', 'dplyr', 
                           'tidyr', 'tibble', 'ggplot2', 'stringr',
                           'lubridate'))
cluster_copy(cluster, c('plot_field', 'extract_significance',
                        'make_coef_table', 'test_within_timepoint',
                        'make_fc_table'))

if(file.exists('../intermediate_files/field_asv_models.rds.gz') & !rerun_models){
  field_asv_models <- read_rds('../intermediate_files/field_asv_models.rds.gz')
} else {
  field_asv_models <- field_data %>%
    mutate(year = as.character(year)) %>%
    select(asv_id, sample_id, health, year, season, site, log2_cpm_norm, n_reads, lib.size) %>%
    nest(data = -c(asv_id)) %>%
    # filter(asv_id == 'ASV25') %>%
    rowwise %>%
    mutate(data = list(mutate(data, log2_cpm_norm = log2_cpm_norm + rnorm(nrow(data), 0, 6.11/100)))) %>%
    partition(cluster) %>%
    
    mutate(model = list(lmer(log2_cpm_norm ~ health * year * season + (1 | site),
                           data = data)),
         plot = list(plot_field(model, data, caption = TRUE) +
                       labs(y = 'log2(CPM)', x = NULL) +
                       theme_classic() +
                       theme(legend.position = 'none',
                             panel.background = element_rect(colour = 'black'),
                             axis.text = element_text(colour = 'black'))),
         extract_significance(model, anova_summary = 'anova'),
         
         average_health_fit = list(emmeans(model, ~health)),
         just_health_fit = list(emmeans(model, ~health, at = list(year = '2017', season = 'W'))),
         all_health_fit = list(emmeans(model, ~health | year * season)),
         
         test_within_timepoint(model, data),
         make_coef_table(~health * year * season, model, data),
         make_fc_table(~health | year * season, model, data)) %>%
    
    
    collect %>%
    ungroup 
  write_rds(field_asv_models, '../intermediate_files/field_asv_models.rds.gz')
}

field_asv_models %>%
  # filter(asv_id == 'ASV25') %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>%
  # mutate(across(starts_with('p.'), ~p.adjust(., 'fdr'))) %>%
  # filter(if_any(c(starts_with('p.'), contains('health')), ~. < alpha)) %>%
  
  left_join(taxonomy, by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(plot + labs(title = str_c(order, family, genus, asv_id, sep = '\n'),
                                 caption = NULL) +
                       theme(plot.title = element_text(size = 10)))) %>%
  ungroup %>%
  pull(plot) %>%
  wrap_plots(byrow = TRUE)
ggsave('../Results/important_asvs_field.png', height = 15, width = 15)


#### Model ASVs in tanks ####
# model<-tank_asv_models$model[[1]]; data <- tank_asv_models$data[[1]]; field_fitted <- tank_asv_models$just_health_fit[[1]]
# field_fitted <- field_asv_models$all_health_fit[[1]]

plot_prePost <- function(model, data, field_fitted){
  initial_emmean <- emmeans(model, ~time_treat)
  treatment_levels <- levels(initial_emmean)$time_treat
  
  
  ref_grid(model) %>%
    add_grouping('treat_out', 'time_treat', 
                 str_remove(treatment_levels, '[0-9]+_')) %>%
    emmeans(~ treat_out) %>%
    broom::tidy(conf.int = TRUE) %>%
    mutate(prePost = if_else(str_detect(treat_out, 'pre'), 'Pre', 'Post'),
           prePost = factor(prePost, levels = c('Field', 'Pre', 'Post')),
           
           exposure = str_extract(treat_out, '^[DN]'),
           health = str_extract(treat_out, '[DH]$')) %>%
    bind_rows(broom::tidy(field_fitted, conf.int = TRUE) %>%
                mutate(prePost = factor('Field', levels = c('Field', 'Pre', 'Post')))) %>%
    ggplot(aes(x = prePost, y = estimate, colour = health, shape = exposure)) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    position = position_dodge2(0.5)) +
    scale_shape_manual(values = c('D' = 'circle', 'N' = 'triangle'), na.value = 'square')
  
}

tank_posthocs <- function(model){
  #Significance between healthy vs disease after exposure
  #significance of diseased compared to before exposure
  #significance of healthy compared to before exposure
  initial_emmean <- emmeans(model, ~time_treat,  
                            lmer.df = 'kenward-roger')
  treatment_levels <- levels(initial_emmean)$time_treat
  
  all_contrasts <- list('DvH' = case_when(str_detect(treatment_levels, '^.*_D$') ~ 1/2, #diseased vs healthy
                                          str_detect(treatment_levels, '^[28].*_H$') ~ -1/4,
                                          TRUE ~ 0),
                        
                        'DvN' = case_when(str_detect(treatment_levels, '^.*_D.*$') ~ 1/4, #diseased vs healthy exposure
                                          str_detect(treatment_levels, '^.*_N.*$') ~ -1/2,
                                          TRUE ~ 0),
                        
                        'DDvDH' = case_when(str_detect(treatment_levels, '^.*_D_D$') ~ 1/2,
                                            str_detect(treatment_levels, '^.*_D_H$') ~ -1/2,
                                            TRUE ~ 0),
                        
                        'DDvNH' = case_when(str_detect(treatment_levels, '^.*_D_D$') ~ 1/2,
                                            str_detect(treatment_levels, '^.*_N_H$') ~ -1/2,
                                            TRUE ~ 0),
                        
                        'DHvNH' = case_when(str_detect(treatment_levels, '^.*_D_H$') ~ 1/2,
                                            str_detect(treatment_levels, '^.*_N_H$') ~ -1/2,
                                            TRUE ~ 0),
                        
                        'PostvPreD' = case_when(str_detect(treatment_levels, '0_pre_H') ~ -1,
                                                str_detect(treatment_levels, '^.*_D$') ~ 1/2,
                                                TRUE ~ 0),
                        
                        'PostvPreH' = case_when(str_detect(treatment_levels, '0_pre_H') ~ -1,
                                                str_detect(treatment_levels, '^.*_H$') ~ 1/4,
                                                TRUE ~ 0))
  
  posthoc_sided <- contrast(initial_emmean, all_contrasts, side = '>', adjust = 'none')
  full_posthoc <- contrast(initial_emmean, all_contrasts, side = '=', adjust = 'none') %>%
    broom::tidy(conf.int = TRUE) %>%
    select(contrast, estimate, std.error, conf.low, conf.high) %>%
    left_join(broom::tidy(posthoc_sided) %>%
                select(contrast, df, statistic, p.value),
              by = 'contrast')
  
  select(full_posthoc, contrast, estimate, p.value) %>%
    rename(fc = estimate, pvalue = p.value) %>%
    pivot_wider(names_from = contrast, values_from = c(fc, pvalue)) %>%
    mutate(tank_posthoc = list(full_posthoc))
}

cluster_copy(cluster, c('plot_prePost', 'tank_posthocs'))

if(file.exists('../intermediate_files/tank_asv_models.rds.gz') & !rerun_models){
  tank_asv_models <- read_rds('../intermediate_files/tank_asv_models.rds.gz')
} else {
  tank_asv_models <- tank_asv_data %>%
    nest_by(asv_id) %>%
    left_join(select(field_asv_models, asv_id, all_health_fit, starts_with('p.withinYear')),
              by = 'asv_id') %>%
    # filter(asv_id == 'ASV25') %>%
    arrange(asv_id) %>%
    rowwise %>%
    # mutate(data = list(mutate(data, log2_cpm_norm = log2_cpm_norm + rnorm(nrow(data), 0, 6.11/100)))) %>%
    partition(cluster) %>%
    mutate(model = list(lmer(log2_cpm_norm ~ time_treat + (1 | geno) + (1 | tank_id),
                             data = data))) %>%
    
    mutate(plot = list(plot_prePost(model, data, field_fitted = all_health_fit) +
                         labs(y = 'log2(CPM)', x = 'Days') +
                         theme_classic() +
                         theme(legend.position = 'none',
                               panel.background = element_rect(colour = 'black'),
                               axis.text = element_text(colour = 'black'),
                               strip.background = element_blank())),
           
           extract_significance(model, anova_summary = 'anova'),
           tank_posthocs(model)) %>%
    collect() %>%
    ungroup
  
  write_rds(tank_asv_models, '../intermediate_files/tank_asv_models.rds.gz')
}

tank_asv_models$plot[[1]]
tank_asv_models %>%
  
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>%
  
  # inner_join(filter(field_asv_models,
  #                   if_all(starts_with('p.withinYear'), ~. < alpha),
  #                   asv_id %in% shap_importance$asv_id) %>%
  #              select(asv_id),
  #            by = 'asv_id') %>%
  
  # select(-contains('Intercept')) %>%
  # mutate(across(starts_with('p.'), ~p.adjust(., 'fdr'))) %>%
  # filter(if_any(c(starts_with('p.') & contains('SX')), ~. < alpha)) %>%
  # filter(if_any(c(starts_with('p.') & contains('SX')))) %>%
  # filter(p.healthSXexpdis < alpha) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(plot + labs(title = str_c(order, family, genus, asv_id, sep = '\n'))+
                       theme(plot.title = element_text(size = 10)))) %>%
  ungroup %>%
  pull(plot) %>%
  wrap_plots()


classified_plots <- tank_asv_models %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(across(contains('p.within'), ~p.adjust(., method = 'fdr'))) %>%
  filter_at(vars(contains('p.within')), all_vars(. < alpha)) %>%
  # filter(p.adjust(p.timetreat, method = 'fdr') < alpha) %>%
  mutate(across(contains('pvalue_'), ~p.adjust(., method = 'fdr'))) %>%
  # filter_at(vars(contains('pvalue_')), any_vars(. < alpha)) %>%

  mutate(likely_type = case_when(pvalue_DDvDH < alpha & 
                                   #pvalue_DDvNH < alpha & 
                                   pvalue_DvH < alpha &
                                   pvalue_PostvPreD < alpha ~ 'Pathogen',
                                 
                                 pvalue_DvN < alpha &
                                   pvalue_PostvPreD < alpha ~ 'Opportunist',
                                 
                                 TRUE ~ 'Commensalist'),
         likely_type = factor(likely_type, levels = c('Pathogen', 'Opportunist', 'Commensalist')),
         asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>%
  
  # select(asv_id, starts_with('pvalue'), likely_type)
  
  left_join(taxonomy, by = 'asv_id') %>%
  rowwise %>%
  mutate(plot = list(plot + labs(title = str_c(order, family, genus, species, asv_id, sep = '\n'))+
                       theme(plot.title = element_text(size = 10)))) %>%
  ungroup %>%
  arrange(likely_type) %>%
  group_by(likely_type) %>%
  summarise(plot_set = list(wrap_plots(plot, nrow = 1, tag_level = 'new') + plot_annotation(title = likely_type))) %>%
  pull(plot_set) %>%
  wrap_plots(nrow = 3, tag_level = 'new')
  
classified_plots
ggsave('../Results/Fig8_important_asvs.png', height = 15, width = 15)

classified_plots +
  plot_annotation(tag_levels = list('A', c(NULL)))


#### P/A % plots ####
raw_pa_data <- tank_asv_models %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(across(contains('p.within'), ~p.adjust(., method = 'fdr'))) %>%
  filter_at(vars(contains('p.within')), all_vars(. < alpha)) %>%
  # filter(p.adjust(p.timetreat, method = 'fdr') < alpha) %>%
  mutate(across(contains('pvalue_'), ~p.adjust(., method = 'fdr'))) %>%
  # filter_at(vars(contains('pvalue_')), any_vars(. < alpha)) %>%
  
  mutate(likely_type = case_when(pvalue_DDvDH < alpha & 
                                   #pvalue_DDvNH < alpha & 
                                   pvalue_DvH < alpha &
                                   pvalue_PostvPreD < alpha ~ 'Pathogen',
                                 
                                 pvalue_DvN < alpha &
                                   pvalue_PostvPreD < alpha ~ 'Opportunist',
                                 
                                 TRUE ~ 'Commensalist'),
         likely_type = factor(likely_type, levels = c('Pathogen', 'Opportunist', 'Commensalist')),
         asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>%
  
  # select(asv_id, starts_with('pvalue'), likely_type)
  
  left_join(taxonomy, by = 'asv_id') %>%
  select(asv_id, domain:species, likely_type, data) %>%
  rename(tank_data = data) %>%
  left_join(select(field_asv_models, asv_id, data) %>%
              rename(field_data = data),
            by = 'asv_id') %>%
  rowwise(asv_id:likely_type) %>%
  mutate(field_data = list(mutate(field_data, treatment = str_c(health, year, season, sep = '_')) %>%
                             select(sample_id, n_reads, treatment)),
         tank_data = list(select(tank_data, sample_id, n_reads, treatment))) %>%
  
  reframe(bind_rows(tank_data, field_data)) %>%
  group_by(across(c(asv_id:likely_type, treatment))) %>%
  summarise(n_present = sum(n_reads > 0),
            n_missing = sum(n_reads == 0),
            total = n(),
            .groups = 'rowwise') %>% 
  reframe(broom::tidy(binom.test(x = n_present, n = total))) %>%
  mutate(health = case_when(str_detect(treatment, '201[67]') ~ str_extract(treatment, '[HD]'),
                            treatment == 'D_D' ~ 'D',
                            TRUE ~ 'H'),
         x_val = case_when(str_detect(treatment, '201[67]') ~ 'Field',
                           str_detect(treatment, '[Pp]re') ~ 'Pre',
                           TRUE ~ 'Post'),
         x_val = factor(x_val, levels = c('Field', 'Pre', 'Post')),
         exposure = case_when(x_val %in% c('Field', 'Pre') ~ 'none',
                              treatment == 'N_H' ~ 'H',
                              TRUE ~ 'D'))

raw_pa_data %>%
  mutate(asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  ggplot(aes(x = x_val, y = estimate, ymin = conf.low, ymax = conf.high,
             colour = health, group = treatment, shape = exposure)) +
  geom_pointrange(position = position_dodge(0.5)) +
  facet_wrap(likely_type~asv_id)
ggsave('../Results/percent_presence_important_asvs.png', height = 15, width = 15)

#### Potential pathogens which samples are missing them ####
the_pathogens <- tank_asv_models %>%
  filter(asv_id %in% shap_importance$asv_id) %>%
  mutate(across(contains('p.within'), ~p.adjust(., method = 'fdr'))) %>%
  filter_at(vars(contains('p.within')), all_vars(. < alpha)) %>%
  # filter(p.adjust(p.timetreat, method = 'fdr') < alpha) %>%
  mutate(across(contains('pvalue_'), ~p.adjust(., method = 'fdr'))) %>%
  # filter_at(vars(contains('pvalue_')), any_vars(. < alpha)) %>%
  
  mutate(likely_type = case_when(pvalue_DDvDH < alpha & 
                                   #pvalue_DDvNH < alpha & 
                                   pvalue_DvH < alpha &
                                   pvalue_PostvPreD < alpha ~ 'Pathogen',
                                 
                                 pvalue_DvN < alpha &
                                   pvalue_PostvPreD < alpha ~ 'Opportunist',
                                 
                                 TRUE ~ 'Commensalist'),
         likely_type = factor(likely_type, levels = c('Pathogen', 'Opportunist', 'Commensalist')),
         asv_id = factor(asv_id, levels = levels(shap_importance$asv_id))) %>%
  arrange(asv_id) %>%
  
  # select(asv_id, starts_with('pvalue'), likely_type)
  
  left_join(taxonomy, by = 'asv_id') %>%
  select(asv_id, domain:species, likely_type, data) %>%
  rename(tank_data = data) %>%
  left_join(select(field_asv_models, asv_id, data) %>%
              rename(field_data = data),
            by = 'asv_id') %>%
  rowwise(asv_id:likely_type) %>%
  mutate(field_data = list(mutate(field_data, treatment = str_c(health, year, season, sep = '_')) %>%
                             select(sample_id, n_reads, lib.size, treatment)),
         tank_data = list(select(tank_data, sample_id, n_reads, lib.size, treatment))) %>%
  
  reframe(bind_rows(tank_data, field_data)) %>% #lib.size
  mutate(health = case_when(str_detect(treatment, '201[67]') ~ str_extract(treatment, '[HD]'),
                            treatment == 'D_D' ~ 'D',
                            TRUE ~ 'H')) %>%
  filter(health == 'D',
         likely_type == 'Pathogen')

the_pathogens %>%
  filter(asv_id == 'ASV25') %>%
  glmer(as.integer(n_reads > 0) ~ scale(lib.size) + (1 | treatment), data = .,
      family = binomial(link = 'logit')) %>%
  emmeans(~lib.size, at = list(lib.size = seq(0, 3000, length.out = 1000)),
          type = 'response') %>%
  as_tibble() %>%
  ggplot(aes(x = lib.size, y = response)) +
  geom_line()
  
  
  car::Anova()
  ggplot(aes(x = lib.size, y = as.integer(n_reads > 0), colour = treatment)) +
  # stat_summary_bin(bins = 100)
  geom_point()

#### Output asv coefficients ####
asv_coef <- full_join(select(field_asv_models, asv_id, starts_with('mu'), starts_with('fcDH')),
                      select(tank_asv_models, asv_id, starts_with('mu'), starts_with('fc')),
                      by = 'asv_id') %>%
  select(asv_id, starts_with('fcDH'), starts_with('mu_'))

write_csv(asv_coef, '../Results/asv_coefficients.csv.gz') 

the_taxa_complete