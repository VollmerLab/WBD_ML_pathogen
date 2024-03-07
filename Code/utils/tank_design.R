library(tidyverse)


#### Data ####
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


#### Document tank outcomes ####
n_distinct(tank_data$geno)
count(tank_data, treatment)

mutate(tank_data,
       treatment = str_remove(treatment, '_[DH]$')) %>%
  group_by(time, treatment, geno) %>%
  summarise(n_frag = n_distinct(fragment_id),
            .groups = 'drop_last') %>%
  summarise(n_geno = n_distinct(geno),
            mean_frag = mean(n_frag),
            se_frag = sd(n_frag) / sqrt(n()),
            n_frag = sum(n_frag))

mutate(tank_data,
       tank = str_remove(tank, '_.*$'),
       treatment = str_remove(treatment, '_[DH]$')) %>%
  group_by(tank, time, treatment, health) %>%
  summarise(n_geno = n_distinct(geno), 
            .groups = 'drop') %>%
  pivot_wider(names_from = health,
              values_from = n_geno,
              values_fill = 0L) %>%
  mutate(genos = str_c(H, 'H / ', D, 'D'), 
         .keep = 'unused') %>%
  pivot_wider(names_from = time,
              values_from = genos)

tank_data %>%
  filter(treatment == 'D_D') %>% 
  count(geno)


tank_data %>%
  filter(treatment %in% c('D_D', 'D_H')) %>%
  count(geno, health) %>%
  pivot_wider(names_from = health,
              values_from = n,
              values_fill = 0L)

tank_data %>%
  count(geno, time, health, fragment_id) %>%
  pivot_wider(names_from = time, values_from = n) %>% View
