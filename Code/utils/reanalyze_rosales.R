library(tidyverse)

#### Data ####
asv_of_interest <- tibble(asv_id = c('5b79cf6d5a5a9bf0bb866aed449eff44',
                                     '7eb68c2ff12bb8a0a46d036c37f8f26e'),
                          our_id = c('ASV25', 'ASV8'))



metadata <- read_tsv('../../Data/rosales/DR_MB1_meta_seqID_run1.txt', 
                     show_col_types = FALSE) %>% 
  select(SampleID, scientific_name, ASV_Frequency, 
         Genotype, replicate, Genotype_replicate, 
         Outcome, Treatment) %>%
  rename(sample_id = SampleID)

acer_asv_counts <- read_tsv('../../Data/rosales/table-fil-BacArc_DRMB1AC.txt', 
                            show_col_types = FALSE) %>%
  rename(asv_id = SampleID) %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'n_read') %>%
  mutate(total_reads = sum(n_read), .by = 'sample_id')  %>%
  inner_join(asv_of_interest,
             by = 'asv_id')

apalm_asv_counts <- read_tsv('../../Data/rosales/table-fil-BacArc_DRMB1AP.txt', 
                             show_col_types = FALSE) %>%
  rename(asv_id = SampleID) %>%
  pivot_longer(cols = -asv_id,
               names_to = 'sample_id',
               values_to = 'n_read') %>%
  mutate(total_reads = sum(n_read), .by = 'sample_id') %>%
  inner_join(asv_of_interest,
             by = 'asv_id')


#### Join Metadata & ASV data ####
full_data <- full_join(metadata,
          bind_rows(acer_asv_counts,
                    apalm_asv_counts),
          by = 'sample_id') %>%
  filter(!is.na(n_read))


#### Simple Plot ####
full_data %>%
  ggplot(aes(x = our_id, y = n_read, colour = scientific_name)) +
  geom_point(position = position_dodge(0.5))


#### split by species & asv id and plot by outcome ####
full_data %>%
  ggplot(aes(x = Outcome, y = n_read / total_reads)) +
  stat_summary(fun.data = mean_se) +
  facet_grid(our_id ~ scientific_name,
             scales = 'free')


#### Analysis ####
library(lme4)
library(glmmTMB)
library(emmeans)

all_models <- full_data %>%
  nest_by(scientific_name, our_id) %>%
  mutate(model = list(glmmTMB(n_read ~ Outcome + (1 | Genotype) + 
                              offset(log(total_reads)),
                            family = 'nbinom2',
                            data = data)),
         emout = list(emmeans(model, ~Outcome, type = 'response', 
                              offset = log(1000)) ),
         pairwise_posthoc = list(update(emout, type = 'link') %>%
                                   contrast('pairwise')))


map(all_models$model, ~car::Anova(.))


#### plot model outcome ####
all_models %>%
  reframe(as_tibble(emout)) %>%
  ggplot(aes(x = Outcome, y = response)) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE)) +
  geom_point() +
  facet_grid(our_id ~ scientific_name,
             scales = 'free') +
  labs(y = '# / 1,000 Reads')

