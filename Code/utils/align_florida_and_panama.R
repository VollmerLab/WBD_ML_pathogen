library(Biostrings)
library(bioseq)
library(phyloseq)
library(tidyverse)


#### Get Sequences ####
panama_sequences <- read_rds('../../intermediate_files/prepped_microbiome.rds.gz') %>%
  refseq
names(panama_sequences) <- str_c('PA_', names(panama_sequences))

florida_sequences <- read_rds('../../../Emily_Microbiome/intermediate_files/preprocess_microbiome.rds') %>%
  refseq
names(florida_sequences) <- str_c('FL_ASV', 1:length(names(florida_sequences)), sep = '_')


c(panama_sequences, florida_sequences) %>%
  Biostrings::writeXStringSet(filepath = '../../intermediate_files/florida_panama_asvs.fasta')

#### Go to Discovery for cd-hit ####
# srun -t 24:00:00 --nodes=1 --cpus-per-task=24 --mem=100G --pty /bin/bash
# source activate cdhit
# cd /scratch/j.selwyn/asv_classification/
# cd-hit-est -i florida_panama_asvs.fasta -o fl_pa_clust.fasta -c 0.95 -n 10 -d 999 -M 0 -T ${SLURM_CPUS_PER_TASK}

#### Read in clustering results ####
# cd_clust <- read_lines('../../intermediate_files/fl_pa_clust.fasta.clstr')
parse_cdhit <- function(cd_clust){
  cluster_starts <- str_which(cd_clust, '^>Cluster [0-9]+')
  
  tibble(start = cluster_starts,
         end = lead(start, default = length(cd_clust) + 1),
         cluster_id = str_extract(cd_clust[start], 'Cluster [0-9]+')) %>%
    rowwise(cluster_id) %>%
    mutate(lines = list(cd_clust[(start + 1):(end - 1)])) %>%
    reframe(asv_id = str_extract(lines, '(FL|PA)_ASV[_0-9]+'))
  
}

clustered_asvs <- read_lines('../../intermediate_files/fl_pa_clust_95.fasta.clstr') %>%
  parse_cdhit()

filter(clustered_asvs, asv_id %in% c('PA_ASV8', 'PA_ASV25', 'FL_ASV_65'))


#### Multiple alignment of PA_ASV25 & FL_ASV_65 ####
c(panama_sequences['PA_ASV25'], florida_sequences['FL_ASV_65']) %>%
  as.character() %>%
  as_dna() %>%
  as_DNAbin() %>%
  dist(method = 'maximum')


c(panama_sequences['PA_ASV25'], florida_sequences['FL_ASV_65']) %>%
  as.character() %>%
  enframe(value = 'sequence') %>%
  rowwise(name) %>%
  reframe(str_split(sequence, '') %>%
           unlist %>%
           tibble(sequence = .) %>%
           mutate(position = row_number())) %>%
  pivot_wider(names_from = name,
              values_from = sequence) %>%
  mutate(same_base = PA_ASV25 == FL_ASV_65) %>%
  filter(!same_base)
  # summarise(total_base = n(),
  #           same = sum(same_base))

c(panama_sequences['PA_ASV25'], florida_sequences['FL_ASV_65']) %>%
  Biostrings::writeXStringSet(filepath = '../../intermediate_files/super_interesting_fl_pa.fasta')
