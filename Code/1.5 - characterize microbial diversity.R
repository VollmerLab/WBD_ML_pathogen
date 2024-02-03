library(tidyverse)
library(microshades)
library(phyloseq)
library(speedyseq)
library(microshades)
library(cowplot)

#### Read in data ####
# updated_taxonomy <- read_csv('../intermediate_files/update_taxonomy.csv', show_col_types = FALSE)
  
asv_data <- bind_rows(read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                               show_col_types = FALSE),
                      
                      read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                               show_col_types = FALSE)) 

unique_asvs <- asv_data %>%
  pull(asv_id) %>%
  unique()


#### Plot Microbial Diversity ####
microbiome_data <- read_rds('../intermediate_files/prepped_microbiome.rds.gz') %>%
  prune_samples(sample_data(.) %>% rownames %in% unique(asv_data$sample_id), .) 


# tax_table(microbiome_data) <- column_to_rownames(updated_taxonomy, 'asv_id') %>% as.matrix

# data <- microbiome_data
# grouping_levels <- c('dataset', 'year', 'season', 'health')
# top_choices <- rev(c('Verrucomicrobiales', 'Alteromonadales', 'Vibrionales', 'Campylobacterales', 'Oceanospirillales'))

make_microbe_plot <- function(data, grouping_levels, highest_taxon = 'order', lower_taxon = 'genus',
                              top_choices = NULL, ...){
  update_data <- data %>%
    tax_glom(lower_taxon) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    as_tibble() %>% 
    mutate(health = str_replace_all(health, c('H' = 'Healthy', 'D' = 'Diseased')),
           # health = factor(health, levels = c('Healthy', 'Diseased')),
           season = str_replace_all(season, c('S' = 'July', 'W' = 'Jan'))) %>%
    group_by(!!!syms(grouping_levels)) %>%
    mutate(total = sum(Abundance)) %>%
    ungroup() %>%
    group_by(!!!syms(grouping_levels), !!sym(lower_taxon)) %>%
    reframe(!!!syms(grouping_levels), 
            !!sym(highest_taxon), !!sym(lower_taxon), 
            total, rel_abun = sum(Abundance)/total) %>%
    distinct() %>%
    rename(Abundance = rel_abun) %>%
    mutate(Sample = str_c(!!!syms(grouping_levels), sep = '_'),
           Sample = if_else(dataset == 'tank', health, str_remove(Sample, 'field_')),
           Sample = fct_relevel(Sample, 'Diseased', after = Inf),
           health = if_else(dataset == 'tank', 'Tank', health),
           health = factor(health, levels = c('Healthy', 'Diseased', 'Tank'))) %>%
    as.data.frame()
  
  if(is.null(top_choices)){
    top_choices <- group_by(update_data, !!sym(highest_taxon)) %>%
      summarise(count = sum(Abundance),
                n_low = n_distinct(!!sym(lower_taxon))) %>%
      arrange(-count) %>%
      filter(n_low > 3) %>%
      slice(1:5) %>%
      pull(!!sym(highest_taxon))
  }
  
  color_objs_ordergenus <- create_color_dfs(update_data, group_level = highest_taxon, 
                                            subgroup_level = lower_taxon,
                                            selected_groups = top_choices, cvd = TRUE)
  mdf_ordergenus <- color_objs_ordergenus$mdf
  cdf_ordergenus <- color_objs_ordergenus$cdf
  
  legend_ordergenus <- custom_legend(mdf_ordergenus, cdf_ordergenus,
                                     group_level = highest_taxon,
                                     subgroup_level = lower_taxon, 
                                     ...)
  
  plot_ordergenus_prelim <- plot_microshades(mdf_ordergenus, cdf_ordergenus) + 
    scale_x_discrete(labels = ~str_remove_all(., '_(Healthy|Diseased)') %>% 
                       str_replace_all('_', '\n') %>% 
                       str_to_title() %>%
                       str_replace_all(c('S' = 'July', 'W' = 'January'))) +
    scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
    facet_grid(~health, scales = "free", space = "free") +
    theme_bw() +
    theme(legend.position = "none", plot.margin = margin(6,20,6,6),
          axis.title.x = element_blank(),
          axis.text = element_text(colour = 'black', size = 12),
          strip.background = element_blank(),
          strip.text = element_text(colour = 'black', size = 16),
          axis.title.y = element_text(colour = 'black', size = 16),
          plot.background = element_rect(colour = 'white'))
  
  plot_out <- plot_grid(plot_ordergenus_prelim, legend_ordergenus,  rel_widths = c(1, .25))
  
  list(plot = plot_out,
       legend = legend_ordergenus,
       color_palette = cdf_ordergenus,
       asv_clumping = as_tibble(mdf_ordergenus) %>%
         select(genus, order, Top_order:group) %>%
         distinct)
  
  
}

make_microbe_plot(microbiome_data,
                  grouping_levels = c('dataset', 'year', 'season', 'health'),
                  highest_taxon = 'order', lower_taxon = 'genus',
                  legend_key_size = 1, legend_text_size = 14) 

microbial_diversity <- microbiome_data %>%
  subset_taxa(rownames(tax_table(microbiome_data)) %in% unique_asvs) %>%
  make_microbe_plot(grouping_levels = c('dataset', 'year', 'season', 'health'),
                    highest_taxon = 'order', lower_taxon = 'genus',
                    # top_choices = rev(c('Verrucomicrobiales', 'Alteromonadales', 
                    #                     'Vibrionales', 'Campylobacterales', 'Oceanospirillales')),
                    legend_key_size = 1, legend_text_size = 14) 
microbial_diversity$plot
ggsave('../Results/Fig2_microbial_diversity.png', height = 10, width = 10, bg = 'white')

write_rds(microbial_diversity[-1], '../intermediate_files/asv_colors.rds')

make_microbe_plot(microbiome_data,
                  grouping_levels = c('dataset', 'year', 'season', 'health'),
                  highest_taxon = 'class', lower_taxon = 'family',
                  legend_key_size = 1, legend_text_size = 14) 




make_microbe_plot(microbiome_data,
                  grouping_levels = c('dataset', 'year', 'season', 'site', 'health'),
                  highest_taxon = 'order', lower_taxon = 'genus')


make_microbe_plot(microbiome_data,
                  grouping_levels = c('dataset', 'health'),
                  highest_taxon = 'order', lower_taxon = 'genus')


#### Summary Stats ####
asv_data %>%
  select(domain:genus, species, asv_id) %>%
  distinct %>%
  pivot_longer(cols = everything()) %>%
  filter(!is.na(value)) %>%
  group_by(name) %>%
  summarise(n_distinct(value))

#### Plot ####
asv_plot_data <- asv_data %>%
  group_by(dataset, season, year, health, sample_id, class, order) %>%
  summarise(log2_cpm_norm = sum(log2_cpm_norm),
            .groups = 'drop') %>%
  group_by(sample_id) %>%
  mutate(pct = log2_cpm_norm / sum(log2_cpm_norm)) %>%
  ungroup %>%
  mutate(combination_variable = paste(class, order, sep = "-"))

n_distinct(asv_plot_data$class)

hex_values <-c(microshades_palette("micro_green", 3, lightest = FALSE), 
               microshades_palette("micro_blue", 3, lightest = FALSE), 
               microshades_palette("micro_purple", 3, lightest = FALSE))

asv_plot_data %>%
  ggplot(aes(x = sample_id, y = pct)) +
  geom_col(aes(fill = combination_variable)) +
  scale_color_manual(values = hex_values, na.translate = FALSE) +
  labs(x = NULL,
       y = "Abundance (%)",
       color = "Class - Order") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white", color = NA),
        plot.title.position = "plot",
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.caption.position = "plot")
