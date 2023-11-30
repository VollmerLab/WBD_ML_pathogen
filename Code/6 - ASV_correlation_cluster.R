##TODO - better plotting of clusters (tsne/umap/nmds/pca etc...)

#### Libraries ####
library(tidyverse)
library(magrittr)
library(umap)
library(Rtsne)
library(patchwork)
library(multidplyr)

library(igraph)
library(tidygraph)
library(ggraph)

cluster <- new_cluster(parallel::detectCores() - 1)

alpha <- 0.05

#### Data ####
taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz', show_col_types = FALSE)

asv_coefs <- read_csv('../Results/asv_coefficients.csv.gz', show_col_types = FALSE)

asv_ranks <- read_csv('../Results/asv_importance.csv.gz', show_col_types = FALSE) %>%
  # filter(p_adjust < alpha) %>%
  left_join(taxonomy, by = 'asv_id')

ora_gsea <- read_csv('../Results/ASV_oraGSEA.csv.gz', show_col_types = FALSE) %>%
  filter(ora_p.value < alpha | gsea_p.value < alpha)

#### Cluster based on mean expressions ####
# asv_dat <- select(asv_coefs, asv_id, starts_with('fc'))
# asv_dat <- select(asv_coefs, asv_id, starts_with('mu'))

cluster_asvs <- function(asv_dat, boot = 50){
  data <- column_to_rownames(asv_dat, 'asv_id')
  
  all_clust <- map(1:max_clust, ~kmeans(data, centers = .x, iter.max = 20, nstart = 20))
  
  boot_clust <- data %>% as.matrix %>% 
    cluster::clusGap(., kmeans, K.max = max_clust, B = boot, verbose = TRUE)
  
  cluster_output_metrics <-  as_tibble(boot_clust$Tab) %>%
    mutate(n_clust = row_number(),
           .before = 1) %>%
    full_join(map_dbl(all_clust, ~.x$tot.withinss) %>%
                tibble(wss = .) %>%
                mutate(n_clust = row_number()),
              by = 'n_clust') %>%
    
    full_join(map_dfr(all_clust[-1], ~cluster::silhouette(.x$cluster, vegdist(data, method = dist_metric)) %>%
                      as_tibble() %>%
                      summarise(mean_sil = mean(sil_width),
                                se_sil = sd(sil_width) / sqrt(n()))) %>%
              mutate(n_clust = row_number() + 1) ,
            by = 'n_clust')
  
  # asv_umap <- umap(data, n_components = 2, random_state = 15)
  
  asv_tsne <- Rtsne(data)
  
  out_df <- map_dfr(1:max_clust, ~enframe(all_clust[[.x]]$cluster, 
                                       name = 'asv_id', value = 'cluster'),
          .id = 'n_clust') %>%
    mutate(n_clust = str_c('clust_', n_clust)) %>%
    pivot_wider(names_from = n_clust, values_from = cluster) %>%
    left_join(set_colnames(asv_tsne$Y, c('tSNE1', 'tSNE2')) %>% 
                as_tibble() %>%
                mutate(asv_id = rownames(data), .before = 1), 
              .,
              by = 'asv_id')
  
  list(out_df = out_df, cluster_metrics = cluster_output_metrics)
}


mean_cluster_groups <- select(asv_coefs, asv_id, starts_with('mu')) %>%
  cluster_asvs(boot = 100)


cluster_stats_plot <- mean_cluster_groups$cluster_metrics %>%
  rename(mean_wss = wss,
         mean_gap = gap,
         se_gap = SE.sim) %>%
  select(n_clust, starts_with('mean'), starts_with('se')) %>%
  pivot_longer(cols = c(starts_with('mean'), starts_with('se')),
               names_to = c('.value', 'metric'),
               names_pattern = '(.*)_(.*)') %>%
  filter(!is.na(mean)) %>%
  ggplot(aes(x = n_clust, y = mean, ymin = mean - se, ymax = mean + se)) +
  geom_line() +
  geom_pointrange() +
  facet_wrap(~metric, scales = 'free_y', nrow = 1) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank())

mean_cluster_tsne <- mean_cluster_groups$out_df %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(cluster = str_c('c', clust_3)) %>%
  ggplot(aes(x = tSNE1, y = tSNE2)) +
  geom_text(aes(label = asv_id,
                colour = cluster), show.legend = FALSE) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank())

mean_cluster_tsne / cluster_stats_plot

mean_cluster_groups$out_df %>%
  filter(clust_4 == clust_4[asv_id == 'ASV25']) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  select(-starts_with('clust'))

str_c(order, family, genus, asv_id, sep = '; ')


#### Cluster based on FC ####
fc_cluster_groups <- select(asv_coefs, asv_id, starts_with('fc')) %>%
  cluster_asvs(boot = 100)


fc_cluster_stats_plot <- fc_cluster_groups$cluster_metrics %>%
  rename(mean_wss = wss,
         mean_gap = gap,
         se_gap = SE.sim) %>%
  select(n_clust, starts_with('mean'), starts_with('se')) %>%
  pivot_longer(cols = c(starts_with('mean'), starts_with('se')),
               names_to = c('.value', 'metric'),
               names_pattern = '(.*)_(.*)') %>%
  filter(!is.na(mean)) %>%
  ggplot(aes(x = n_clust, y = mean, ymin = mean - se, ymax = mean + se)) +
  geom_line() +
  geom_pointrange() +
  facet_wrap(~metric, scales = 'free_y', nrow = 1) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank())

fc_cluster_tsne <- fc_cluster_groups$out_df %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(cluster = str_c('c', clust_4)) %>%
  ggplot(aes(x = tSNE1, y = tSNE2)) +
  geom_text(aes(label = asv_id,
                colour = cluster), show.legend = FALSE) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank())

fc_cluster_tsne / fc_cluster_stats_plot

fc_cluster_groups$out_df %>%
  filter(clust_4 == clust_4[asv_id == 'ASV25']) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  select(-starts_with('clust'))

#### Correlation based on Mean Expression ####
plot_corr_heatmap <- function(corr_data, alpha = 1){
  the_tree <- select(corr_data, asv_id1, asv_id2, estimate) %>%
    arrange(asv_id1, asv_id2) %>%
    pivot_wider(names_from = asv_id2,
                values_from = estimate) %>%
    column_to_rownames('asv_id1') %>%
    as.matrix() %>%
    scale %>%
    dist() %>%
    hclust
  
  corr_data %>%
    filter(p.value < alpha) %>%
    mutate(across(starts_with('asv_id'), 
                  ~factor(., levels = the_tree$labels[the_tree$order]))) %>%
    ggplot(aes(x = asv_id1, y = asv_id2, fill = estimate)) +
    geom_raster() +
    scale_fill_gradient2(limits = c(-1, 1))
}

corr_asv_mean <- select(asv_coefs, asv_id, starts_with('mu')) %>%
  pivot_longer(cols = -asv_id, 
               names_to = 'timepoint',
               values_to = 'log2cpm') %>%
  nest(data = -c(asv_id)) %>%
  expand(nesting(asv_id1 = asv_id, data1 = data), 
         nesting(asv_id2 = asv_id, data2 = data)) %>%
  rowwise(asv_id1, asv_id2) %>%
  partition(cluster) %>%
  summarise(broom::tidy(cor.test(data1$log2cpm, data2$log2cpm))) %>%
  collect() %>%
  ungroup()
  
plot_corr_heatmap(corr_asv_mean, alpha = 0.05)

as_tbl_graph(corr_asv_mean, directed = TRUE) %>%
  activate(edges) %>%
  filter(from != to) %>%
  filter(from > to) %>%
  mutate(p_adjust = p.adjust(p.value, method = 'fdr')) %>%
  filter(p_adjust < 0.05,
         abs(estimate) > 0.9) %>%
  ggraph(layout = 'nicely') +
  geom_edge_link(aes(colour = estimate > 0)) +
  geom_node_point(aes(colour = name == 'ASV25')) +
  scale_colour_manual(values = c('TRUE' = 'green', 'FALSE' = 'black'))

as_tbl_graph(corr_asv_mean, directed = TRUE) %>%
  activate(edges) %>%
  filter(from != to) %>%
  filter(from > to) %>%
  mutate(p_adjust = p.adjust(p.value, method = 'fdr')) %>%
  filter(p_adjust < 0.05,
         abs(estimate) > 0.9) %>%
  activate(nodes) %>%
  mutate(cluster = group_edge_betweenness()) %>%
  ggraph(layout = 'nicely') +
  geom_edge_link() +
  geom_node_point(aes(colour = as.character(cluster),
                      size = name == 'ASV25'),
                  show.legend = FALSE) +
  scale_size_manual(values = c('TRUE' = 3, 'FALSE' = 1))


as_tbl_graph(corr_asv_mean, directed = TRUE) %>%
  activate(edges) %>%
  filter(from != to) %>%
  filter(from > to) %>%
  mutate(p_adjust = p.adjust(p.value, method = 'fdr')) %>%
  filter(p_adjust < 0.05,
         abs(estimate) > 0.9) %>%
  activate(nodes) %>%
  mutate(cluster = group_edge_betweenness()) %>%
  filter(cluster == cluster[name == 'ASV25']) %>%
  as_tibble() %>%
  select(name) %>%
  rename(asv_id = name) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  count(across(class:family))



as_tbl_graph(corr_asv_mean, directed = TRUE) %>%
  activate(edges) %>%
  filter(from != to) %>%
  filter(from > to) %>%
  mutate(p_adjust = p.adjust(p.value, method = 'fdr')) %>%
  filter(p_adjust < 0.05,
         abs(estimate) > 0.9) %>%
  activate(nodes) %>%
  left_join(taxonomy, by = c('name' = 'asv_id')) %>%
  ggraph(layout = 'nicely') +
  geom_edge_link() +
  geom_node_point(aes(colour = family,
                      size = name == 'ASV25'),
                  show.legend = FALSE) +
  scale_size_manual(values = c('TRUE' = 3, 'FALSE' = 1))

#### Correlation based on FC ####
corr_asv_fc <- select(asv_coefs, asv_id, starts_with('fc')) %>%
  pivot_longer(cols = -asv_id, 
               names_to = 'timepoint',
               values_to = 'log2cpm') %>%
  nest(data = -c(asv_id)) %>%
  expand(nesting(asv_id1 = asv_id, data1 = data), 
         nesting(asv_id2 = asv_id, data2 = data)) %>%
  rowwise(asv_id1, asv_id2) %>%
  partition(cluster) %>%
  summarise(broom::tidy(cor.test(data1$log2cpm, data2$log2cpm))) %>%
  collect() %>%
  ungroup()

plot_corr_heatmap(corr_asv_fc, alpha = 1)


as_tbl_graph(corr_asv_fc, directed = TRUE) %>%
  activate(edges) %>%
  filter(from != to) %>%
  filter(from > to) %>%
  mutate(p_adjust = p.adjust(p.value, method = 'fdr')) %>%
  filter(p_adjust < 0.1) %>%
  ggraph(layout = 'nicely') +
  geom_edge_link(aes(colour = estimate > 0)) +
  geom_node_point(aes(colour = name == 'ASV25'))

#### Cluster ASVs ####

ora_gsea
inner_join(asv_ranks, 
           ora_gsea,
           by = 'family')
