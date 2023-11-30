#### SHAP Family Importance ####
individual_shap_values %>%
  # filter(test_train == 'train') %>%
  # left_join(taxonomy, by = c('asv_id')) %>%
  group_by(wflow_id, asv_id) %>%
  summarise(shap_importance = mean(abs(shap)),
            .groups = 'drop') %>%
  left_join(taxonomy, by = 'asv_id') %>%
  group_by(wflow_id, family) %>%
  summarise(shap_importance = sum(shap_importance),
            n_asv = n_distinct(asv_id),
            .groups = 'drop_last') %>%
  mutate(taxa_rank = rank(-shap_importance)) %>%
  ungroup %>%
  select(wflow_id, family, shap_importance, taxa_rank, n_asv) %>%
  
  filter(taxa_rank <= 5) %>%
  mutate(family = fct_reorder(family, taxa_rank)) %>% 
  rename(model = wflow_id) %>%
  
  ggplot(aes(x = model, y = shap_importance, 
             color = family, 
             group = family)) +
  
  
  
  geom_bump() +
  geom_point(aes(size = n_asv)) +
  
  # geom_text(data = . %>% filter(model == levels(model)[1]),
  #           aes(label = asv_id, x = model),
  #           vjust = -1) +
  
  # coord_cartesian(ylim=c(MAX_ASV, 1)) +
  # scale_y_reverse(labels = 1:MAX_ASV, breaks = 1:MAX_ASV) +
  guides(colour = guide_legend(ncol = 3)) + 
  labs(x = NULL,
       y = 'SHAP Importance Ranking',
       colour = NULL) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')

#### Directional Shap Importance ####
individual_shap_values %>%
  group_by(wflow_id, health, asv_id) %>%
  summarise(shap = mean(shap),
            .groups = 'drop') %>%
  inner_join(filter(model_important_asvs, asv_rank <= MAX_ASV),
             by = c('wflow_id', 'asv_id')) %>%
  mutate(higher_taxonomy = fct_reorder(higher_taxonomy, asv_rank)) %>% 
  rename(model = wflow_id) %>%
  
  ggplot(aes(x = model, y = shap, 
             color = higher_taxonomy, 
             group = interaction(asv_id, health))) +
  
  
  
  geom_bump() +
  geom_point() +
  
  # geom_text(data = . %>% filter(model == levels(model)[1]),
  #           aes(label = asv_id, x = model),
  #           vjust = -1) +
  
  # coord_cartesian(ylim=c(MAX_ASV, 1)) +
  # scale_y_reverse(labels = 1:MAX_ASV, breaks = 1:MAX_ASV) +
  guides(colour = guide_legend(ncol = 3)) + 
  labs(x = NULL,
       y = 'SHAP Importance Ranking',
       colour = NULL) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')
