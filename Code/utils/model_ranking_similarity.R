library(tidyverse)
library(magrittr)
library(gespeR)
library(ggdendro)

all_asv_ranks <- read_csv('../../intermediate_files/model_shaps.csv.gz', 
                      show_col_types = FALSE) %>%
  group_by(wflow_id, train_test) %>%
  mutate(across(ends_with('value'), ~scale(.)[,1])) %>%
  ungroup %>%
  pivot_longer(cols = starts_with('ASV'),
               names_to = c('asv_id', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  group_by(wflow_id, asv_id) %>%
  summarise(shap_importance = mean(abs(shap)),
            .groups = 'drop') %>%
  pivot_wider(names_from = wflow_id, 
              values_from = shap_importance) %>%
  mutate(across(where(is.numeric), ~rank(-.)))


asv_ranks <- read_csv('../../Results/asv_importance.csv.gz',
                      show_col_types = FALSE) %>%
  rename(base_median = median_rank,
         base_model = response) %>%
  filter(!is.na(base_model)) %>%
  left_join(select(all_asv_ranks, asv_id, base_knn, base_pls, base_null),
            by = 'asv_id') %>%
  select(asv_id, starts_with('base')) %>%
  arrange(base_model)


rbo_permute <- function(r1, r2, ...){
  names(r1) <- sample(names(r1))
  # names(r2) <- sample(names(r2))
  rbo(r1, r2, ...)
}

my_tau <- function(list1, list2, perm = FALSE){
  if(perm){
    names(list1) <- sample(names(list1))
  }
  
  A <- fct_inorder(names(list1))
  B <- factor(names(list2), levels = levels(A))
  
  out <- Kendall::Kendall(A, B)
  tibble(tau = out$tau, tau_p = out$sl)
}


ranking_similarities <- asv_ranks %>%
  pivot_longer(cols = starts_with('base'),
               names_to = 'model',
               values_to = 'rank',
               names_prefix = 'base_') %>%
  # filter(rank <= 5) %>%
  nest_by(model) %>%
  mutate(rank_order = list(data %$%
                             set_names(rank, asv_id) %>%
                             sort())) %>%
  ungroup %>%
  expand(nesting(model.1 = model, rank.1 = rank_order), 
         nesting(model.2 = model, rank.2 = rank_order)) %>%
  rowwise(model.1, model.2) %>%
  mutate(#broom::tidy(cor.test(rank.1, rank.2, method = 'kendall')),
         my_tau(rank.1, rank.2, perm = FALSE),
         # tau_p = mean(c(tau, replicate(10 - 1, my_tau(rank.1, rank.2, perm = TRUE)))),
         
         # across(starts_with('rank'), ~list(.[. < max(.)])),
         
         rbo = rbo(rank.1, rank.2, p = 150/min(length(rank.1), length(rank.2)), 
                   k = min(100, length(rank.1), length(rank.2)), side = 'bottom', 
                   uneven.lengths = FALSE),
         
         rbo_p = mean(c(rbo, replicate(1000 - 1, rbo_permute(rank.1, rank.2, p = 150/min(length(rank.1), length(rank.2)), 
                                                        k = min(100, length(rank.1), length(rank.2)), side = 'bottom', 
                                                        uneven.lengths = FALSE))) >= rbo)) %>%
  
  # dplyr::select(-method, -alternative) %>%
  ungroup

#### Model Similarities ####
ranking_similarities %>%
  select(model.1, model.2, rbo, rbo_p, tau, tau_p) %>%
  
  filter(!model.1 %in% c('knn', 'pls', 'median'),
         !model.2 %in% c('knn', 'pls', 'median')) %>%
  
  filter(model.1 > model.2) %>%
  filter(!model.1 %in% c('model'),
         !model.2 %in% c('model')) %>%
  filter(model.1 != 'null',
         model.2 != 'null') %>%
  summarise(across(c(rbo, tau), c(~mean(.), ~sd(.) / sqrt(n()))))

ranking_similarities %>%
  select(model.1, model.2, rbo, rbo_p, tau, tau_p) %>%
  
  filter(!model.1 %in% c('knn', 'pls', 'median'),
         !model.2 %in% c('knn', 'pls', 'median')) %>%
  
  filter(model.1 > model.2) %>%
  filter(model.1 %in% c('model', 'null')) %>%
  group_by(model.1) %>%
  summarise(across(c(rbo, tau), c(~mean(.), ~sd(.) / sqrt(n()))))


test_rankings <- ranking_similarities %>%
  select(model.1, model.2, rbo, rbo_p, tau, tau_p) %>%
  
  filter(!model.1 %in% c('knn', 'pls', 'median'),
         !model.2 %in% c('knn', 'pls', 'median'),
         model.1 != model.2) %>%
  
  filter(model.1 %in% c('null', 'model') | 
           (model.1 > model.2 & 
              !model.1 %in% c('null', 'model') &
              !model.2 %in% c('null', 'model'))) %>%
  filter(!(model.1 == 'model' & model.2 == 'null' | model.1 == 'null' & model.2 == 'model')) %>%
  mutate(grouping = case_when(model.1 %in% c('null', 'model') ~ model.1,
                              TRUE ~ 'ML')) %>%
  mutate(tau = if_else(tau < 0, 1/100, tau))

library(betareg)
library(emmeans)
rbo_beta <- betareg(rbo ~ grouping, data = test_rankings)
car::Anova(rbo_beta)

emmeans(rbo_beta, ~grouping) %T>%
  print %>%
  contrast('pairwise')

tau_beta <- betareg(tau ~ grouping, data = test_rankings)
car::Anova(tau_beta)

emmeans(tau_beta, ~grouping) %T>%
  print %>%
  contrast('pairwise')

#### Plot clustering ####
ranking_similarities %>%
  select(model.1, model.2, rbo) %>%
  filter(!model.1 %in% c('model', 'median'),
         !model.2 %in% c('model', 'median')) %>%
  
  filter(!model.1 %in% c('knn', 'pls'),
         !model.2 %in% c('knn', 'pls')) %>%
  mutate(rbo = 1 - rbo) %>%
  pivot_wider(names_from = model.2,
              values_from = rbo) %>%
  column_to_rownames('model.1') %>%
  as.matrix %>%
  as.dist() %>%
  hclust()  %>%
  ggdendrogram(rotate = TRUE, size = 2)
ggsave('../../Results/model_similarity.png', height = 5, width = 5)

ranking_similarities %>%
  select(model.1, model.2, tau) %>%
  filter(!model.1 %in% c('model', 'median'),
         !model.2 %in% c('model', 'median')) %>%
  filter(!model.1 %in% c('knn', 'pls'),
         !model.2 %in% c('knn', 'pls')) %>%
  mutate(tau = 1 - tau) %>%
  pivot_wider(names_from = model.2,
              values_from = tau) %>%
  column_to_rownames('model.1') %>%
  as.matrix %>%
  as.dist() %>%
  hclust()  %>%
  ggdendrogram(rotate = TRUE, size = 2)


list1 <- tmp2$rank.1[[1]]
list2 <- tmp2$rank.2[[1]]



list1 <- tmp2$rank.1[[3]][1:10]
list2 <- tmp2$rank.2[[3]][1:10]

my_tau(tmp2$rank.1[[3]][1:10], tmp2$rank.2[[3]][1:10])

tmp2$rank.2[[7]]


tmp2 %>%
  ggplot(aes(x = model.1, y = model.2, fill = rbo)) +
  geom_raster()


tmp2 %>%
  select(model.1, model.2, estimate) %>%
  mutate(estimate = 1 - estimate) %>%
  pivot_wider(names_from = model.2,
              values_from = estimate) %>%
  column_to_rownames('model.1') %>%
  as.matrix %>%
  as.dist() %>%
  hclust()  %>%
  ggdendrogram(rotate = TRUE, size = 2)

