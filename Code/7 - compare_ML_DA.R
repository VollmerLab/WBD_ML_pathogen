
library(multcomp)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(ggVennDiagram)
library(ComplexUpset)
library(cowplot)


#### Data ####
taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz',
                     show_col_types = FALSE) %>%
  mutate(across(where(is.character), 
                ~str_replace_na(., replacement = '')))

ml_ranking <- read_csv('../Results/asv_importance.csv.gz', show_col_types = FALSE)

full_data <- list.files(path = '../intermediate_files', 
                        pattern = 'asv_groups.csv', full.names = TRUE) %>%
  tibble(files = .) %>%
  rowwise(files) %>%
  reframe(read_csv(files, show_col_types = FALSE)) %>%
  mutate(analysis = str_extract(files, 'DA|ML'),
         # analysis = case_when(analysis == 'DA' ~ 'Differential Abundance',
         #                      analysis == 'ML' ~ 'Machine Learning'),
         .keep = 'unused') %>%
  mutate(important = `Disease Associated` | `Healthy Associated`,
         .keep = 'all')

#### Rate of Change ####
consistency_change <- full_data %>%
  summarise(important = sum(important),
            consistent = sum(consistent),
            .by = c(`Disease Associated`, analysis)) %>%
  mutate(association = if_else(`Disease Associated`, 'd', 'h'),
         inconsistent = important - consistent)

change_model <- glm(cbind(consistent, important - consistent) ~ analysis * association, 
      family = binomial(),
      data = consistency_change)
car::Anova(change_model)

emmeans(change_model, ~analysis, type = 'response') %>% contrast('revpairwise')
emmeans(change_model, ~association, type = 'response') %>% 
  contrast('pairwise') 

emmeans(change_model, ~ analysis * association, type = 'response') %>%
  broom::tidy() %>%
  ggplot(aes(x = analysis, y = prob, 
             ymin = prob - std.error, 
             ymax = prob + std.error,
             colour = association)) +
  geom_errorbar(width = 0.1, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5))

#### Important & Consistent ####
transformed_data <- full_data %>%
  select(-contains('Associated')) %>%
  pivot_wider(names_from = analysis, 
              values_from = c(important, consistent),
              values_fill = FALSE) %>%
  mutate(importance_group = case_when(important_DA & important_ML ~ 'Both',
                                      important_ML ~ 'ML',
                                      important_DA ~ 'DA',
                                      TRUE ~ 'Neither'),
         consistent_group = case_when(consistent_DA & consistent_ML ~ 'Both',
                                      consistent_ML ~ 'ML',
                                      consistent_DA ~ 'DA',
                                      TRUE ~ 'Neither'),
         .keep = 'all') %>%
  mutate(joint_group = str_c('Important: ', importance_group,
                             '\nConsistent: ', consistent_group)) %>%
  left_join(select(ml_ranking, asv_id, starts_with('base')), 
            by = 'asv_id') %>%
  pivot_longer(cols = starts_with('base'),
               names_prefix = 'base_',
               names_to = 'ML_model',
               values_to = 'rank') %>%
  
  mutate(n_asv = n_distinct(asv_id),
         .by = joint_group)

transformed_data %>%
  filter(n_asv == 1)


rank_model <- lmer(log(rank) ~ ML_model * (importance_group / consistent_group) + 
                     (1 | asv_id),
      data = filter(transformed_data, n_asv > 5))

anova(rank_model, ddf = 'Kenward-Roger')

emmeans(rank_model, ~importance_group + consistent_group, type = 'response') %>%
  regrid() %>%
  contrast('trt.vs.ctrl', reverse = FALSE)

emmeans(rank_model, ~importance_group + consistent_group, type = 'response') %>%
  regrid() %>%
  contrast('consec', reverse = FALSE)

emmeans(rank_model, ~ML_model | importance_group + consistent_group, type = 'response') %>%
  regrid() %>%
  contrast('pairwise')

emmeans(rank_model, ~ML_model | importance_group + consistent_group, type = 'response',
        at = list(importance_group = 'DA', consistent_group = 'DA')) %>%
  regrid() %>%
  contrast(method = list('lasso - others' = c(-1/3, 1, -1/3, -1/3)))

emmeans(rank_model, ~ML_model * (importance_group + consistent_group), type = 'response',
        at = list(importance_group = c('DA', 'Neither'), consistent_group = c('DA', 'Neither')), 
        nesting = NULL) %>%
  regrid() %>%
  contrast(method = list('DA.DA: lasso - others' = c(-1/3, 1, -1/3, -1/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         'DA.DA.lasso - N.N.all' = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/4, -1/4, -1/4, -1/4)))


average_model_rank_df <- emmeans(rank_model, ~(importance_group + consistent_group), type = 'response') %>%
  cld(Letters = LETTERS) %>%
  broom::tidy(conf.int = TRUE) %>%
  inner_join(distinct(transformed_data, importance_group, consistent_group, joint_group),
             by = c('importance_group', 'consistent_group')) %>%
  mutate(.group = str_trim(.group))


rank_plot <- emmeans(rank_model, ~ML_model * (importance_group + consistent_group), type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  right_join(distinct(transformed_data, importance_group, consistent_group, joint_group),
             by = c('importance_group', 'consistent_group')) %>%
  # mutate(joint_group = factor(joint_group) %>%
  #          fct_relevel('Important: Both\nConsistent: ML', 'Important: ML\nConsistent: Neither', 
  #                      after = Inf)) %>%
  
  mutate(joint_group = factor(joint_group) %>%
           fct_relevel('Important: ML\nConsistent: Neither',
                       after = 3)) %>%
  
  ggplot(aes(x = joint_group, y = response, colour = ML_model)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                position = position_dodge(0.5),
                width = 0.1, show.legend = FALSE) +
  geom_point(position = position_dodge(0.5)) +
  
  geom_point(data = filter(transformed_data, n_asv < 5),
             inherit.aes = FALSE,
             aes(x = joint_group, y = rank, colour = ML_model,
                 group = asv_id),
             position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.height = 0,
                                             jitter.width = 0.2)) +
  geom_point(data = filter(transformed_data, n_asv < 5) %>%
               summarise(rank = mean(rank), 
                         .by = c(asv_id, joint_group)),
             inherit.aes = FALSE, size = 2.5,
             aes(x = joint_group, y = rank, group = asv_id),
             position = position_dodge(0.5)) +
  
  
  geom_errorbar(data = average_model_rank_df, colour = 'black',
                aes(ymin = conf.low,
                    ymax = conf.high), width = 0.25, linewidth = 1) +
  geom_point(data = average_model_rank_df, colour = 'black', size = 2.5) +
  geom_text(data = average_model_rank_df, colour = 'black', 
            aes(label = .group, y = Inf), vjust = 1.5) +
  
  scale_color_discrete(na.translate = F, labels = c('forest' = 'RF', 'lasso' = 'Lasso',
                                                    'mlp' = 'MLP', 'svmLinear' = 'SVM')) +
  # guides(colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(size = 5))) +
  guides(colour = guide_legend(position = 'inside', override.aes = list(size = 5))) +
  labs(tag = 'B',
       y = 'Subcomponent\nRank',
       colour = 'ML Sub-model') +
  theme_classic() +
  theme(legend.position.inside = c(0.25, 0.75),
        legend.key = element_blank(),
        panel.background = element_rect(colour = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        plot.tag = element_text(colour = 'black', size = 18, 
                                face = 'bold', hjust = 0))

rank_plot

#### Extra side panel ####
set_size_plot <- full_data %>%
  pivot_longer(cols = c(consistent, important),
               names_to = 'type',
               values_to = 'keep') %>%
  filter(keep) %>%
  mutate(group = str_c(analysis, type, sep = ': ')) %>%
  summarise(n_asv = n_distinct(asv_id),
            n_healthy = n_distinct(asv_id[`Healthy Associated`]),
            n_disease = n_distinct(asv_id[`Disease Associated`]),
            .by = 'group') %>%
  select(-n_asv) %>%
  pivot_longer(cols = c(n_healthy, n_disease),
               names_to = 'type', values_to = 'n') %>%
  mutate(print_value = if_else(type == 'n_disease', sum(n), n),
         .by = group)  %>%
  mutate(group = factor(group, levels = c('DA: consistent', 'DA: important',
                                          'ML: consistent', 'ML: important'))) %>%
  ggplot(aes(y = group, x = n)) +
  # geom_col(data = . %>% summarise(n = sum(n), .by = 'group')) +
  geom_col(aes(fill = type),
           position = position_dodge(0.5), width = 0.45) +
  geom_text(position = position_dodge(0.5),
            aes(label = n, colour = type), 
            hjust = 0, size = 3, fontface = 'bold') +

  scale_fill_manual(values = c('n_healthy' = '#3B9AB2',
                               'n_disease' = '#F21A00')) +
  scale_colour_manual(values = c('n_healthy' = '#3B9AB2',
                               'n_disease' = '#F21A00')) +
  labs(tag = 'E',
       x = 'Set Size') +
  theme_classic() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black'),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.x = element_text(angle = 0,  colour = 'black', size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.tag = element_text(colour = 'black', size = 18, 
                                face = 'bold', hjust = 0))
set_size_plot

#### Correlation Analysis ####
top_asvs <- filter(transformed_data, if_all(where(is.logical), ~.)) %>%
  select(asv_id) %>%
  distinct %>%
  pull(asv_id)

correlation_data <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                             show_col_types = FALSE) %>%
  select(asv_id, sample_id, log2_cpm_norm, health) %>%
  full_join(., .,
            by = c('sample_id', 'health'),
            relationship = "many-to-many") %>%
  filter(health == 'D') %>%
  group_by(asv_id.x, asv_id.y, health) %>%
  reframe(r = cor(log2_cpm_norm.x, log2_cpm_norm.y)) %>%
  filter(asv_id.x != asv_id.y) %>%
  filter(asv_id.y %in% top_asvs) %>%
  inner_join(distinct(transformed_data, asv_id, 
                      importance_group, consistent_group,
                      joint_group, n_asv),
             by = c('asv_id.x' = 'asv_id'))

ggplot(correlation_data, aes(x = r)) +
  geom_histogram() +
  facet_wrap(~importance_group + consistent_group, scales = 'free_y')

correlation_model <- lmer(r ~ importance_group / consistent_group + 
                            (1 | asv_id.x) + (1 | asv_id.y),
                          data = filter(correlation_data, n_asv > 5))

anova(correlation_model, ddf = 'Kenward-Roger')

correlation_means <- emmeans(correlation_model, ~consistent_group, 
                             pbkrtest.limit = 5000, infer = TRUE) 
correlation_means %>%
  contrast('identity')

correlation_means %>%
  contrast('trt.vs.ctrl')

correlation_plot <- correlation_means %>%
  cld(Letters = LETTERS) %>%
  broom::tidy(conf.int = TRUE) %>%
  mutate(.group = str_trim(.group)) %>%
  right_join(distinct(transformed_data, importance_group, 
                      consistent_group, joint_group),
             by = c('importance_group', 'consistent_group')) %>%
  mutate(joint_group = factor(joint_group) %>%
           fct_relevel('Important: ML\nConsistent: Neither',
                       after = 3),
         .group = c('A', 'A', 'A', 'A', 'B', NA_character_, NA_character_)) %>%
  
  ggplot(aes(x = joint_group, y = estimate, #colour = health,
             ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_errorbar(width = 0.1, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5)) +
  
  geom_point(data = filter(correlation_data, n_asv < 5) %>%
               summarise(r = mean(r), 
                         .by = c('asv_id.x', 'joint_group')),
             inherit.aes = FALSE, position = position_dodge(0.5),
             aes(x = joint_group, y = r, group = asv_id.x)) +
  
  geom_text(aes(label = .group, y = Inf), vjust = 1.5) +
  labs(tag = 'A',
       y = 'Correlation in\nDiseased Corals') +
  theme_classic() +
  theme(legend.position = 'top',
        legend.key = element_blank(),
        panel.background = element_rect(colour = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        plot.tag = element_text(colour = 'black', size = 18, 
                                face = 'bold', hjust = 0))


correlation_plot

#
#### Interaction Plot ####
intersection_plot <- distinct(transformed_data, asv_id, joint_group) %>%
  left_join(select(full_data, asv_id, `Disease Associated`, `Healthy Associated`) %>%
              summarise(disease_associated = any(`Disease Associated`),
                        .by = 'asv_id'),
            by = 'asv_id') %>%
  count(joint_group, disease_associated) %>%
  mutate(type = case_when(joint_group == 'Important: Neither\nConsistent: Neither' ~ 'none',
                          disease_associated ~ 'disease',
                          TRUE ~ 'healthy')) %>%
  mutate(joint_group = factor(joint_group) %>%
           fct_relevel('Important: ML\nConsistent: Neither',
                       after = 3)) %>%
  ggplot(aes(x = joint_group, y = n, fill = type)) +
  geom_col(position = position_dodge(0.5), 
           width = 0.45) +
  geom_text(position = position_dodge(0.5),
            aes(label = n, colour = type), 
            vjust = -0.5, size = 3, 
            fontface = 'bold') +
  
  scale_fill_manual(values = c('healthy' = '#3B9AB2',
                               'disease' = '#F21A00',
                               'none' = 'gray65')) +
  scale_colour_manual(values = c('healthy' = '#3B9AB2',
                                 'disease' = '#F21A00',
                                 'none' = 'gray65')) +
  scale_y_continuous(labels=scales::comma_format(), 
                     expand = expansion(mult = c(0.01, 0.05))) +
  labs(tag = 'C',
       y = 'Intersection\nSize') +
  theme_classic() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        plot.tag = element_text(colour = 'black', size = 18, 
                                face = 'bold', hjust = 0))

#### Composition Set Up ####
colour_options <- read_rds('../intermediate_files/asv_colors.rds')
microbe_colors <- set_names(colour_options$color_palette$hex,
                            colour_options$color_palette$group)
levels(colour_options$asv_clumping$Top_order)

composition_data <- transformed_data %>%
  select(-ML_model, -rank) %>%
  distinct %>%
  select(asv_id, 
         `DA: Consistent` = consistent_DA,
         `DA: Important` = important_DA, 
         `ML: Consistent` = consistent_ML, 
         `ML: Important` = important_ML, 
         everything()) %>%
  left_join(taxonomy, by = 'asv_id') %>%
  left_join(mutate(colour_options$asv_clumping,
                   genus = str_remove_all(genus, '\\<i\\>|\\</i\\>')),
            by = c('order', 'genus')) %>%
  mutate(group = case_when(!is.na(group) ~ group,
                           is.na(group) & order == 'Oceanospirillales' ~ 'Oceanospirillales-Other',
                           is.na(group) & order == 'Vibrionales' ~ 'Vibrionales-Other',
                           is.na(group) & order == 'Alteromonadales' ~ 'Alteromonadales-Other',
                           is.na(group) & order == 'Verrucomicrobiales' ~ 'Verrucomicrobiales-Other',
                           is.na(group) & order == 'Campylobacterales' ~ 'Campylobacterales-Other',
                           TRUE ~ 'Other-Other'),
         group = factor(group, levels = levels(colour_options$asv_clumping$group))) %>%
  # filter(!is.na(group)) %>%
  select(-Top_order, -Top_genus) %>%
  rename(the_colour = group) %>%
  left_join(select(full_data, asv_id, `Disease Associated`, `Healthy Associated`) %>%
              summarise(disease_associated = any(`Disease Associated`),
                        .by = 'asv_id'),
            by = 'asv_id') %>%
  mutate(type = case_when(joint_group == 'Important: Neither\nConsistent: Neither' ~ 'none',
                          disease_associated ~ 'disease',
                          TRUE ~ 'healthy'))

composition_data %>%
  summarise(n = n(),
            .by = c(joint_group, the_colour, type)) %>%
  quitte::ggplot_bar_stacked_dodged(aes(x = type, 
                                        y = n,
                                        fill = the_colour,
                                        dodge = joint_group)) +
  scale_fill_manual(values = microbe_colors) +
  labs(y = 'Relative\nAbundance',
       tag = 'C') +
  theme_classic() +
  theme(legend.position = 'none')

composition_data %>%
  ggplot(aes(x = joint_group, fill = the_colour)) +
  geom_bar(fill = 'white', colour = 'black', stat = 'count', 
           position = 'fill', show.legend = FALSE) +
  geom_bar(stat = 'count', position = 'fill', show.legend = FALSE,
           colour = 'black') +
  scale_y_continuous(labels=scales::percent_format(), 
                     expand = expansion(mult = c(0.01, 0.05))) +
  scale_fill_manual(values = microbe_colors) +
  labs(y = 'Relative\nAbundance',
       tag = 'C') +
  theme_classic() 

composition_data %>%
  ggplot(aes(x = as.numeric(interaction(joint_group, type)), 
             fill = the_colour)) +
  geom_bar(fill = 'white', colour = 'black', stat = 'count', 
           position = 'fill', show.legend = FALSE) +
  geom_bar(stat = 'count', position = 'fill', 
           show.legend = FALSE,
           colour = 'black') +
  scale_y_continuous(labels=scales::percent_format(), 
                     expand = expansion(mult = c(0.01, 0.05))) +
  scale_x_continuous(breaks = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5),
                     labels = unique(composition_data$joint_group)) +
  scale_fill_manual(values = microbe_colors) +
  labs(y = 'Relative\nAbundance',
       tag = 'C') +
  theme_classic() 

composition_data %>%
  ggplot(aes(x = type, fill = the_colour)) +
  geom_bar(fill = 'white', colour = 'black', stat = 'count', 
           position = 'fill', show.legend = FALSE) +
  geom_bar(stat = 'count', position = 'fill', show.legend = FALSE,
           colour = 'black') +
  scale_y_continuous(labels=scales::percent_format(), 
                     expand = expansion(mult = c(0.01, 0.05))) +
  scale_fill_manual(values = microbe_colors) +
  facet_grid(~ joint_group, scales = 'free_x', space = 'free_x') +
  labs(y = 'Relative\nAbundance',
       tag = 'C') +
  theme_classic() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black'),
        axis.ticks.y = element_line(colour = 'black'),
        axis.minor.ticks.y.left = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10),
        plot.tag = element_text(colour = 'black', size = 18, 
                                face = 'bold', hjust = 0))
  

#### Full Analysis Plot ####
full_plot <- composition_data %>%
  select(-disease_associated) %>%
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        sort_sets = FALSE, 
        name = NULL,
        
        sort_intersections=FALSE,
        intersections=list(c('ML: Important', 'DA: Important', 'ML: Consistent', 'DA: Consistent'),
                           c('ML: Important', 'DA: Important', 'ML: Consistent'),
                           c('ML: Important', 'DA: Important'),
                           'ML: Important',
                           c('DA: Consistent', 'DA: Important'),
                           'DA: Important',
                           'Outside of known sets'),
        
        base_annotations = list(
          'Intersection size' = intersection_size(
            counts = TRUE
          ) + 
            scale_y_continuous(labels=scales::comma_format(), 
                               expand = expansion(mult = c(0.01, 0.05))) +
            labs(tag = 'C',
                 y = 'Intersection\nSize') +
            theme_classic() +
            theme(legend.position = 'none',
                  panel.background = element_rect(colour = 'black'),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title = element_text(colour = 'black', size = 14),
                  axis.text = element_text(colour = 'black', size = 10),
                  plot.tag = element_text(colour = 'black', size = 18, 
                                          face = 'bold', hjust = 0))
        ),
        
        annotations = list(
          'Correlation' = ggplot(mapping = aes(y = n_asv)) +
            geom_point() +
            labs(tag = 'A') +
            theme_classic(), 
          
          'Rank' = ggplot(mapping = aes(y = n_asv)) +
            geom_point() +
            labs(tag = 'B') +
            theme_classic()
          
          # 'Order' = ggplot(mapping = aes(fill = the_colour)) +
          #   geom_bar(fill = 'white', colour = 'black', stat = 'count', 
          #            position = 'fill', show.legend = FALSE) +
          #   geom_bar(stat = 'count', position = 'fill', show.legend = FALSE) +
          #   scale_y_continuous(labels=scales::percent_format(), 
          #                      expand = expansion(mult = c(0.01, 0.05))) +
          #   scale_fill_manual(values = microbe_colors) +
          #   labs(y = 'Relative\nAbundance',
          #        tag = 'C') +
          #   theme_classic() +
          #   theme(legend.position = 'none',
          #         panel.background = element_rect(colour = 'black'),
          #         axis.ticks.y = element_line(colour = 'black'),
          #         axis.minor.ticks.y.left = element_blank(),
          #         axis.text.x = element_blank(),
          #         axis.ticks.x = element_blank(),
          #         axis.title.x = element_blank(),
          #         axis.title = element_text(colour = 'black', size = 14),
          #         axis.text = element_text(colour = 'black', size = 10),
          #         plot.tag = element_text(colour = 'black', size = 18, 
          #                                 face = 'bold', hjust = 0))
        ),
        
        set_sizes=(
          upset_set_size(geom = geom_bar(fill = c('gray65', '#F21A00', '#3B9AB2', 'green')),
                         position = 'right') + #'gray65'
            geom_text(aes(label = after_stat(count)), hjust = 1.1, stat = 'count', 
                      colour = 'white', fontface = 'bold', size = 4) +
            # + annotate(geom='text', label='@', x='Drama', y=850, color='white', size=3) +
            # expand_limits(y = 200) +
            labs(tag = 'E') +
            theme_classic() +
            theme(panel.background = element_rect(colour = 'black'),
                  axis.ticks.x = element_line(colour = 'black'),
                  axis.title = element_text(colour = 'black', size = 14),
                  axis.text.x = element_text(angle = 0,  colour = 'black', size = 10),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.tag = element_text(colour = 'black', size = 18, 
                                          face = 'bold', hjust = 0))
        ),
        
        matrix = (
          intersection_matrix() +
            labs(tag = 'D') +
            theme_classic() 
        ),
        
        themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(color = c('black'),
                                                                                                    size = 14, face = 'bold'),
                                                                         panel.background = element_rect(colour = 'black'),
                                                                         panel.border = element_rect(colour = 'black', 
                                                                                                     fill = 'transparent'),
                                                                         plot.tag = element_text(colour = 'black', size = 18, 
                                                                                                 face = 'bold', hjust = 0),
                                                                         axis.ticks = element_line(colour = 'black'),
                                                                         axis.line = element_line(colour = "black", 
                                                                                                  linewidth = rel(1))))),
        upset_stripes(colors = NA)
  )

full_plot
full_plot[[3]] <- rank_plot + expand_limits(y = c(NA, 250)) + 
  theme(legend.position.inside = c(0.1, 0.65),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.title.y = element_text(vjust = -25))
full_plot[[1]] <- correlation_plot + expand_limits(y = c(NA, 0.12)) +
  theme(axis.title.y = element_text(vjust = -25))
# full_plot[[5]] <- full_plot[[5]] +
#   theme(axis.title.y = element_text(vjust = -25))
# full_plot[[7]] <- full_plot[[7]] +
#   theme(axis.title.y = element_text(vjust = -25))
full_plot[[5]] <- intersection_plot +
  expand_limits(y = c(NA, 200)) +
  theme(axis.title.y = element_text(vjust = -25))

full_plot[[8]] <- set_size_plot + expand_limits(x = c(NA, 200))
  

# ggdraw(full_plot) +
#   draw_plot(colour_options$legend, 
#             x = 0.75, y = 0.15, width = 0.25, height = 0.85)
full_plot
ggsave('../Results/Fig5_DAvML.png', 
       height = 12, width = 10, bg = 'white')
ggsave('../Results/Fig5.svg',
       height = 12, width = 10, bg = 'white')
