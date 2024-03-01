
#### All Random Model ####
random_models <- field_data %>%
  partition(cluster) %>%
  mutate(model = list(lmer(log2_cpm_norm ~ health + (1 | timepoint) + (1 | site),
                           data = data)),
         anova_table = list(anova(model, ddf = "Kenward-Roger")),
         d_h = as_tibble(emmeans(model, 'pairwise'~health)$contrast)$estimate, #Diseased - Healthy
         nDF = anova_table$NumDF,
         dDF = anova_table$DenDF,
         fvalue = anova_table$`F value`,
         pvalue = anova_table$`Pr(>F)`) %>%
  collect %>%
  ungroup


#### Get Significance ####
random_models %>%
  mutate(padjust = p.adjust(pvalue, method = 'fdr')) %>%
  filter(padjust < alpha) %>%
  group_by(more_disease = d_h > 0) %>%
  summarise(n = n(),
            across(c(domain:asv_id), n_distinct))

random_models %>%
  mutate(padjust = p.adjust(pvalue, method = 'fdr')) %>%
  filter(padjust < alpha) %>%
  group_by(more_disease = d_h > 0, across(domain:genus)) %>%
  summarise(n = n()) %>%
  ungroup %>%
  mutate(higher_taxa = str_c(class, order, family, genus, sep = '; '),
         n = if_else(more_disease, n, -1 * n),
         higher_taxa = fct_reorder(higher_taxa, abs(n), .fun = sum)) %>%
  ggplot(aes(y = higher_taxa, x = n, fill = more_disease)) +
  geom_col() 



random_models %>%
  mutate(padjust = p.adjust(pvalue, method = 'fdr')) %>%
  filter(padjust < alpha) %>%
  group_by(more_disease = d_h > 0, across(domain:genus)) %>%
  summarise(n = n()) %>%
  ungroup %>%
  mutate(higher_taxa = str_c(class, order, family, genus, sep = '; '),
         higher_taxa = fct_reorder(higher_taxa, abs(n), .fun = sum)) %>%
  ggplot(aes(y = higher_taxa, x = n, fill = more_disease)) +
  geom_linerange(aes(xmin = 0, xmax = n), position = position_dodge(0.5),
                 linetype = 'dashed') +
  geom_point(shape = 21, position = position_dodge(0.5), size = 5) +
  scale_x_continuous(trans = scales::log10_trans())

select(random_models, -where(is.list)) %>%
  mutate(fdr = p.adjust(pvalue, method = 'fdr'),
         fdr = fdr < alpha,
         positive = d_h > 0) %>%
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        annotations = list(
          'Order' = ggplot(mapping = aes(fill = order)) +
            geom_bar(stat = 'count', position = 'fill') +
            scale_y_continuous(labels=scales::percent_format())
        )
  )


#### Fixed Model ####
fixed_models <- field_data %>%
  partition(cluster) %>%
  mutate(model = list(gls(log2_cpm_norm ~ health * timepoint * site,
                          data = data)),
         
         anova(model, type = 'marginal') %>%
           as_tibble(rownames = 'term') %>%
           select(term, `p-value`) %>%
           filter(term != '(Intercept)') %>%
           pivot_wider(names_from = term, values_from = `p-value`) %>%
           rename_with(~str_replace_all(., ':', 'X'))) %>%
  collect %>%
  ungroup


## Upset Plot ##


?upset

select(fixed_models, -where(is.list)) %>%
  mutate(across(health:healthXtimepointXsite, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepointXsite, ~. < alpha)) %>%
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        annotations = list(
          'Order' = ggplot(mapping = aes(fill = order)) +
            geom_bar(stat = 'count', position = 'fill') +
            scale_y_continuous(labels=scales::percent_format())
        )
  )

select(fixed_models, -where(is.list)) %>%
  mutate(across(health:healthXtimepointXsite, ~p.adjust(., method = 'fdr')),
         across(health:healthXtimepointXsite, ~. < alpha)) %>%
  filter(if_any(contains('health'))) 
