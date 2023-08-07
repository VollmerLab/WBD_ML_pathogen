setwd("/Users/emilytrytten/Desktop/Screenshots/Career/Vollmer Lab/GitHub/Panama_Tank_Field/Code")


pan_field_dat <- read_csv("../intermediate_files/normalized_field_asv_counts.csv")

pan_tank_dat <- read_csv("../intermediate_files/normalized_tank_asv_counts.csv")

#### FIELD DATA ####

f_dat_trimmed <- pan_field_dat %>%
  select(-c(dataset, cpm_norm, cpm, n_reads, lib.size, norm.factors)) %>%
  mutate(sampling_time = case_when(
    year == 2016 & season == "W" ~ 1,
    year == 2016 & season == "S" ~ 2,
    year == 2017 & season == "W" ~ 3,
    year == 2017 & season == "S" ~ 4
  ), .after = site) %>%
  mutate(sampling_time = factor(sampling_time))

pan_field_models <- f_dat_trimmed %>%
  #filter(asv_id %in% c("ASV1", "ASV100", "ASV10")) %>%
  nest_by(across(c('asv_id', domain:genus))) %>%
  mutate(model = list(lmer(log2_cpm_norm ~ health*sampling_time + (1 | site),
                   data))) %>%
  ungroup()


pan_field_plots <- pan_field_models %>%
  #filter(asv_id %in% c("ASV1", "ASV100", "ASV10")) %>%
  rowwise %>%
  mutate(plot_info = list(emmeans(model, ~health*sampling_time) %>%
                            broom::tidy(conf.int = TRUE) %>%
                            group_by(sampling_time) %>%
                            mutate(sig_placement = max(conf.high)) %>%
                            ungroup()
                          )) %>%
  mutate(pw_plot_info = list(emmeans(model, ~health | sampling_time) %>%
                               pairs() %>%
                            broom::tidy(conf.int = TRUE) %>%
                              mutate(fdr_p = p.adjust(p.value, method = "fdr"),
                                     sig = ifelse(fdr_p < 0.05, "sig", "nonsig")) %>%
                              left_join(plot_info %>% select(sampling_time, sig_placement) %>% distinct(),
                                        by = join_by(sampling_time))
                            )) %>%
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = sampling_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
      geom_pointrange(aes(col = health), position = position_dodge(width = 0.5)) +
      scale_color_manual(values = c("D" = "#F10A0A",
                                    "H" = "#29A5B2")) +
      scale_x_discrete(breaks=c("1","2","3", "4"),
                       labels=c("Winter 2016", "Summer 2016", "Winter 2017", "Summer 2017")) +
      geom_point(data = pw_plot_info, aes(x = sampling_time, y = sig_placement + 1, alpha = sig), pch = "*", size = 5, inherit.aes = FALSE) +
      scale_alpha_manual(values = c("nonsig" = 0,
                                    "sig" = 1),
                         guide = "none") +
      theme_bw() +
      ggtitle(paste(asv_id, "-", family, genus, sep = " "))
    ))
      
pan_field_plots %>% filter(asv_id == "ASV8") %>% pull(plot)   

#### TANK DATA ####

t_dat_trimmed <- pan_tank_dat %>%
  select(-c(dataset, cpm_norm, cpm, n_reads, lib.size, norm.factors, site, year, season)) %>%
  mutate(treatment = paste(time, anti_health, sep = "_")) %>%
  nest_by(across(c('asv_id', domain:genus))) %>%
  mutate(model = list(lmer(log2_cpm_norm ~ treatment + (1 | tank) + (1 | plate) + (1 | geno),
                           data))) %>%
  ungroup()

remove_graph <- t_dat_trimmed %>% unnest(data) %>% group_by(asv_id, treatment)  %>%
  mutate(tot = n()) %>% filter(log2_cpm_norm < 6.12) %>% reframe(zeroes = n(), tot) %>% 
  mutate(perc = 100*zeroes/tot) %>% filter(perc > 95) %>% ungroup() %>% group_by(asv_id) %>%
  summarize(num_gone = n()) %>% filter(num_gone > 90) %>% pull(asv_id)


pan_tank_plots <- t_dat_trimmed %>%
  filter(!asv_id %in% remove_graph) %>% #lots of 0s, probs won't graph - recheck specific ASVs if applicable
  rowwise %>%
  mutate(plot_info = list(emmeans(model, ~treatment) %>%
                            broom::tidy(conf.int = TRUE) %>%
                            separate(treatment, into = c('time', 'timecat', 'antibiotic', 'health')) %>%
                            mutate(graph_cat = ifelse(time == "0", NA, 
                                                      paste(antibiotic, health, sep = "_"))) %>%
                            {. ->> intermed } %>%
                            mutate(graph_cat = ifelse(time == "0", paste(antibiotic, "D", sep = "_"), 
                                                      graph_cat)) %>%
                            dplyr::slice(rep(1:2, 1)) %>%
                            rbind(intermed) %>%
                            mutate(graph_cat = ifelse(is.na(graph_cat), paste(antibiotic, "N", sep = "_"), 
                                                      graph_cat)) %>%
                            mutate(c_time = parse_number(time)) %>%
                            mutate(facet_lab = "Experimental") %>%
                            mutate(c_time = ifelse(time == "0", ifelse(antibiotic == "A", c_time - 0.15, c_time + 0.15),
                                                   case_when(graph_cat == "A_D" ~ c_time - 0.35,
                                                             graph_cat == "N_D" ~ c_time - 0.15,
                                                             graph_cat == "A_N" ~ c_time + 0.15,
                                                             graph_cat == "N_N" ~ c_time + 0.35)))
         )) %>%
  mutate(plot = list(
    ggplot(data = plot_info, aes(x = c_time, y = estimate, ymin = conf.low, ymax = conf.high)) +
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("A_D", "N_D"))), aes(colour1 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_line(data = (plot_info %>% filter(graph_cat %in% c("A_N", "N_N"))), aes(colour2 = graph_cat, linetype = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("A_D", "N_D"))), width = 0, aes(colour1 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_errorbar(data = (plot_info %>% filter(graph_cat %in% c("A_N", "N_N"))), width = 0, aes(colour2 = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("A_D", "N_D"))), size = 3, aes(colour1 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour1"))) + 
      (geom_point(data = (plot_info %>% filter(graph_cat %in% c("A_N", "N_N"))), size = 3, aes(colour2 = graph_cat, pch = graph_cat)) %>%
         rename_geom_aes(new_aes = c("colour" = "colour2"))) +
      scale_color_manual(aesthetics = "colour1", values = c("A_D" = "#ED023B", "N_D" = "#FF7ABD"), guide = "legend", 
                         name = "Diseased", labels = c("Antibiotics", "Control")) +
      scale_shape_manual(values = c("A_D" = 15, "A_N" =  15, "N_D" = 16, "N_N" = 16), guide = "none") +
      scale_color_manual(aesthetics = "colour2", values = c("A_N" = "#03975B", "N_N" = "#40E8A4"), guide = "legend", 
                         name = "Healthy", labels = c("Antibiotics", "Control")) +
      guides(colour1 = guide_legend(
        override.aes=list(linetype = c(1, 6), shape = c(15, 16))),
        colour2 = guide_legend(
          override.aes=list(linetype = c(1, 6), shape = c(15, 16))),
        ) +
      scale_x_continuous(breaks=c(0, 2, 8)) +
      scale_linetype_manual(values = c("A_D" = 1, "A_N" =  1, "N_D" = 6, "N_N" = 6), guide = "none") +
      theme_bw() +
      xlab("Time") +
      ylab(expression("Normalized log"[2]*" (cpm)")) +
      labs(title = str_c(family, " ", genus, " (", asv_id, ")", sep = ""))
  ))

pan_tank_plots %>% filter(asv_id == "ASV25") %>% pull(plot)

pan_tank_plots %>% filter(asv_id == "ASV25") %>% pull(plot_info)

#write_rds(pan_tank_plots, "../intermediate_files/pan_tank_plots.rds")

pan_field_plots %>% filter(family == "Francisellaceae")
