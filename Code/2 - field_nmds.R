library(multcomp)
library(multcompView)
library(tidyverse)
library(magrittr)
library(patchwork)
library(vegan)
library(phyloseq)
library(lubridate)
library(nlme)
library(emmeans)

#### Functions ####

diag.plots <- function(model, col.nos=c(), data=NULL, std=T, force.RvF=F) {
  library(nlme); library(mgcv); library(lmtest)
  
  get.r <- function(model.class, model) {
    if (model.class=="lm" | model.class=="aov") {
      r <- resid(model)
    }
    if (model.class=="gls" | model.class=="lme" | model.class=="gnls") {
      r <- residuals(model, type="normalized")
    }
    if (model.class=="gam") {
      r <- resid(model)
      names(r) <- rownames(model$model)
    }
    if (model.class=="nls") {
      r <- resid(model)
      names(r) <- rownames(model$model)
    }
    # if(model.class=='merModLmerTest'){
    if(model.class=='lmerMod'){
      r<-resid(model)
    }
    if(model.class == 'glmmTMB'){
      r<-resid(model)
    }
    
    return(r)
  }
  
  check.ANOVA <- function(model.class, model) {
    #Modified by JDS due to issue with being unable to identify if gls is an anova without the dataframe used present in the global environment
    #This was a problem for rendering PDFs from student assignments in a loop
    if (model.class=="lm" | model.class=="aov") {
      mm <- model$model
      is.ANOVA <- T
      for (j in 2:ncol(mm)) {if (!is.factor(mm[,j])) {is.ANOVA <- F}}
    }
    if (model.class=="gls") {
      if (length(model$contrasts)==0) {is.ANOVA <- F}
      else {
        #mm <- getData(model) #JDS Commented out - may break...who knows
        input.names <- names(model$parAssign)[-1]
        if (length(input.names)==0) {is.ANOVA <- T}
        else {
          is.ANOVA <- T
          # for (a.name in input.names) {  # assumes main terms present for interaction #JDS Commented out - may break...who knows
          #   j <- grep(a.name, colnames(mm)) #JDS Commented out - may break...who knows
          #   if (length(j) > 0 && !is.factor(mm[,j])) {is.ANOVA <- F} #JDS Commented out - may break...who knows
          # }
        }
      }
    }
    if (model.class=="lme") {
      mm <- getData(model)
      input.names <- attr(model$terms, "term.labels")
      if (length(input.names)==0) {is.ANOVA <- T}
      else {
        is.ANOVA <- T
        for (a.name in input.names) {  # assumes main terms present for interaction
          j <- grep(a.name, colnames(mm))
          if (length(j) > 0 && !is.factor(mm[,j])) {is.ANOVA <- F}
        }
      }
    }
    if (model.class=="gam" | model.class == "nls" | model.class=="gnls" | model.class=="lmerMod" | model.class == 'glmmTMB') {is.ANOVA <- F}
    return(is.ANOVA)
  }
  
  plot.continuous <-function(f, r, f.name, r.name, model.class, data=NULL, model=NULL) {
    f <- round(f, 6)
    #  if (!is.null(data)) {f <- f[rownames(data) %in% names(r)]}
    if (!is.null(model) & (model.class=="lm" | model.class=="aov")) {
      bp.p.value <- bptest(model)$p.value
    }
    else {
      whee <- lm(r~f)
      bp.p.value <- bptest(whee)$p.value
    }
    sub.text <- "too few unique x's to fit a spline"
    if (length(unique(f))>4) {
      if (length(unique(f)) < 10) {k <- round(length(unique(f))/2)} else {k <- 10}
      rf.model<-gam(r~s(f, k=k), data=data.frame(r=as.double(r),f=as.double(f)), se.fit=T)
      sub.text <- paste("spline p-val =", signif(summary(rf.model)$s.table[1,4],4))
    }
    plot(r~f, xlab=f.name, ylab=r.name, sub=sub.text,
         main=paste("B-P p-val = ", signif(bp.p.value, 4)))
    abline(h=0, col=4, lwd=2)
    if (length(unique(f))>4) {
      f.seq<-seq(min(f),max(f),(max(f)-min(f))/100)
      r.seq<-predict(rf.model,newdata=data.frame(f=f.seq), se.fit=T)
      lines(r.seq$fit ~ f.seq,col="red", lwd=2)
      lines(r.seq$fit+2*r.seq$se.fit ~ f.seq, col="red", lty=2, lwd=2)
      lines(r.seq$fit-2*r.seq$se.fit ~ f.seq, col="red", lty=2, lwd=2)    
    }
  }
  
  qqnorm.with.sim.bounds<-function(the.data,sw=T,robust=T,main=NA,legend=F) {
    
    if (is.na(main)) {main<-"Normal Q-Q Plot"}
    n<-length(the.data)
    y.mx<-matrix(0,5000,n)
    #  robust--use median and mad; !robust--use mean and sd
    if (robust) {x.center<-median(the.data)} else {x.center<-mean(the.data)}
    if (robust) {s.center<-mad(the.data)} else {s.center<-sd(the.data)}
    for (i in 1:5000) {y.mx[i,]<-sort(rnorm(n,x.center,s.center))}
    lo.bds<-apply(y.mx,2,lo.bd<-function(x) {quantile(x,0.025)})
    hi.bds<-apply(y.mx,2,hi.bd<-function(x) {quantile(x,0.975)})
    lo.lo.bds<-apply(y.mx,2,lo.lo.bd<-function(x) {quantile(x,0.005)})
    hi.hi.bds<-apply(y.mx,2,hi.hi.bd<-function(x) {quantile(x,0.995)})
    meds<-apply(y.mx,2,median)
    ideal.x<-qnorm(ppoints(n))
    if (sw) {
      plot(ideal.x,sort(the.data),ylim=c(min(c(lo.lo.bds,the.data)),max(c(the.data,hi.hi.bds))),main=main,xlab="Theoretical Quantiles",ylab="Sample Quantiles",sub=paste("SW p-value",signif(shapiro.test(the.data)$p.value,4),' - Length',length(the.data)))
    }
    else {plot(ideal.x,sort(the.data),ylim=c(min(c(lo.lo.bds,the.data)),max(c(the.data,hi.hi.bds))),main=main,xlab="Theoretical Quantiles",ylab="Sample Quantiles")}
    abline(thud<-lm(meds~ideal.x),col="blue")
    lines(ideal.x,lo.bds,col="green",lty=2)
    lines(ideal.x,hi.bds,col="green",lty=2)
    lines(ideal.x,lo.lo.bds,col="red",lty=2)
    lines(ideal.x,hi.hi.bds,col="red",lty=2)
    if (legend) legend("topleft",legend=c("median","95% bounds","99% bounds"),lty=1,col=4:2)
  }
  
  bfl.test<-function(formula,data,...) {
    yf<-all.vars(formula);y<-data[,yf[1]];f<-data[,yf[2]]
    medians<-tapply(y,f,median)
    z<-abs(y-medians[f])
    z.aov<-aov(z~f)
    boxplot(resid(aov(y~f))~f,ylab="Residuals",
            sub=paste("BFL p-value =",signif(anova(z.aov)[1,5],4)),...)
  }
  
  
  if (length(col.nos)>0 & is.null(data)) {
    cat("You must specify the name of the dataset used in the model, for example:")
    cat("\n")
    cat("diag.plots(my.model, col.nos=c(2, 4), data=my.data)")
    cat("\n"); cat("\n")
    return()
  }
  model.class <- class(model)[1]
  ##  Use model.class to dig out necessary data
  r <- get.r(model.class, model)
  is.ANOVA <- check.ANOVA(model.class, model)
  #mm <- getData(model)
  #is.ANOVA<-TRUE
  ##  Create space for plots
  n.plots <- length(col.nos)
  if (std) {
    n.plots <- n.plots + 2
    if (model.class=="lm" | model.class=="aov") {n.plots <- n.plots + 1}
    if (is.ANOVA & !force.RvF) {n.plots <- n.plots - 1}
  }
  n.plot.cols<-trunc((n.plots+1)/2)
  par(mfcol=c(2, n.plot.cols))
  ##  If desired, create standard plots
  if (std) {
    f <- fitted(model)
    if (!is.ANOVA | force.RvF) {plot.continuous(f, r, "Fitted", "Residuals", model.class, model=model)}
    #    print("here")
    if (length(r)>=5000) {qqnorm(r)} else {qqnorm.with.sim.bounds(r)}
    if (model.class=="lm") {plot(model, which=4)}
  }
  ##  Create non-standard plots
  if (length(col.nos) > 0) {
    data<-as.data.frame(data)
    for (j in col.nos) {
      if(is.character(data[,j])){data[,j]<-as.factor(data[,j])}
      if (is.factor(data[,j])) {
        f <- data[,j]
        df.for.bfl <- data.frame(r, f)
        bfl.test(r~f, data=df.for.bfl, xlab=colnames(data)[j])
      }
      else {plot.continuous(data[,j], r, colnames(data)[j], "Residuals", model.class, data=data)}
    }
  }
}

#### Data ####
field_data_cpm <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                       show_col_types = FALSE) %>%
  select(asv_id, sample_id, log2_cpm_norm) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm)


if(!file.exists('../intermediate_files/taxonomy.csv.gz')){
  taxonomy <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
           show_col_types = FALSE) %>%
    select(asv_id, domain:species) %>%
    distinct %>%
    mutate(across(everything(), str_replace_na)) %>%
    write_csv('../intermediate_files/taxonomy.csv.gz')
} else {
  taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz', 
                       show_col_types = FALSE) %>%
    mutate(across(everything(), str_replace_na))
}

field_metadata <- read_csv('../intermediate_files/normalized_field_asv_counts.csv', 
                           show_col_types = FALSE) %>%
  select(sample_id, year, season, site, health, lib.size) %>%
  distinct %>%
  mutate(timepoint = str_c(year, season, sep = '_'),
         across(where(is.character), as.factor)) %>%
  select(sample_id, site, health, timepoint, lib.size)

#### Alpha Diversity Metrics ####
microbiome_data <- read_rds("../intermediate_files/prepped_microbiome.rds.gz")

microbiome::summarize_phyloseq(microbiome_data)

alpha_metrics <- microbiome_data %>%
  subset_samples(is.na(tank)) %>%
  rarefy_even_depth(verbose = TRUE, rngseed = 5222024) %>%
  microbiome::alpha(zeroes = FALSE) %>%
  as_tibble(rownames = 'sample_id') 

alpha_metric_analysis <- alpha_metrics %>%
  right_join(field_metadata,
             by = 'sample_id') %>%
  select(sample_id, site, health, timepoint, lib.size, observed, 
         diversity_shannon, diversity_inverse_simpson) %>%
  rename(richness = observed) %>%
  pivot_longer(cols = c(richness, diversity_shannon, diversity_inverse_simpson),
               names_to = 'metric', 
               values_to = 'value') %>%
  mutate(log_value = log(value),
         sqrt_value = sqrt(value),
         inv_value = 1/value) %>%
  nest_by(metric)

## Richness ##
richness_aov <- lme(value ~ health * timepoint, 
                    random = ~1 | site,
                    data = alpha_metric_analysis$data[[3]])
diag.plots(richness_aov, col.nos = c(2:4), data = alpha_metric_analysis$data[[3]])


#Fix site/timepoint variance
richness_gls <- lme(value ~ health * timepoint, 
                    random = ~1 | site, 
                    data = alpha_metric_analysis$data[[3]],
                    weights = varComb(varIdent(form = ~1|health), 
                                      varIdent(form = ~1|timepoint), 
                                      varIdent(form = ~1|site)))

diag.plots(richness_gls, col.nos = c(2:4), data = alpha_metric_analysis$data[[3]])
car::Anova(richness_gls, type = '2')

#Try log 
richness_gls_sqrt <- lme(sqrt_value ~ health * timepoint, 
                    random = ~1 | site, 
                    data = alpha_metric_analysis$data[[3]],
                    weights = varComb(varIdent(form = ~1|health), 
                                      varIdent(form = ~1|timepoint), 
                                      varIdent(form = ~1|site)))
diag.plots(richness_gls_sqrt, col.nos = c(2:4), data = alpha_metric_analysis$data[[3]])
anova(richness_gls_sqrt)
car::Anova(richness_gls_sqrt, type = '2')


## Shannon ##
shannon_aov <- lme(value ~ health * timepoint, 
                   random = ~1 | site, 
                   data = alpha_metric_analysis$data[[2]])
diag.plots(shannon_aov, col.nos = c(2:4), data = alpha_metric_analysis$data[[2]])

shannon_gls <- lme(value ~ health * timepoint, 
                   random = ~1 | site, 
                   data = alpha_metric_analysis$data[[2]],
                    weights = varComb(varIdent(form = ~1|health), 
                                      varIdent(form = ~1|timepoint), 
                                      varIdent(form = ~1|site)))

diag.plots(shannon_gls, col.nos = c(2:4), data = alpha_metric_analysis$data[[2]])


#Best is no transformation


## Simpson ##
simpson_aov <- lme(value ~ health * timepoint, 
                   random = ~1 | site, 
                   data = alpha_metric_analysis$data[[1]])
diag.plots(simpson_aov, col.nos = c(2:4), data = alpha_metric_analysis$data[[1]])


simpson_gls_log <- lme(log_value ~ health * timepoint, 
                       random = ~1 | site,
                       data = alpha_metric_analysis$data[[1]],
                   weights = varComb(varIdent(form = ~1|health), 
                                     varIdent(form = ~1|timepoint), 
                                     varIdent(form = ~1|site)))
diag.plots(simpson_gls_log, col.nos = c(2:4), data = alpha_metric_analysis$data[[1]])

#### Plot/Outputs ####
as_tibble(car::Anova(simpson_gls_log, type = '2'), rownames = 'term')

bind_rows(
  richness = as_tibble(car::Anova(richness_gls_sqrt, type = '2'), rownames = 'term'),
  shannon = as_tibble(car::Anova(shannon_gls, type = '2'), rownames = 'term'),
  simpson = as_tibble(car::Anova(simpson_gls_log, type = '2'), rownames = 'term'),
  .id = 'metric'
) %>%
  mutate(`Chisq` = sprintf(Chisq, fmt = '%#.3f'),
         `p-value` = scales::pvalue(`Pr(>Chisq)`),
         .keep = 'unused') %>%
  filter(term != '(Intercept)') %>%
  write_csv('../Results/table2_alpha_diversity_anova.csv')


significant_metrics <- bind_rows(
  Richness = emmeans(richness_gls_sqrt, ~health | timepoint) %>%
    contrast('pairwise') %>%
    as_tibble(),
  
  Shannon = emmeans(shannon_gls, ~health | timepoint) %>%
    contrast('pairwise') %>%
    as_tibble(),
  
  `Inverse Simpson` = emmeans(simpson_gls_log, ~health | timepoint) %>%
    contrast('pairwise') %>%
    as_tibble(),
  
  .id = 'metric'
) %>%
  filter(p.value < 0.05) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '-'),
         timepoint = ymd(str_c(timepoint, '1', sep = '-'))) %>%
  select(metric, timepoint)

significant_metrics_time <- bind_rows(
  Richness = emmeans(richness_gls_sqrt, ~timepoint) %>%
    cld(Letters = LETTERS) %>%
    as_tibble(),
  
  Shannon = emmeans(shannon_gls, ~timepoint) %>%
    cld(Letters = LETTERS) %>%
    as_tibble(),
  
  `Inverse Simpson` = emmeans(simpson_gls_log, ~timepoint) %>%
    cld(Letters = LETTERS) %>%
    as_tibble(),
  
  .id = 'metric'
) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '-'),
         timepoint = ymd(str_c(timepoint, '1', sep = '-')),
         .group = str_trim(.group)) %>%
  select(metric, timepoint, .group)



ref_grid(richness_gls_sqrt) %>%
  update(tran = poisson(link = 'sqrt')) %>%
  regrid(transform = 'response') %>%
  emmeans(~health) %>%
  contrast('pairwise')

bind_rows(
  Richness = ref_grid(richness_gls_sqrt) %>%
    update(tran = poisson(link = 'sqrt')) %>%
    emmeans(~health | timepoint, type = 'response') %>%
    as_tibble() %>%
    rename(emmean = response),
  
  `Shannon` = emmeans(shannon_gls, ~health | timepoint) %>%
    as_tibble(),
  
  `Inverse Simpson` = ref_grid(simpson_gls_log) %>%
    update(tran = poisson(link = 'log')) %>%
    emmeans(~health | timepoint, type = 'response') %>%
    as_tibble() %>%
    rename(emmean = response),
  
  .id = 'metric'
) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '-'),
         timepoint = ymd(str_c(timepoint, '1', sep = '-')),
         health = str_replace_all(health, c('D' = 'Diseased', 'H' = 'Healthy'))) %>%
  ggplot(aes(x = timepoint, y = emmean, ymin = lower.CL, ymax = upper.CL,
             colour = health)) +
  geom_pointrange(position = position_dodge2(50)) +
  geom_text(data = significant_metrics,
            aes(x = timepoint, y = Inf, label = '*'), 
            size = 10, inherit.aes = FALSE,
            vjust = 1) +
  scale_colour_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, type = "continuous"),
                                         c('Healthy', 'Diseased'))) +
  # facet_wrap(~metric, scales = 'free_y')
  scale_x_date(breaks = ymd(c('2016-01-01', '2016-07-01', 
                              '2017-01-01', '2017-07-01')), 
               date_labels = '%Y\n%b') +
  facet_grid(metric ~ ., scales = 'free_y', 
             switch = 'y') +
  labs(x = NULL,
       y = NULL,
       colour = 'Disease\nState') +
  theme_classic() +
  theme(strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 14),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank(),
        # legend.position = c(0.9, 0.9),
        legend.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 12))
ggsave('../Results/alpha_diversity.png', height = 7, width = 7)


bind_rows(
  Richness = ref_grid(richness_gls_sqrt) %>%
    update(tran = poisson(link = 'sqrt')) %>%
    emmeans(~timepoint, type = 'response') %>%
    as_tibble() %>%
    rename(emmean = response),
  
  `Shannon` = emmeans(shannon_gls, ~timepoint) %>%
    as_tibble(),
  
  `Inverse Simpson` = ref_grid(simpson_gls_log) %>%
    update(tran = poisson(link = 'log')) %>%
    emmeans(~timepoint, type = 'response') %>%
    as_tibble() %>%
    rename(emmean = response),
  
  .id = 'metric'
) %>%
  mutate(timepoint = str_replace_all(timepoint, c('S' = 'July', 'W' = 'Jan')),
         timepoint = str_replace_all(timepoint, '_', '-'),
         timepoint = ymd(str_c(timepoint, '1', sep = '-'))) %>% filter(timepoint == ymd('2017-07-01'))
  
#### PERMANOVA ####
if(file.exists('../intermediate_files/field_other_adonis.rds.gz')){
  the_adonis <- read_rds('../intermediate_files/field_adonis.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  
  adon_y_data <- column_to_rownames(field_data_cpm, 'sample_id')
  
  the_adonis <- adonis2(vegdist(adon_y_data, binary = FALSE, method = 'bray') ~ 
                          health * site * timepoint, 
                        permutations = 9999, by = 'terms',
                        data = field_metadata, parallel = clus)
  
  stopCluster(cl = clus)
  write_rds(the_adonis, '../intermediate_files/field_adonis.rds.gz')
}

the_adonis

#### NMDS ####
if(file.exists('../intermediate_files/field_tank_nmds.rds.gz')){
  the_nmds <- read_rds('../intermediate_files/field_tank_nmds.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  the_nmds <- column_to_rownames(field_data_cpm, 'sample_id') %>%
    metaMDS(distance = 'bray', binary = FALSE, trymax = 1000, 
            parallel = clus)
  stopCluster(cl = clus)
  write_rds(the_nmds, '../intermediate_files/field_tank_nmds.rds.gz')
}

#### Plot ####
colony_points_nmds <- scores(the_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(field_metadata, by = 'sample_id')


colony_points_nmds %>% 
  mutate(timepoint = str_replace_all(timepoint, c('W' = 'Jan', 'S' = 'Jul')),
         health = if_else(health == 'D', 'Diseased', 'Healthy')) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(shape = site, colour = health, fill = health), size = 1.5) +
  # facet_wrap(~timepoint, labeller = labeller(timepoint = ~str_replace(., '_', ' - '))) +
  facet_wrap(~timepoint, labeller = labeller(timepoint = c('2016_Jan' = 'A', '2016_Jul' = 'B',
                                                           '2017_Jan' = 'C', '2017_Jul' = 'D'))) +
  scale_shape_manual(values = c('CK14' = 'circle filled', 'CK4' = 'square filled', 'HS' = 'diamond filled', 
                                'SB' = 'triangle filled', 'TS' = 'triangle down filled')) +
  scale_colour_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, type = "continuous"),
                                         c('Healthy', 'Diseased'))) +
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, type = "continuous"),
                                         c('Healthy', 'Diseased'))) +
  guides(fill = 'none',
         shape = guide_legend(override.aes = list(fill = 'black', size = 4)),
         colour = guide_legend(override.aes = list(size = 4))) +
  labs(colour = 'Disease\nState',
       shape = 'Site') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(hjust = 0, colour = 'black', size = 18, face = 'bold'),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 12))
ggsave('../Results/FigS1_NMDS_field.png', height = 10, width = 10)
