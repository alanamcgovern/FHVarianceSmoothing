## results script for settings that use same surface for each simulation

library(tidyverse)
library(readr)
library(knitr)
library(ggpubr)
library(ggh4x)

######### SUMMARY OF DIAGNOSTICS --------------------
diag_full <- NULL
for(setting in c(1:5,7)){
#for(setting in c(1,'1_new')){
  setwd(paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/Sim',setting))
  tmp <- read_csv('Summary_diag.csv')
  tmp$setting <- setting
  diag_full <- rbind(diag_full,tmp)
}

diag_full %>% group_by(setting) %>% summarise(sum(divergences))

diag_full %>% group_by(model) %>% summarise(sum(max_treedepth))

diag_full %>% filter(max_treedepth>0)

diag_full %>% filter(!(model %in% c('Mod4','Mod4a'))) %>% group_by(model) %>% summarise(mean(max_rhat > 1.01))


sum(diag_full$divergences)
diag_full %>% filter(model %in% c('Mod1','Mod2','Mod1a','Mod2a')) %>% group_by(setting,model) %>% 
  summarise(minb = min(bulk_ess),avgb = mean(bulk_ess),mint = min(tail_ess),avgt = mean(tail_ess)) %>% 
  summarise(min(minb),min(avgb),min(mint),min(avgt))

diag_full%>% filter(model %in% c('Mod1','Mod2','Mod1a','Mod2a')) %>% ggplot() + 
  geom_histogram(aes(bulk_ess)) +
  facet_wrap(~model)

diag_full%>% filter(model %in% c('Mod1','Mod2','Mod1a','Mod2a')) %>% ggplot() + 
  geom_histogram(aes(tail_ess)) +
  facet_wrap(~model)

diag_full%>% filter(model %in% c('Mod1','Mod2','Mod1a','Mod2a')) %>% ggplot() + 
  geom_histogram(aes(min_efbmi)) +
  facet_wrap(~model)

########### SUMMARY OF SAMPLES -----------
setting <- c(1,4)[2]
setwd(paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/Sim',setting))
load(file=paste0('sampled_clusters.rda'))

sample_info_list <- lapply(1:100,function(j){
  tmp <- sampled_clusters[[j]] %>% group_by(admin2) %>% summarise(m_u = sum(urban),m_r = sum(1-urban),n_u = sum(urban*n),n_r = sum((1-urban)*n))
  tmp <- left_join(data.frame(admin2=1:300),tmp,by='admin2') %>% mutate(across(where(is.numeric), ~replace_na(., 0)))
  tmp$sim <- j
  return(tmp)
})
sample_info <- do.call(rbind,sample_info_list)

# sampled clusters and observations per admin2 (urban and rural)
sample_info %>% #filter(m_u + m_r >0) %>% 
  summarise(mean(m_u),
                                                     sd(m_u),
                                                     mean(m_r),
                                                     sd(m_r),
                                                     mean(n_u),
                                                     sd(n_u),
                                                     mean(n_r),
                                                     sd(n_r))
# average number of admin2 with no sampled clusters
sample_info %>% group_by(sim) %>% summarise(x = sum(m_u+m_r==0)) %>% summarise(mean(x),sd(x))
# average number of admin2 with less than 5 sampled clusters
sample_info %>% group_by(sim) %>% summarise(x = sum(m_u+m_r < 5)) %>% summarise(mean(x),sd(x))


########## ACCURACY OF APPROXIMATE DISTRIBUTIONS ----------------

# from Comparison_Models.R
load('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/moments.rda') 

model_labs <- c('draws_srs' = 'SRS','draws_CSS' = "CSS", 'draws_CSM' = "CSM", 'draws_emp'='Empirical')
model_colors <- c('draws_srs' = 'firebrick','draws_CSS' = "#1f78b4", 'draws_CSM' = "#6a3d9a",'draws_emp'="#A1D99B")

setting_key <- data.frame(setting = c(1:5),
                          setting_labels = c('1: correctly specified',
                                             '2: varying urban effects',
                                             '3: t-distributed data',
                                             '1a: 5x sampled clusters',
                                             '4: within cluster correlation'))
setting_levels <- c('1: correctly specified',
                    '1a: 5x sampled clusters',
                    '2: varying urban effects',
                    '3: t-distributed data',
                    '4: within cluster correlation')

first_moment <- merge(setting_key, moments_list[[1]])
second_moment <- merge(setting_key, moments_list[[2]])
third_moment <- merge(setting_key, moments_list[[3]])

  m1 <- first_moment %>% filter(setting ==1) %>%
    pivot_longer(cols = c(draws_srs,draws_CSS,draws_CSM,draws_emp),names_to = 'model',values_to = 'value') %>%
    mutate(model = factor(model, levels = c('draws_srs', 'draws_CSS', 'draws_CSM', 'draws_emp'))) %>%
    ggplot(aes(model,value)) + 
    geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
    scale_y_continuous(trans='log', breaks = c(0.005,0.01,0.025,0.05,0.1,0.2)) +
    scale_fill_manual(name = 'Distribution',values = model_colors,labels=model_labs) + 
    labs(x=NULL,y=NULL) + theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=8)) +
    ggtitle('A. Mean')
  
  m2 <- second_moment %>% filter(setting ==1) %>%
    pivot_longer(cols = c(draws_srs,draws_CSS,draws_CSM,draws_emp),names_to = 'model',values_to = 'value') %>%
    mutate(model = factor(model, levels = c('draws_srs', 'draws_CSS', 'draws_CSM', 'draws_emp'))) %>%
    ggplot(aes(model,value)) + 
    geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
    scale_y_continuous(trans='log', breaks = c(1e-6,1e-5,1e-4,1e-3,1e-2)) +
    scale_fill_manual(name = 'Distribution',values = model_colors,labels=model_labs) + 
    labs(x=NULL,y=NULL) + theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=8)) +
    ggtitle('B. Variance')
  
  m3 <- third_moment %>% filter(setting ==1) %>%
    pivot_longer(cols = c(draws_srs,draws_CSS,draws_CSM,draws_emp),names_to = 'model',values_to = 'value') %>%
    mutate(model = factor(model, levels = c('draws_srs', 'draws_CSS', 'draws_CSM', 'draws_emp'))) %>%
    ggplot(aes(model,value)) + 
    geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
    scale_y_continuous(trans='log', breaks = c(1,2,5,10,20,50)) +
    scale_fill_manual(name = 'Distribution',values = model_colors,labels=model_labs) + 
    labs(x=NULL,y=NULL) + theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=8)) +
    ggtitle('C. Skewness')
  
  
  # 875 x 350
  ggarrange(plotlist = list(m1,m2,m3), ncol=3, common.legend = T, legend='bottom')
  
  
# appendix

  m1 <- first_moment %>% 
    pivot_longer(cols = c(draws_srs,draws_CSS,draws_CSM,draws_emp),names_to = 'model',values_to = 'value') %>%
    mutate(model = factor(model, levels = c('draws_srs', 'draws_CSS', 'draws_CSM', 'draws_emp'))) %>%
    ggplot(aes(model,value)) + 
    geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
    scale_y_continuous(trans='log', breaks = c(2e-3,5e-3,2e-2,5e-2,2e-1)) +
    scale_fill_manual(name = 'Distribution',values = model_colors,labels=model_labs) + 
    labs(x=NULL,y=NULL) + theme_bw() +
    facet_wrap(~factor(setting_labels,
                       levels = setting_levels),nrow=1) + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=8)) +
    ggtitle('A. Mean')

  m2 <- second_moment %>% 
    pivot_longer(cols = c(draws_srs,draws_CSS,draws_CSM,draws_emp),names_to = 'model',values_to = 'value') %>%
    mutate(model = factor(model, levels = c('draws_srs', 'draws_CSS', 'draws_CSM', 'draws_emp'))) %>%
    ggplot(aes(model,value)) + 
    geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
    scale_y_continuous(trans='log', breaks = c(1e-8,1e-6,1e-4,1e-2)) +
    scale_fill_manual(name = 'Distribution',values = model_colors,labels=model_labs) + 
    labs(x=NULL,y=NULL) + theme_bw() +
    facet_wrap(~factor(setting_labels,
                       levels = setting_levels),nrow=1) + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=8)) +
    ggtitle('B. Variance')
  
  m3 <- third_moment %>% 
    pivot_longer(cols = c(draws_srs,draws_CSS,draws_CSM,draws_emp),names_to = 'model',values_to = 'value') %>%
    mutate(model = factor(model, levels = c('draws_srs', 'draws_CSS', 'draws_CSM', 'draws_emp'))) %>%
    ggplot(aes(model,value)) + 
    geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
    scale_y_continuous(trans='log', breaks = c(0.5,2,5,20,50)) +
    scale_fill_manual(name = 'Distribution',values = model_colors,labels=model_labs) + 
    labs(x=NULL,y=NULL) + theme_bw() +
    facet_wrap(~factor(setting_labels,
                       levels = setting_levels),nrow=1) + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=8)) +
    ggtitle('C. Skewness')


# 950 x 615
ggarrange(plotlist = list(m1,m2,m3), nrow=3, common.legend = T, legend='bottom')


########## ACCURACY OF EMPIRICAL DESIGN VARIANCE ----------

mc_boot <- NULL
for(setting in c(1:5)){
  setwd(paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/Sim',setting))
  
  load(file=paste0('params.rda'))
  load(file=paste0('objects.rda'))
  load(file=paste0('direct_all.rda'))
  n_admin2 <- length(params$admin2_mean)
  
  mc_err <- sapply(1:n_admin2, function(i){
  m_tmp <- unlist(lapply(direct_all, function(x){x$mean[i]}))
  m_tmp <- m_tmp[!is.na(m_tmp)]
  sd(replicate(1000,var(sample(m_tmp,length(m_tmp),replace=T))))
  })
  adm2.dir.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct_all, function(x){x$mean[i]})),na.rm = T)})
  
  mc_boot <- rbind(mc_boot, data.frame(setting = setting, 
                                       area = 1:n_admin2,
                                       mc_err, 
                                       adm2.dir.variance,
                                       cv = mc_err/adm2.dir.variance))
}

mc_boot <- merge(mc_boot,setting_key)

# 875 x 395
mc_boot %>% ggplot + geom_boxplot(aes(y=100*cv),outlier.size = 0.5, median.linewidth = 0.5,
                                  box.linewidth = 0.5,whisker.linewidth = 0.5) +
  facet_wrap(~factor(setting_labels,
                     levels = setting_levels),nrow=1) +
  scale_y_continuous(trans = 'log',breaks = c(1,2,5,10)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab('Percent error')

############ MODEL PERFORMANCE -------------------------------

df_metrics <- var.results <- NULL
for(setting in c(1:5,7)){
#for(setting in c(1,'1_new')){
  setwd(paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/Sim',setting))
  
  load(file=paste0('params.rda'))
  load(file=paste0('objects.rda'))
  load(file=paste0('direct.rda'))
  load(file=paste0('sampled_clusters.rda'))
  n_admin2 <- length(params$admin2_mean)
  adm2.dir.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct, function(x){x$mean[i]})),na.rm = T)})
  
  # mean estimates
  results.cluster <- read_csv('Summary.csv')
  results <- results.cluster[,c('model','sim','area','mean','lower','upper')]
  
  results$true_mean <- params$admin2_mean[results$area]
  
  results <- results %>% mutate(ciwidth = upper - lower,
                                se = (mean-true_mean)^2,
                                cov = lower<true_mean & upper>true_mean,
                                i.score = upper - lower + 20*(lower-true_mean)*I(true_mean<lower) + 20*(true_mean-upper)*I(true_mean>upper))
  
  df_rmse <- results %>% 
    group_by(model, 
             area) %>% 
    summarise(rmse = sqrt(mean(se)), .groups = "drop")
  
  df_cov <- results %>% 
    group_by(model, 
             area) %>% 
    summarise(cov = mean(cov), .groups = "drop")
  
  df_width <- results %>% 
    group_by(model, 
             area) %>% 
    summarise(mean_width = mean(upper - lower), .groups = "drop")
  
  df_iscore <- results %>% 
    group_by(model, 
             area) %>% 
    summarise(i.score = mean(i.score), .groups = "drop")
  
  df_metrics_tmp <- merge(merge(df_rmse,df_cov),merge(df_width,df_iscore))
  df_metrics_tmp$setting <- setting
  df_metrics <- rbind(df_metrics,df_metrics_tmp)
  
  # variance estimation

  ### combine design variance estimates
  var.results_tmp <- read.csv('Summary_Var.csv') %>% rename(v_est = median) %>% dplyr::select(-c(mean))
  standard_est <- lapply(1:100,function(k){
    tmp <- direct[[k]][,c('admin2','variance')] %>% rename(v_est=variance,area=admin2)
    tmp$sim <- k
    return(tmp)
  })
  standard_est <- do.call('rbind',standard_est)
  standard_est <- standard_est[!is.na(standard_est$v_est),]
  standard_est <- standard_est[standard_est$v_est > 1e-10,]
  standard_est$model <- 'Mod0'
  var.results_tmp <- rbind(var.results_tmp,standard_est[,c('sim','area','model','v_est')])

  ### add sigma2 estimates
  hyper.results <- read.csv('Summary_hyper.csv')
  sig2.results_tmp <- hyper.results %>% filter(str_detect(hyper.results$param,'sig2')) %>% mutate(median = exp(median))
  sig2.results_tmp$area <- as.numeric((gsub("log_sig2\\[|\\]", "",sig2.results_tmp$X)))
  sig2.results_tmp <- sig2.results_tmp[,c('sim','area','model','median')] %>% rename(sig2_est = median)

  var.results_tmp <- merge(var.results_tmp,sig2.results_tmp,by=c('sim','area','model'),all=T)

  ### back calculate constant
  var.results_tmp$cons <- var.results_tmp$v_est/var.results_tmp$sig2_est

  ### add true values
  var.results_tmp$sig2_truth <- params$sig2_pop[var.results_tmp$area]
  var.results_tmp$v_truth <- adm2.dir.variance[var.results_tmp$area]
  
  ### add theoretical variance
  var.results_tmp$v_theo <- var.results_tmp$sig2_truth*var.results_tmp$cons

  var.results_tmp$setting <- setting
  var.results <- rbind(var.results,var.results_tmp)

}

setting_key <- data.frame(setting = c(1:5,7),
                          setting_labels = c('1: correctly specified',
                                             '2: varying urban effects',
                                             '3: t-distributed data',
                                             '1a: 5x sampled clusters',
                                             '4: within cluster correlation',
                                             '5: re-enumeration'))
setting_levels <- c('1: correctly specified',
                    '1a: 5x sampled clusters',
                    '2: varying urban effects',
                    '3: t-distributed data',
                    '4: within cluster correlation',
                    '5: re-enumeration')

df_metrics <- merge(df_metrics,setting_key)
var.results <- merge(var.results,setting_key)

model_labs <- c(
  Mod1 = "CSS-struct",
  Mod2 = "CSM-struct",
  Mod1a = "CSS-unstruct",
  Mod2a = "CSM-unstruct",
  Mod3 = "SRS-struct",
  Mod3a = "SRS-unstruct",
  Mod0_oracle = "Oracle",
  Mod0= "Standard")

model_colors <- c(
  "Mod1"   = "#1f78b4",  # dark blue
  "Mod1a"  = "#a6cee3",  # light blue
  
  "Mod2"   = "#6a3d9a",  # dark orange
  "Mod2a"  = "#cab2d6",  # light orange
  
  "Mod3"   = "firebrick",  # dark orange
  "Mod3a"  = "darksalmon",  # light orange
  
  "Mod0_oracle" = "#A1D99B",  # dark purple
  "Mod0"    = "gold2"   # light purple
)

model_order <- c('Mod2','Mod2a','Mod1','Mod1a','Mod3','Mod3a','Mod0','Mod0_oracle')
var.results$model <- factor(var.results$model, levels=model_order)
df_metrics$model <- factor(df_metrics$model, levels=model_order)

## plots  ------

# 910 x 400
common_theme <- theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8),
        legend.position = 'bottom') 

# compare estimated to theoretical to true design variance ------
# v1 <- var.results %>% filter(setting %in% 1:5) %>% 
#   group_by(setting_labels,model,area) %>% 
#   summarise(val = mean(v_theo/v_truth, na.rm=T)) %>% 
#   ggplot(aes(model,val)) +
#   geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
#   stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
#   scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
#   geom_hline(yintercept = 1,col='red',lty=5) +
#   facet_wrap(~factor(setting_labels,
#                      levels = setting_levels),nrow=1) +
#   scale_y_continuous(breaks = c(0.25,0.5,0.75,1)) +
#   labs(x=NULL,y=NULL)+ common_theme + ggtitle('A. Theoretical-to-truth ratio')
# 
# v2 <- var.results %>% filter(setting %in% 1:5) %>% 
#   group_by(setting_labels,model,area) %>% 
#   summarise(val = mean(v_est/v_theo,na.rm=T))  %>% 
#   ggplot(aes(model,val)) +
#   geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
#   stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
#   scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
#   geom_hline(yintercept = 1,col='red',lty=5) +
#   facet_wrap(~factor(setting_labels,
#                      levels = setting_levels),nrow=1) +
#   scale_y_continuous(trans='log',breaks = c(1,2.5,5,10)) +
#   labs(x=NULL,y=NULL)+ common_theme + ggtitle('B. Estimate-to-theoretical ratio')

# 875 x 350
var.results %>% filter(setting %in% 1:5)%>% group_by(setting_labels,model,area,v_truth) %>% 
   summarise(v_ratio=mean(v_est/v_truth,na.rm=T)) %>% ggplot(aes(model,v_ratio)) +
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  geom_hline(yintercept = 1,col='red',lty=5) +
   facet_wrap(~factor(setting_labels,
                      levels = setting_levels),nrow=1) + 
   scale_y_continuous(trans='log',breaks = c(0.1,0.25,0.5,1,2.5,5,10)) +
  coord_cartesian(ylim = c(0.1, 10)) +
  labs(x=NULL,y=NULL)+ common_theme

# 860 x 685
# my_legend <- get_legend(v3)
# ggarrange(plotlist = list(v1,v2,v3),nrow = 3,common.legend = T,
#           legend.grob = my_legend, legend = 'bottom')

# compare mean and interval estimates (Fig 5) ---------

d1 <- df_metrics%>% filter(setting %in% 1:5) %>% ggplot(aes(model,rmse)) + ggtitle('A. RMSE')+ 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks=c(0.05,0.1,0.25,0.5,1)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  facet_wrap(~factor(setting_labels,
                    levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),axis.text.x = element_blank(),
        axis.text.y = element_text(size=6))

d2 <- df_metrics%>% filter(setting %in% 1:5) %>% ggplot(aes(model,100*cov)) +  
  ggtitle('B. Coverage of 90% credible intervals (clipped at 50%)') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  geom_hline(yintercept = 90,col='red',lty=5) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  ylim(c(50,100)) +
#  coord_cartesian(ylim = c(50, 100)) +
  facet_wrap(~factor(setting_labels,
                    levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),axis.text.x =  element_blank(),
        axis.text.y = element_text(size=6))

d3 <- df_metrics%>% filter(setting %in% 1:5) %>% ggplot(aes(model,mean_width)) + ggtitle('C. Average 90% credible interval width')+
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  labs(x=NULL,y=NULL) + theme_bw() +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  facet_wrap(~factor(setting_labels,
                    levels = setting_levels),nrow=1) + 
  coord_cartesian(ylim = c(0.15,1.5)) +
  theme(plot.title = element_text(size=12),
    axis.text.x = element_blank(),
        axis.text.y = element_text(size=6))

d4 <- df_metrics%>% filter(setting %in% 1:5) %>% ggplot(aes(model,i.score)) + ggtitle('D. Average interval score for 90% credible intervals') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks = c(0.25,0.5,1,2,4)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL)+ theme_bw() +
  coord_cartesian(ylim = c(0.15,10)) +
  facet_wrap(~factor(setting_labels,
                    levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),axis.text.x =element_blank(),
        axis.text.y = element_text(size=6))

#860 x 915
ggarrange(plotlist = list(d1,d2,d3,d4),nrow=4,common.legend = T,legend = "bottom")

## plots for appendix -----------
setting_key <- data.frame(setting = c(1,7),
                          setting_labels_new = c('1: listed cluster sizes are correct',
                                             '5: cluster sizes are re-enumerated'))

df_metrics <- merge(df_metrics,setting_key)

d1 <- df_metrics%>% filter(setting %in% c(1,7)) %>% ggplot(aes(model,rmse)) + ggtitle('A. RMSE')+ 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks=c(0.05,0.1,0.25,0.5,1)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  facet_wrap(~setting_labels) +
  theme(plot.title = element_text(size=12),axis.text.x = element_blank(),
        axis.text.y = element_text(size=6))

d2 <- df_metrics%>% filter(setting %in% c(1,7)) %>% ggplot(aes(model,100*cov)) +  
  ggtitle('B. Coverage of 90% credible intervals (clipped at 50%)') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  geom_hline(yintercept = 90,col='red',lty=5) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  ylim(c(50,100)) +
  #  coord_cartesian(ylim = c(50, 100)) +
  facet_wrap(~setting_labels) +
  theme(plot.title = element_text(size=12),axis.text.x =  element_blank(),
        axis.text.y = element_text(size=6))

d3 <- df_metrics%>% filter(setting %in% c(1,7)) %>% ggplot(aes(model,mean_width)) + ggtitle('C. Average 90% credible interval width')+
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  labs(x=NULL,y=NULL) + theme_bw() +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  facet_wrap(~setting_labels) +
  coord_cartesian(ylim = c(0.15,1.5)) +
  theme(plot.title = element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6))

d4 <- df_metrics%>% filter(setting %in% c(1,7)) %>% ggplot(aes(model,i.score)) + ggtitle('D. Average interval score for 90% credible intervals') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks = c(0.25,0.5,1,2,4)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL)+ theme_bw() +
  coord_cartesian(ylim = c(0.15,10)) +
  facet_wrap(~setting_labels) +
  theme(plot.title = element_text(size=12),axis.text.x =element_blank(),
        axis.text.y = element_text(size=6))

#850 x 550
ggarrange(plotlist = list(d1,d2,d3,d4),nrow=2,ncol=2,common.legend = T,legend = "bottom")

#################################
results %>% group_by(model,area,true_mean) %>% summarise(sim_mean = mean(mean)) %>% ggplot(aes(true_mean,sim_mean)) + 
  geom_point() + facet_wrap(~model) + geom_abline(intercept = 0,slope =1)


low_cov_areas <- (df_cov %>% filter(cov<0.25) %>% group_by(area) %>% summarise(mean(cov)))$area

var.results %>% filter(model=='Mod0') %>% group_by(area) %>% reframe(n=n(),m=mean(v_est),s=sd(v_est)) %>% 
  group_by(area %in% low_cov_areas) %>% summarise(mean(n),mean(m),mean(s))

data.frame(area = 1:300, true_mean = params$admin2_mean) %>% mutate(sq_diff = sqrt((true_mean + 0.9)^2)) %>% 
  group_by(area %in% low_cov_areas) %>% summarise(mean(sq_diff))


data.frame(area = 1:300, true_var = adm2.dir.variance) %>% 
  mutate(sq_diff = sqrt((mean(log(adm2.dir.variance)) - log(true_var))^2)) %>% 
  group_by(area %in% low_cov_areas) %>% summarise(mean(sq_diff))

var.results %>% filter(model=='Mod1') %>% group_by(area) %>% 
  summarise(v_ratio=mean(v_est/v_truth,na.rm=T)) 


