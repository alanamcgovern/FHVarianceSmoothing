## results script for settings that use same surface for each simulation

library(tidyverse)
library(readr)
library(knitr)
library(ggpubr)
library(ggh4x)

######### SUMMARY OF DIAGNOSTICS --------------------
diag_full <- NULL
for(setting in 1:5){
  setwd(paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/Sim',setting))
  tmp <- read_csv('Summary_diag.csv')
  tmp$setting <- setting
  diag_full <- rbind(diag_full,tmp)
}

diag_full %>% filter(!(model %in% c('Mod4','Mod4a'))) %>% group_by(model) %>% summarise(sum(max_treedepth))

diag_full %>% filter(!(model %in% c('Mod4','Mod4a'))) %>% filter(max_treedepth>0)

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

dist_metrics_full <- NULL
for(setting in c(1:5)){
  setwd(paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/Sim',setting))
  load('approx_dist_metrics.rda') # from Comparison_Models.R
  dist_metrics_full <- rbind(dist_metrics_full,dist_metrics)
}

model_labs <- c(n = "Simple", sw = "SASW")
model_colors <- c("n" = "#1f78b4", "sw" = "#6a3d9a")

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

dist_metrics_full <- merge(setting_key,dist_metrics_full)

m1 <- dist_metrics_full %>% group_by(setting_labels,area) %>% 
  summarise(n = mean(Wass_Naive),sw = mean(Wass_Satt)) %>% 
  pivot_longer(cols = c(n,sw),names_to = 'model',values_to = 'mean_wass') %>%
  ggplot(aes(model,mean_wass)) + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
 # scale_y_continuous(trans='log',breaks = c(5e-4,5e-3,0.05,0.2)) +
  scale_fill_manual(name="Distribution",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  facet_wrap(~factor(setting_labels,
                     levels = setting_levels),nrow=1) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8)) + 
  ggtitle('A. Average Wasserstein distance of order 2')

m2 <- dist_metrics_full %>% group_by(setting_labels,area) %>% 
  summarise(n = mean(Ratio_Naive),sw = mean(Ratio_Satt)) %>% 
  pivot_longer(cols = c(n,sw),names_to = 'model',values_to = 'mean_ratio') %>%
  ggplot(aes(model,mean_ratio)) + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
 # scale_y_continuous(trans='log',breaks=c(0.5,1,2.5,5,10),limits = c(0,2)) + 
  ylim(c(0,2)) +
  geom_hline(yintercept = 1,col='red')+ 
  scale_fill_manual(name="Distribution",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  facet_wrap(~factor(setting_labels,
                     levels = setting_levels),nrow=1) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8)) +
  ggtitle('B. Average ratio of means')

m3 <- dist_metrics_full %>% group_by(setting_labels,area) %>% 
  summarise(n = mean(Diff_Naive),sw = mean(Diff_Satt)) %>% 
  pivot_longer(cols = c(n,sw),names_to = 'model',values_to = 'mean_diff') %>%
  ggplot(aes(model,mean_diff)) + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  geom_hline(yintercept = 0,col='red')+ 
  scale_fill_manual(name="Distribution",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  facet_wrap(~factor(setting_labels,
                     levels = setting_levels),nrow=1,scales='free_y') + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8)) +
  ggtitle('B. Average difference in means')

# 940 x 505
ggarrange(plotlist = list(m1,m3), nrow=2, common.legend = T, legend='bottom')


############ MODEL PERFORMANCE -------------------------------

df_metrics <- var.results <- sig2.results <- NULL
for(setting in c(1:7)){
  setwd(paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/Sim',setting))
  
  load(file=paste0('params.rda'))
  load(file=paste0('objects.rda'))
  load(file=paste0('direct.rda'))
  load(file=paste0('sampled_clusters.rda'))
  n_admin2 <- length(params$admin2_mean)
  adm2.dir.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct, function(x){x$mean[i]})),na.rm = T)})
  
  # mean estimates
  load(paste0('admin2.fh.comparisons.rda'))
  results.cluster <- read_csv('Summary.csv')
  
  results <- rbind(results.adm2[,c('model','sim','area','mean','lower','upper')],
                   results.cluster[,c('model','sim','area','mean','lower','upper')])
  
  results$true_mean <- params$admin2_mean[results$area]
  
  results <- results %>% mutate(ciwidth = upper - lower,
                                se = (mean-true_mean)^2,
                                cov = lower<true_mean & upper>true_mean,
                                i.score = upper - lower + 20*(lower-true_mean)*I(true_mean<lower) + 20*(true_mean-upper)*I(true_mean>upper))
  
  df_rmse <- results %>% 
    group_by(model, area, true_mean) %>% 
    summarise(rmse = sqrt(mean(se)), .groups = "drop")
  
  df_cov <- results %>% 
    group_by(model, area, true_mean) %>% 
    summarise(cov = mean(cov), .groups = "drop")
  
  df_width <- results %>% 
    group_by(model, area, true_mean) %>% 
    summarise(mean_width = mean(upper - lower), .groups = "drop")
  
  df_iscore <- results %>% 
    group_by(model, area, true_mean) %>% 
    summarise(i.score = mean(i.score), .groups = "drop")
  
  df_metrics_tmp <- merge(merge(df_rmse,df_cov),merge(df_width,df_iscore))
  df_metrics_tmp$setting <- setting
  df_metrics <- rbind(df_metrics,df_metrics_tmp)
  
  # design variance estimates
  var.results_tmp <- read.csv('Summary_Var.csv')
  standard_est <- lapply(1:100,function(k){
    tmp <- direct[[k]][,c('admin2','variance')] %>% rename(mean=variance,area=admin2)
    tmp$sim <- k
    return(tmp)
  })
  standard_est <- do.call('rbind',standard_est)
  standard_est <- standard_est[!is.na(standard_est$mean),]
  standard_est <- standard_est[standard_est$mean > 1e-10,]
  standard_est$model <- 'std'
  var.results_tmp <- rbind(var.results_tmp,standard_est[,c('sim','area','model','mean')])
  var.results_tmp$true_var <- adm2.dir.variance[var.results_tmp$area]
  var.results_tmp$setting <- setting
  
  var.results <- rbind(var.results,var.results_tmp)
  
  hyper.results <- read.csv('Summary_hyper.csv')
  sig2.results_tmp <- hyper.results %>% filter(str_detect(hyper.results$param,'sig2')) %>% mutate(mean = exp(mean))
  sig2.results_tmp$area <- as.numeric((gsub("log_sig2\\[|\\]", "",sig2.results_tmp$X)))
  sig2.results_tmp$true_value <- params$sig2_pop[sig2.results_tmp$area]
  sig2.results_tmp$setting <- setting
  
  load('var_constants.rda')
  sig2.results_tmp <- left_join(sig2.results_tmp,var_constants,by=c('sim','area'))
  sig2.results_tmp <- sig2.results_tmp %>% mutate(cons = ifelse(model %in% c('Mod1','Mod1a'),cons_naive,cons_exact))
  
  sig2.results <- rbind(sig2.results,sig2.results_tmp[,c('sim','area','model','mean','true_value','cons','setting')])
  
}

# chisq param estiamtes
# chisq.params <- read.csv('Summary_chisq_params.csv')
# load('correct_chisq_params.rda')
# chisq.params <- left_join(chisq.params,chisq_correct_params,by=c('sim','area'))

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
sig2.results <- merge(sig2.results,setting_key)


model_labs <- c(
  Mod1 = "Simple-struct",
  Mod2 = "SASW-struct",
  Mod4 = "Debias-SASW-struct",
  Mod1a = "Simple-unstruct",
  Mod2a = "SASW-unstruct",
  Mod4a = "Debias-SASW-unstruct",
  oracle = "Oracle",
  std = "Standard")

model_colors <- c(
  "Mod1"   = "#1f78b4",  # dark blue
  "Mod1a"  = "#a6cee3",  # light blue
  
  "Mod2"   = "#6a3d9a",  # dark orange
  "Mod2a"  = "#cab2d6",  # light orange
  
  "Mod4"   = "#e31a1c",
  "Mod4a"  = "#fb9a99",
  
  "oracle" = "#A1D99B",  # dark purple
  "std"    = "#E6AB02"   # light purple
)



## summary ----
df_metrics %>% group_by(setting,model) %>% summarise(mean(rmse),
                                                     mean(cov),
                                                     mean(mean_width),
                                                     mean(i.score))

#df_metrics <- df_metrics %>% filter(model %in% c('Mod1','Mod2','Mod1a','Mod2a','oracle','std'))

## plots  ------

# 910 x 400
common_theme <- theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8),
        legend.position = 'bottom') 
  
var.results <- var.results %>% filter(setting!=6)
sig2.results <- sig2.results %>% filter(setting!=6)


tmp <- sig2.results %>% mutate(var_model = cons*true_value) %>% 
  dplyr::select(setting_labels,sim,area,model,var_model)

# compare theoretical to true design variance (Fig 3) ------
# 940 x 350
merge(tmp,var.results) %>% filter(setting %in% 1:5,(model %in% c('Mod1','Mod2'))) %>% 
  group_by(setting_labels,model,area) %>% 
  summarise(val = mean(var_model/true_var))  %>% 
  ggplot(aes(model,val)) +
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_fill_manual(name="Model",values = model_colors,labels=c(
    Mod1 = "Simple",
    Mod2 = "SASW")) + 
  geom_hline(yintercept = 1,col='red',lty=5) +
  facet_wrap(~factor(setting_labels,
                     levels = setting_levels),nrow=1) +
  scale_y_continuous(trans='log',breaks = c(0.1,0.25,0.5,1,2.5,5,10)) +
  labs(x=NULL,y=NULL)+ common_theme

# compare estiamted to true design variance (Fig 4) ---- 

v1 <- sig2.results %>% filter(setting %in% 1:5,!(model %in% c('Mod4','Mod4a'))) %>% group_by(setting_labels,model,area) %>% 
  summarise(val = mean(mean/true_value))  %>% 
  ggplot(aes(model,val)) +
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  geom_hline(yintercept = 1,col='red',lty=5) +
  facet_wrap(~factor(setting_labels,
                     levels = setting_levels),nrow=1) +
  scale_y_continuous(trans='log',breaks = c(0.1,0.25,0.5,1,2.5,5,10)) +
  labs(x=NULL,y=NULL)+ common_theme +
  theme(legend.position = 'none') +
  ggtitle('A. Within-stratum superpopulation variance')

v2 <-  var.results %>% filter(setting %in% 1:5,!(model %in% c('Mod4','Mod4a')) )%>% group_by(setting_labels,model,area,true_var) %>% 
   summarise(v_ratio=mean(mean/true_var)) %>% ggplot(aes(model,v_ratio)) +
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  geom_hline(yintercept = 1,col='red',lty=5) +
   facet_wrap(~factor(setting_labels,
                      levels = setting_levels),nrow=1) + 
   scale_y_continuous(trans='log',breaks = c(0.1,0.25,0.5,1,2.5,5,10),limits = c(0.05,14)) +
  labs(x=NULL,y=NULL)+ common_theme + 
   ggtitle('B. Design variance')

# 830 x 515
ggarrange(plotlist = list(v1,v2),nrow = 2,heights = c(1,1.25))

# compare mean and interval estimates (Fig 5) ---------

d1 <- df_metrics%>% filter(setting %in% 1:5,!(model %in% c('Mod4','Mod4a'))) %>% ggplot(aes(model,rmse)) + ggtitle('A. RMSE')+ 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks=c(0.05,0.1,0.25,0.5,1)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  facet_wrap(~factor(setting_labels,
                    levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),axis.text.x = element_blank(),
        axis.text.y = element_text(size=6))

d2 <- df_metrics%>% filter(setting %in% 1:5,!(model %in% c('Mod4','Mod4a'))) %>% ggplot(aes(model,100*cov)) +  
  ggtitle('B. Coverage of 90% credible intervals (clipped at 50%)') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  geom_hline(yintercept = 90,col='red',lty=5) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  ylim(c(50,100)) +
  facet_wrap(~factor(setting_labels,
                    levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),axis.text.x =  element_blank(),
        axis.text.y = element_text(size=6))

d3 <- df_metrics%>% filter(setting %in% 1:5,!(model %in% c('Mod4','Mod4a'))) %>% ggplot(aes(model,mean_width)) + ggtitle('C. Average 90% credible interval width')+
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  labs(x=NULL,y=NULL) + theme_bw() +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  facet_wrap(~factor(setting_labels,
                    levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),
    axis.text.x = element_blank(),
        axis.text.y = element_text(size=6))

d4 <- df_metrics%>% filter(setting %in% 1:5,!(model %in% c('Mod4','Mod4a'))) %>% ggplot(aes(model,i.score)) + ggtitle('D. Average interval score for 90% credible intervals') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks = c(0.25,0.5,1,2,4)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL)+ theme_bw() +
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
setting_levels <- c('1: listed cluster sizes are correct',
                    '5: cluster sizes are re-enumerated')

df_metrics <- merge(df_metrics,setting_key)


a1 <- df_metrics%>% filter(setting %in% c(1,7),!(model %in% c('Mod4','Mod4a'))) %>% ggplot(aes(model,rmse)) + ggtitle('A. RMSE')+ 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks=c(0.05,0.1,0.25,0.5,1)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  facet_wrap(~factor(setting_labels_new,
                     levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),axis.text.x = element_blank(),
        axis.text.y = element_text(size=6))

a2 <- df_metrics%>% filter(setting %in% c(1,7),!(model %in% c('Mod4','Mod4a'))) %>% ggplot(aes(model,100*cov)) +  
  ggtitle('B. Coverage of 90% credible intervals (clipped at 50%)') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  geom_hline(yintercept = 90,col='red',lty=5) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  ylim(c(50,100)) +
  facet_wrap(~factor(setting_labels_new,
                     levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),axis.text.x =  element_blank(),
        axis.text.y = element_text(size=6))

a3 <- df_metrics%>% filter(setting %in% c(1,7),!(model %in% c('Mod4','Mod4a'))) %>% ggplot(aes(model,mean_width)) + ggtitle('C. Average 90% credible interval width')+
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  labs(x=NULL,y=NULL) + theme_bw() +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  facet_wrap(~factor(setting_labels_new,
                     levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6))

a4 <- df_metrics%>% filter(setting %in% c(1,7),!(model %in% c('Mod4','Mod4a'))) %>% ggplot(aes(model,i.score)) + ggtitle('D. Average interval score for 90% credible intervals') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks = c(0.25,0.5,1,2,4)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL)+ theme_bw() +
  facet_wrap(~factor(setting_labels_new,
                     levels = setting_levels),nrow=1) + 
  theme(plot.title = element_text(size=12),axis.text.x =element_blank(),
        axis.text.y = element_text(size=6))

#850 x 550
ggarrange(plotlist = list(a1,a2,a3,a4),nrow=2,ncol=2,common.legend = T,legend = "bottom")




