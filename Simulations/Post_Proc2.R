## results script for settings that use same surface for each simulation

library(tidyverse)
library(readr)
library(knitr)


df_metrics <- var.results <- sig2.results <- NULL
for(setting in 1:3){
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
  
  sig2.results <- rbind(sig2.results,sig2.results_tmp[,c('sim','area','model','mean','true_value','setting')])
  
}

# chisq param estiamtes
# chisq.params <- read.csv('Summary_chisq_params.csv')
# load('correct_chisq_params.rda')
# chisq.params <- left_join(chisq.params,chisq_correct_params,by=c('sim','area'))


model_labs <- c(
  Mod1 = "Naive-struct",
  Mod2 = "Satt-struct",
  Mod1a = "Naive-unstruct",
  Mod2a = "Satt-unstruct",
  oracle = "Oracle",
  std = "Standard")

model_colors <- c(
  # Pair 1 – blue
  "Mod1"   = "#1f78b4",  # dark blue
  "Mod1a"  = "#a6cee3",  # light blue
  
  # Pair 2 – orange
  "Mod2"   = "#6a3d9a",  # dark orange
  "Mod2a"  = "#cab2d6",  # light orange
  
  # Pair 3 – purple
  "oracle" = "#A1D99B",  # dark purple
  "std"    = "#E6AB02"   # light purple
)



## summary ----
df_rmse <- results %>% 
  group_by(model, area, true_mean) %>% 
  summarise(rmse = sqrt(mean(se)), .groups = "drop")
df_rmse %>% group_by(model) %>% summarise(mean(rmse))

df_cov <- results %>% #filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
  group_by(model, area, true_mean) %>% 
  summarise(cov = mean(cov), .groups = "drop")
df_cov %>% group_by(model) %>% summarise(mean(cov))

df_width <- results %>% 
  group_by(model, area, true_mean) %>% 
  summarise(mean_width = mean(upper - lower), .groups = "drop")
df_width%>% group_by(model) %>% summarise(mean(mean_width))

df_iscore <- results %>% 
  group_by(model, area, true_mean) %>% 
  summarise(i.score = mean(i.score), .groups = "drop")
df_iscore %>% group_by(model) %>% summarise(mean(i.score))

## final plots  ------

# 910 x 400

# v1 <- sig2.results %>% group_by(setting,model,area,true_value) %>% summarise(ratio=mean(mean/true_value)) %>% ggplot(aes(model,ratio)) +
#   geom_boxplot(fill='grey90',outlier.size = 0.5) + 
#   stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "blue") + 
#   geom_hline(yintercept = 1,col='red',lty=5) +
#   facet_grid(~setting) + scale_y_continuous(trans='log',breaks = c(1,3,5,10)) +
#   scale_x_discrete(labels=model_labs) + labs(x=NULL,y=NULL)+ theme_bw() +
#   theme(axis.text.x = element_text(size=10,angle = 45, hjust = .75,vjust=0.8),
#         axis.text.y = element_text(size=8))

 var.results %>% group_by(setting,model,area,true_var) %>% summarise(v_ratio=mean(mean/true_var)) %>% ggplot(aes(model,v_ratio)) +
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  geom_hline(yintercept = 1,col='red',lty=5) +
  facet_grid(~setting) + scale_y_continuous(trans='log',breaks = c(0.1,0.25,0.5,1,2.5,5,10)) +
  labs(x=NULL,y=NULL)+ theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8),
        legend.position = 'bottom')

#ggarrange(plotlist = list(v1,v2),nrow = 2)

# 910 x 1050

d1 <- df_metrics %>% ggplot(aes(model,rmse)) + ggtitle('A. RMSE')+ 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks=c(0.05,0.1,0.25,0.5,1)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  facet_grid(~setting) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8))

d2 <- df_metrics %>% ggplot(aes(model,100*cov)) +  ggtitle('B. Coverage of 90% credible interval') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  geom_hline(yintercept = 90,col='red',lty=5) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL) + theme_bw() +
  facet_grid(~setting) + 
  theme(axis.text.x =  element_blank(),
        axis.text.y = element_text(size=8))

d3 <- df_metrics %>% ggplot(aes(model,mean_width)) + ggtitle('C. Average 90% credible interval width')+
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  labs(x=NULL,y=NULL) + theme_bw() +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  facet_grid(~setting) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8))

d4 <- df_metrics %>% ggplot(aes(model,i.score)) + ggtitle('D. Average interval score') + 
  geom_boxplot(aes(fill=model),outlier.size = 0.15, median.linewidth = 0.5,box.linewidth = 0.25,whisker.linewidth = 0.25) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.9, lwd=0.25, color = "grey60") + 
  scale_y_continuous(trans='log',breaks = c(0.5,1,3,7)) +
  scale_fill_manual(name="Model",values = model_colors,labels=model_labs) + 
  labs(x=NULL,y=NULL)+ theme_bw() +
  facet_grid(~setting) + 
  theme(axis.text.x =element_blank(),
        axis.text.y = element_text(size=8))

ggarrange(plotlist = list(d1,d2,d3,d4),ncol=1,common.legend = T,legend = "bottom")


#pdf(paste0('Sim',setting,' Plots.pdf'))

### scratch -----

# compare estimates across models

# mean
results %>% group_by(model,area,true_mean) %>% summarise(val = mean(mean)) %>%
  ggplot()+ geom_point(aes(true_mean,val)) + facet_wrap(~ model,nrow=3,labeller = as_labeller(model_labs)) +
  geom_abline(intercept = 0,slope = 1)


logsig2.results %>% group_by(model,area,true_value) %>% summarise(est=mean(mean)) %>% ggplot() +
  geom_point(aes(true_value,est)) + ylab('Average estimate') + facet_wrap(~model,labeller = as_labeller(model_labs)) + 
  geom_abline(intercept = 0,slope=1) + ggtitle('Estimate of log(sig2)')

chisq.params %>% filter(model %in% c('Mod1','Mod2')) %>% dplyr::select(sim,area,model,df) %>% 
  pivot_wider(id_cols = sim:area,names_from = model,values_from = df) %>%
  ggplot() + geom_point(aes(Mod1,Mod2)) + geom_abline(intercept = 0,slope = 1,col='red')

chisq.params %>% mutate(t_mean = scale*df) %>% filter(model %in% c('Mod1','Mod2')) %>% dplyr::select(sim,area,model,t_mean) %>% 
  pivot_wider(id_cols = sim:area,names_from = model,values_from = t_mean) %>%
  ggplot() + geom_point(aes(Mod1,Mod2),size=0.5) + geom_abline(intercept = 0,slope = 1,col='red')

chisq.params %>% group_by(model,area,dof_correct) %>% summarise(val = mean(df)) %>%
  ggplot()+ geom_point(aes(dof_correct,val),size=0.5) + facet_wrap(~ model,nrow=3,labeller = as_labeller(model_labs)) +
  geom_abline(intercept = 0,slope = 1)


  # RMSE
  df_label <- df_rmse %>%
    group_by(model) %>% 
    summarise(mean_rmse = mean(rmse))
  
  g <- df_rmse %>% ggplot(aes(true_mean, rmse)) + geom_point() +
    facet_wrap(~ model,nrow=3,labeller = as_labeller(model_labs)) +
    geom_label(data = df_label,
               aes( x = -Inf, y = Inf, label = paste0("RMSE = ", round(mean_rmse, 3))),
               hjust = -0.1, vjust = 1.1) +
    ggtitle("RMSE")
  
  print(g)
  
  
  # Coverage
  df_label <- df_cov %>% #filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
    group_by(model) %>% 
    summarise(cov = mean(cov))
  
  g <- df_cov %>% ggplot(aes(true_mean, cov)) + geom_point() +
    facet_wrap(~ model, nrow=3,labeller = as_labeller(model_labs)) +
    geom_hline(yintercept = 0.9,col='red') +
    geom_label(data = df_label,
               aes( x = -Inf, y = Inf, label = paste0("Coverage = ", round(cov, 3))),
               hjust = -0.1, vjust = 1.1) +
    ggtitle("Coverage")
  
  print(g)
  
  # Interval widths
  df_label <- df_width %>% 
    group_by(model) %>% 
    summarise(mean_width = mean(mean_width))
  
  g <- df_width %>% ggplot(aes(true_mean, mean_width)) + geom_point() +
    facet_wrap(~ model,labeller = as_labeller(model_labs)) +
    geom_label(data = df_label,
               aes( x = -Inf, y = Inf, label = paste0("Mean CI width = ", round(mean_width, 3))),
               hjust = -0.1, vjust = 1.1) +
    ggtitle("CI Width")
  print(g)
  
  # Interval score
  df_label <- df_iscore %>% 
    group_by(model) %>% 
    summarise(ais = mean(i.score))
  
  g <- df_iscore %>% ggplot(aes(true_mean, i.score)) + geom_point() +
    facet_wrap(~ model,labeller = as_labeller(model_labs)) +
    geom_label(data = df_label,
               aes( x = -Inf, y = Inf, label = paste0("AIS = ", round(ais, 3))),
               hjust = -0.1, vjust = 1.1) +
    ggtitle("Average Interval Scores")
  
  print(g)
  

  par(mfrow=c(1,1))
  plot(results[results$model=='Mod1',]$df,results[results$model=='Mod2',]$df,cex=0.5,xlab='Naive',ylab='Satt (no kappa)',main='Degrees of freedom')
  abline(0,1)
  

hyper.results %>% filter(param == 'sig_u[1]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'sig_u[2]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'phi2') %>% group_by(model) %>% summarise(mean(mean))

hyper.results %>% filter(param == 'beta[1]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'beta[2]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'beta[3]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'beta[4]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'beta[5]') %>% group_by(model) %>% summarise(mean(mean))

hyper.results %>% filter(param == 'gamma[1]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'gamma[2]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'gamma[3]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'gamma[4]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param == 'gamma[5]') %>% group_by(model) %>% summarise(mean(mean))

logsig2.results %>% group_by(model,area,true_value) %>% summarise(est=mean(mean)) %>% ggplot() +
  geom_point(aes(true_value,est)) + ylab('Average estimate') + facet_wrap(~model,labeller = as_labeller(model_labs)) + 
  geom_abline(intercept = 0,slope=1) + ggtitle('Estimate of log(sig2)')



results %>% filter(model %in% c('Mod1','Mod2')) %>% dplyr::select(sim,area,model,ciwidth) %>% 
  pivot_wider(id_cols = sim:area,names_from = model,values_from = ciwidth) %>%
  ggplot() + geom_point(aes(Mod1,Mod2)) + geom_abline(intercept = 0,slope = 1,col='red')

var.results %>% filter(model %in% c('Mod1','Mod2')) %>% dplyr::select(sim,area,model,mean) %>% 
  pivot_wider(id_cols = sim:area,names_from = model,values_from = mean) %>%
  ggplot() + geom_point(aes(Mod1,Mod2)) + geom_abline(intercept = 0,slope = 1,col='red')







