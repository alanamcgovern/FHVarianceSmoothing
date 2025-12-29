## results script for settings that use same surface for each simulation

library(tidyverse)
library(readr)
library(knitr)

setting <- 2
setwd(paste0('/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Simulations/Sim',setting))

load(file=paste0('params.rda'))
load(file=paste0('objects.rda'))
load(file=paste0('direct.rda'))
load(file=paste0('sampled_clusters.rda'))
adm2.dir.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct, function(x){x$mean[i]})),na.rm = T)})
adm2.dir.variance.est <- sapply(1:n_admin2, function(i){mean(unlist(lapply(direct, function(x){x$var[i]})),na.rm = T)})

load(paste0('admin2.fh.comparisons.rda'))
results.cluster <- read_csv('Summary.csv')

results <- rbind(results.adm2[,c('model','sim','area','mean','lower','upper')],
                 results.cluster[,c('model','sim','area','mean','lower','upper')])
results <- results %>% filter(model!='oracle')

results$true_mean <- params$admin2_mean[results$area]

results <- results %>% mutate(ciwidth = upper - lower,
                              se = (mean-true_mean)^2,
                              cov = lower<true_mean & upper>true_mean,
                              i.score = upper - lower + 20*(lower-true_mean)*I(true_mean<lower) + 20*(true_mean-upper)*I(true_mean>upper))

var.results <- read.csv('Summary_Var.csv')
var.results$true_var <- adm2.dir.variance[var.results$area]

hyper.results <- read.csv('Summary_hyper.csv')

# add degrees of freedom to results

df_df <- lapply(1:length(sampled_clusters),function(i){
  dir.dat <- direct[[i]]
  change_id <- which(dir.dat$variance < 1e-10)
  if(length(change_id)>0)
    dir.dat[change_id,]$mean <- NA
  data_areas <- which(!is.na(dir.dat$mean))
  
  tmp1 <- sampled_clusters[[i]] %>% filter(admin2 %in% data_areas) %>% group_by(admin2) %>% reframe(n=n(),h=length(unique(urban))) %>% 
    mutate(df = pmax(1,n-h)) %>% dplyr::select(admin2,df)
  out1 <- left_join(data.frame(admin2 = 1:n_admin2),tmp1) %>% mutate(sim=i, model = 'Mod1') %>% rename(area = admin2)
  
  tmp2 <- hyper.results %>% filter(str_detect(param,'df'), sim==i, model=='Mod2') %>% rename(df=mean)
  tmp2$admin2 <- data_areas
  out2 <- left_join(data.frame(admin2 = 1:n_admin2),tmp2[,c('admin2','df')]) %>% mutate(sim=i, model = 'Mod2')%>% rename(area = admin2)
  
  return(rbind(out1,out2))
})
df_df <- do.call(rbind,df_df)

results <- left_join(results,df_df,by=c('sim','model','area'))


model_labs <- c(
  Mod1 = "Naive + structured",
  Mod2 = "Satt (no kappa) + structured",
  Mod1a = "Naive + unstructured",
  Mod2a = "Satt (no kappa) + unstructured",
  oracle = "Oracle",
  std = "No variance smoothing")

## plots  ------
pdf(paste0('Sim',setting,' Plots.pdf'))

{
  
  # t <- results %>% group_by(model,area) %>% summarise(mean = mean(mean))
  # 
  # par(mfrow=c(3,2))
  # 
  # plot(params$admin2_mean,t[t$model=='std',]$mean,xlab='True mean',ylab='Average posterior mean')
  # abline(0,1)
  # 
  # plot(params$admin2_mean,t[t$model=='Mod1',]$mean,xlab='True mean',ylab='Average posterior mean')
  # abline(0,1)
  # 
  # plot(params$admin2_mean,t[t$model=='Mod1a',]$mean,xlab='True mean',ylab='Average posterior mean')
  # abline(0,1)
  # 
  # plot(params$admin2_mean,t[t$model=='Mod2',]$mean,xlab='True mean',ylab='Average posterior mean')
  # abline(0,1)
  # 
  # plot(params$admin2_mean,t[t$model=='Mod2a',]$mean,xlab='True mean',ylab='Average posterior mean')
  # abline(0,1)
  
  
  # RMSE
  df_rmse <- results %>% 
    group_by(model, area, true_mean) %>% 
    summarise(rmse = sqrt(mean(se)), .groups = "drop")
  
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
  df_cov <- results %>% #filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
    group_by(model, area, true_mean) %>% 
    summarise(cov = mean(cov), .groups = "drop")
  
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
  df_width <- results %>% 
    group_by(model, area, true_mean) %>% 
    summarise(mean_width = mean(upper - lower), .groups = "drop")
  
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
  df_iscore <- results %>% 
    group_by(model, area, true_mean) %>% 
    summarise(i.score = mean(i.score), .groups = "drop")
  
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
  
  var.results %>% group_by(model,area,true_var) %>% summarise(v_est=mean(mean)) %>% ggplot() +
    geom_point(aes(true_var,v_est)) + xlab('Empirical') + ylab('Average estimate') + facet_wrap(~model,labeller = as_labeller(model_labs)) + 
    geom_abline(intercept = 0,slope=1) + ggtitle('Estimate of V')
  
  par(mfrow=c(1,1))
  plot(results[results$model=='Mod1',]$df,results[results$model=='Mod2',]$df,cex=0.5,xlab='Naive',ylab='Satt (no kappa)',main='Degrees of freedom')
  abline(0,1)
  
}

dev.off()

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



