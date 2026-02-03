## results script for settings that use different surface for each simulation

library(tidyverse)
library(readr)
library(knitr)

setting <- 13
setwd(paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/Sim',setting))

load(file=paste0('true_vals.rda'))
load(file=paste0('objects.rda'))
load(file=paste0('direct.rda'))
load(file=paste0('sampled_clusters.rda'))

load(paste0('admin2.fh.comparisons.rda'))
results.cluster <- read_csv('Summary.csv')

results <- rbind(results.adm2[,c('model','sim','area','mean','lower','upper')],
                  results.cluster[,c('model','sim','area','mean','lower','upper')])
results <- results %>% filter(model!='oracle')

results$true_mean <- NA
for(j in 1:100){
  results[results$sim==j,]$true_mean <- true_vals[[j]]$admin2_mean[results[results$sim==j,]$area]
}

results <- results %>% mutate(ciwidth = upper - lower,
                              se = (mean-true_mean)^2,
                              cov = lower<true_mean & upper>true_mean,
                              i.score = upper - lower + 20*(lower-true_mean)*I(true_mean<lower) + 20*(true_mean-upper)*I(true_mean>upper))


var.results <- read.csv('Summary_Var.csv')

hyper.results <- read.csv('Summary_hyper.csv')

logsig2.results <- hyper.results %>% filter(str_detect(hyper.results$param,'sig2'))
logsig2.results$area <- as.numeric((gsub("log_sig2\\[|\\]", "",logsig2.results$X)))

logsig2.results$true_value <- NA
for(j in 1:100){
  # ONLY WORKS WHEN VARIANCE IS THE SAME FOR URBAN/RURAL
  logsig2.results[logsig2.results$sim==j,]$true_value <- log(true_vals[[j]]$sig2[logsig2.results[logsig2.results$sim==j,]$area])
}

# add degrees of freedom to results
chisq.params <- read.csv('Summary_chisq_params.csv')
load('correct_chisq_params.rda')
chisq.params <- left_join(chisq.params,chisq_correct_params,by=c('sim','area'))

model_labs <- c(
  Mod1 = "Naive + structured",
  Mod2 = "Satt (no kappa) + structured",
  Mod1a = "Naive + unstructured",
  Mod2a = "Satt (no kappa) + unstructured",
  Mod3 = "Satt (w kappa) + structured",
  oracle = "Oracle",
  std = "No variance smoothing")

## plots  --------

pdf(paste0('Sim',setting,'_Plots.pdf'))
{

gridExtra::grid.table(results %>% group_by(model) %>% reframe(rmse = sqrt(mean(se)),coverage = mean(cov), ais = mean(i.score), width = mean(ciwidth)))

# RMSE
df_rmse <- results %>% 
  group_by(model,sim) %>% 
  summarise(rmse = sqrt(mean(se)), .groups = "drop")

par(mfrow=c(3,2))

plot(df_rmse[df_rmse$model=='std',]$rmse, df_rmse[df_rmse$model=='Mod1',]$rmse,xlab='No smoothing',ylab='Naive + structured',main='RMSE (per simulation)')
abline(0,1)

plot(df_rmse[df_rmse$model=='std',]$rmse, df_rmse[df_rmse$model=='Mod2',]$rmse,xlab='No smoothing',ylab='Satt (no kappa) + structured')
abline(0,1)

plot(df_rmse[df_rmse$model=='Mod1',]$rmse, df_rmse[df_rmse$model=='Mod2',]$rmse,xlab='Naive + structured',ylab='Satt (no kappa) + structured')
abline(0,1)

plot(df_rmse[df_rmse$model=='Mod1',]$rmse, df_rmse[df_rmse$model=='Mod1a',]$rmse,xlab='Naive + structured',ylab='Naive + unstructured')
abline(0,1)

plot(df_rmse[df_rmse$model=='Mod2',]$rmse, df_rmse[df_rmse$model=='Mod2a',]$rmse,xlab='Satt (no kappa) + structured',ylab='Satt (no kappa) + unstructured')
abline(0,1)

# Coverage

df_cov <- results %>% 
  group_by(model, sim) %>% 
  summarise(cov = mean(cov), .groups = "drop")

df_label <- df_cov %>% 
  group_by(model) %>% 
  summarise(cov = mean(cov))

g <- df_cov %>% ggplot(aes(cov)) + geom_histogram() +
  facet_wrap(~ model,nrow = 3,labeller = as_labeller(model_labs)) +
  geom_vline(xintercept = 0.9,col='red') +
  geom_label(data = df_label,
             aes( x = -Inf, y = Inf, label = paste0("Average = ", round(cov, 3))),
             hjust = -0.1, vjust = 1.1) +
  ggtitle("Coverage")

print(g)

# CI width
df_width <- results %>% 
  group_by(model,sim) %>% 
  summarise(mean_width = mean(ciwidth), .groups = "drop")

df_label <- df_width %>% 
  group_by(model) %>% 
  summarise(mean = mean(mean_width))

g <- df_width %>% ggplot(aes(mean_width)) + geom_histogram() +
  facet_wrap(~ model,nrow = 3,labeller = as_labeller(model_labs)) +
  geom_label(data = df_label,
             aes( x = -Inf, y = Inf, label = paste0("Average = ", round(mean, 3))),
             hjust = -0.1, vjust = 1.1) +
  ggtitle("Interval width")

print(g)

# AIS
df_ais <- results %>% 
  group_by(model,sim) %>% 
  summarise(ais = mean(i.score), .groups = "drop")

par(mfrow=c(3,2))

plot(df_ais[df_ais$model=='std',]$ais, df_ais[df_ais$model=='Mod1',]$ais,xlab='No smoothing',ylab='Naive + structured',main='AIS (per simulation)')
abline(0,1)

plot(df_ais[df_ais$model=='std',]$ais, df_ais[df_ais$model=='Mod2',]$ais,xlab='No smoothing',ylab='Satt (no kappa) + structured')
abline(0,1)

plot(df_ais[df_ais$model=='Mod1',]$ais, df_ais[df_ais$model=='Mod2',]$ais,xlab='Naive + structured',ylab='Satt (no kappa) + structured')
abline(0,1)

plot(df_ais[df_ais$model=='Mod1',]$ais, df_ais[df_ais$model=='Mod1a',]$ais,xlab='Naive + structured',ylab='Naive + unstructured')
abline(0,1)

plot(df_ais[df_ais$model=='Mod2',]$ais, df_ais[df_ais$model=='Mod2a',]$ais,xlab='Satt (no kappa) + structured',ylab='Satt (no kappa) + unstructured')
abline(0,1)


# degrees of freedom
chisq.params %>% ggplot() + geom_point(aes(dof_correct,df)) + facet_grid(~model) + geom_abline(intercept = 0,slope=1,col='red')

chisq.params %>% filter(model %in% c('Mod1','Mod2')) %>% dplyr:select(sim,area,model,df) %>%
  pivot_wider(id_cols=c(sim,area,dof_correct,scale_correct),values_from=scale:df,names_from=model) %>%
  ggplot() + geom_point(aes(df_Mod1,df_Mod2)) + geom_abline(intercept = 0, slope=1)


}
dev.off()


## compare estimation of log(sig2)
logsig2.results %>% ggplot() + geom_point(aes(true_value,mean),size=0.25) + facet_grid(~model) + geom_abline(intercept = 0,slope=1,col='red')

## compare how well difference in means is estimated
hyper.results %>% filter(param=='beta[2]') %>% ggplot() + geom_histogram(aes(mean)) + facet_grid(~model) + geom_vline(xintercept = 1,col='red')

## compare chi square params
chisq.params %>% filter(model %in% c('Mod1','Mod2')) %>% 
  pivot_wider(id_cols=c(sim,area,dof_correct,scale_correct),values_from=scale:df,names_from=model) %>%
  ggplot() + geom_point(aes(df_Mod1,df_Mod2)) + geom_abline(intercept = 0, slope=1)

chisq.params %>% filter(model %in% c('Mod1','Mod2')) %>% 
  pivot_wider(id_cols=c(sim,area,dof_correct,scale_correct),values_from=scale:df,names_from=model) %>%
  summarise(mean(df_Mod1>df_Mod2,na.rm=T))

chisq.params %>% ggplot() + geom_point(aes(dof_correct,df)) + facet_grid(~model) + geom_abline(intercept = 0,slope=1,col='red')
chisq.params %>% filter(model %in% c('Mod1','Mod2','Mod1.5')) %>% ggplot() + geom_point(aes(scale_correct,scale)) + facet_grid(~model) + geom_abline(intercept = 0,slope=1,col='red') 

chisq.params %>% ggplot() + geom_point(aes(dof_correct*scale_correct,df*scale)) + facet_grid(~model) + geom_abline(intercept = 0,slope=1,col='red')

chisq.params %>% group_by(model) %>% summarise(mean(df*scale - dof_correct*scale_correct,na.rm=T))
chisq.params %>% group_by(model) %>% summarise(mean(2*df*scale^2 - 2*dof_correct*scale_correct^2,na.rm=T))


chisq.params %>% filter(model%in%c('Mod1','Mod2')) %>% ggplot() + 
  geom_point(aes(dof_correct,df)) + facet_grid(~model) + geom_abline(intercept = 0,slope=1,col='red')

plot(chisq.params[chisq.params$model=='Mod1.5',]$df,chisq.params[chisq.params$model=='Mod1',]$df)
abline(0,1)

chisq.params %>% ggplot() + geom_point(aes(2*dof_correct*scale_correct^2,2*df*scale^2)) + facet_grid(~model) + geom_abline(intercept = 0,slope=1,col='red')

merge(var.results,chisq.params) %>% ggplot() + geom_point(aes(df*scale,mean)) + facet_grid(~model) + geom_abline(intercept = 0,slope=1)
