## results script for settings that use different surface for each simulation

library(tidyverse)
library(readr)
library(knitr)

setting <- 1
setwd(paste0('/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Simulations/Sim',setting))

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

par(mfrow=c(1,1))
plot(results[results$model=='Mod1',]$df,results[results$model=='Mod2',]$df,cex=0.5,xlab='Naive',ylab='Satt (no kappa)',main='Degrees of freedom')
abline(0,1)


}
dev.off()

