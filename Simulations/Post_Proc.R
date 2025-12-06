library(tidyverse)

setting <- 2
setwd(paste0('/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Simulations/Sim',setting))

#load(file=paste0('params.rda'))
load(file=paste0('true_vals.rda'))
load(file=paste0('objects.rda'))
load(file=paste0('direct.rda'))
load(file=paste0('popdata.rda'))

#load(paste0('admin2.fh.comparisons.rda'))
results.cluster <- read_csv('Summary.csv')

 # results <- rbind(results.adm2[,c('model','sim','area','mean','lower','upper')],
 #                  results.cluster[,c('model','sim','area','mean','lower','upper')])
results <- results.cluster[,c('model','sim','area','mean','lower','upper')]

#results$true_mean <- params$admin2_mean[results$area]

results$true_mean <- NA
for(j in 1:100){
  results[results$sim==j,]$true_mean <- true_vals[[j]]$admin2_mean[results[results$sim==j,]$area]
}


results <- results %>% mutate(ciwidth = upper - lower,
                              se = (mean-true_mean)^2,
                              cov = lower<true_mean & upper>true_mean,
                              i.score = upper - lower + 20*(lower-true_mean)*I(true_mean<lower) + 20*(true_mean-upper)*I(true_mean>upper))


#sig2.results.adm1 <- read_csv('Summary_sig2_admin1_level.csv')
#sig2.results.adm2 <- read_csv('Summary_sig2_admin2_level.csv')

#hyper.results <- read_csv('Summary_hyper.csv')


# 
# model_labs <- c(
#   Mod1 = "Complex (1 kappa)",
#   Mod3 = "Complex (multiple kappa)",
#   Mod5 = "Naive",
#   Mod7 = "Hybrid (1 kappa)",
#   Mod9 = "Hybrid (multiple kappa)",
#   std = 'No smoothing',
#   oracle = 'Oracle'
# )

## plots --------

results %>% group_by(model) %>% reframe(coverage = mean(cov), mis = median(i.score), ais = mean(i.score), width = mean(ciwidth))


## only makes sense when area mean is the same for every simulation
pdf(paste0('Sim',setting,' Plots.pdf'))

{

# RMSE
df_rmse <- results %>% filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
  group_by(model, area, true_mean) %>% 
  summarise(rmse = sqrt(mean(se)), .groups = "drop")

df_label <- df_rmse %>% filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
  group_by(model) %>% 
  summarise(mean_rmse = mean(rmse))

g <- df_rmse %>% ggplot(aes(true_mean, rmse)) + geom_point() +
  facet_wrap(~ model) +
  geom_label(data = df_label,
             aes( x = -Inf, y = Inf, label = paste0("RMSE = ", round(mean_rmse, 3))),
             hjust = -0.1, vjust = 1.1) +
  ggtitle("RMSE")

print(g)


# Coverage
df_cov <- results %>% filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
  group_by(model, area, true_mean) %>% 
  summarise(cov = mean(cov), .groups = "drop")

df_label <- df_cov %>% filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
  group_by(model) %>% 
  summarise(cov = mean(cov))

g <- df_cov %>% ggplot(aes(true_mean, cov)) + geom_point() +
  facet_wrap(~ model) +
  geom_hline(yintercept = 0.9,col='red') +
  geom_label(data = df_label,
             aes( x = -Inf, y = Inf, label = paste0("Coverage = ", round(cov, 3))),
             hjust = -0.1, vjust = 1.1) +
  ggtitle("Coverage")

print(g)

# Interval widths
df_width <- results %>% filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
  group_by(model, area, true_mean) %>% 
  summarise(mean_width = mean(upper - lower), .groups = "drop")

df_label <- df_width %>% filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
  group_by(model) %>% 
  summarise(mean_width = mean(mean_width))

g <- df_width %>% ggplot(aes(true_mean, mean_width)) + geom_point() +
  facet_wrap(~ model) +
  geom_label(data = df_label,
             aes( x = -Inf, y = Inf, label = paste0("Mean CI width = ", round(mean_width, 3))),
             hjust = -0.1, vjust = 1.1) +
  ggtitle("CI Width")
print(g)

# directly compare width between two methods

# plot(df_width[df_width$model=='Naive4',]$mean_width,df_width[df_width$model=='Strat4',]$mean_width)
# abline(0,1)
# 
# plot(df_width[df_width$model=='Naive1',]$mean_width,df_width[df_width$model=='Strat1',]$mean_width)
# abline(0,1)
# 
# plot(df_width[df_width$model=='Strat1',]$mean_width,df_width[df_width$model=='Strat4',]$mean_width)
# abline(0,1)

# Interval score
df_iscore <- results %>% filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
  group_by(model, area, true_mean) %>% 
  summarise(i.score = median(i.score), .groups = "drop")

df_label <- df_iscore %>% filter(!(model %in% c('Mod1_1','Mod2_1','Mod3_1'))) %>%
  group_by(model) %>% 
  summarise(ais = median(i.score))

g <- df_iscore %>% ggplot(aes(true_mean, i.score)) + geom_point() +
  facet_wrap(~ model) +
  geom_label(data = df_label,
             aes( x = -Inf, y = Inf, label = paste0("MIS = ", round(ais, 3))),
             hjust = -0.1, vjust = 1.1) +
  ggtitle("Median Interval Scores")

print(g)


## hyperpar plots ------

par(mfrow=c(2,2))

# ## MOD1_1
# t1 <- sig2.results.adm1 %>% filter(model == 'Mod1_1') %>% group_by(model,area) %>% summarise(est = mean(mean))
# t2 <- sim_pop_long %>% group_by(admin1,cluster) %>% summarise(mean = mean(value)) %>% group_by(admin1) %>% summarise(val = var(mean))
# plot(log(t2$val),log(t1$est),main='Mod1_1')
#abline(0,1)


## MOD1_2
t1 <- sig2.results.adm2 %>% filter(model == 'Mod1_2') %>% group_by(model,area) %>% summarise(est = mean(mean))
t2 <- sim_pop_long %>% group_by(admin2,cluster) %>% summarise(mean = mean(value)) %>% group_by(admin2) %>% summarise(val = var(mean))
plot(log(t2$val),log(t1[t2$admin2,]$est),main='Mod1_2',xlab='True log-variance',ylab='Posterior mean of log-variance')
abline(0,1)

## MOD1.5
t1 <- sig2.results.adm2 %>% filter(model == 'Mod1.5_2') %>% group_by(model,area) %>% summarise(est = mean(mean))
plot(log(params$sig2[1:300]),log(t1$est),main='Mod1.5_2')
abline(0,1)


## MOD2_1
# t1 <- sig2.results.adm1 %>% filter(model == 'Mod2_1') %>% group_by(model,area) %>% summarise(est = mean(mean))
# t2 <- sim_pop_long %>% group_by(admin1) %>% summarise(val = var(value))
# plot(log(t2$val),log(t1$est),main='Mod2_1')
# abline(0,1)


## MOD2_2
t1 <- sig2.results.adm2 %>% filter(model == 'Mod2_2') %>% group_by(model,area) %>% summarise(est = mean(mean))
t2 <- sim_pop_long %>% group_by(admin2) %>% summarise(val = var(value))
plot(log(t2$val),log(t1[t2$admin2,]$est),main='Mod2_2',xlab='True log-variance',ylab='Posterior mean of log-variance')
abline(0,1)


## MOD3_2
sig2.results.mod3_2 <- sig2.results.adm2[sig2.results.adm2$model=='Mod3_2',]
dup_id <- duplicated(sig2.results.mod3_2[,c('sim','area')])

rural.sig2 <- sig2.results.mod3_2[!dup_id,]
urban.sig2 <- sig2.results.mod3_2[dup_id,]

t1 <- rural.sig2 %>% group_by(area) %>% summarise(est = mean(mean))
t2 <- sim_pop_long %>% filter(urban==0) %>% group_by(admin2) %>% summarise(val = var(value))

plot(log(t2$val),log(t1[t2$admin2,]$est),main='Mod3_2 -- Rural',xlab='True log-variance',ylab='Posterior mean of log-variance')
abline(0,1)

t1 <- urban.sig2 %>% group_by(area) %>% summarise(est = mean(mean))
t2 <- sim_pop_long %>% filter(urban==1) %>% group_by(admin2) %>% summarise(val = var(value))
plot(log(t2$val),log(1.5*t1[t2$admin2,]$est),main='Mod3_2 -- Urban',xlab='True log-variance',ylab='Posterior mean of log-variance')
abline(0,1)

}

dev.off()

hyper.results %>% filter(param=='gamma[1]') %>% group_by(model) %>% summarise(mean(mean))
hyper.results %>% filter(param=='kappa') %>% group_by(model) %>% summarise(mean(mean))

## kappa estimation not correct..
hist(hyper.results[hyper.results$model=='Mod3_2'&hyper.results$param=='kappa',]$mean)




