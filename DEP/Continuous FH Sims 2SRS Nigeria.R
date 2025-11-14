# Simulates continuous responses for the population (using height-for-age in Nigeria as a guide)
# takes two stage stratified random samples 
# fits FH models 1) under the correct sampling model 2) with correct variance but assuming normality 3) naively (with estimated variance and assuming normality)

library(sf)
library(tidyverse)
library(INLA)
library(ggpubr)
library(locfit)
library(surveyPrev)
library(raster)
library(spdep)
library(data.table)

source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")

# load geometry ------
setwd('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm')
poly.adm2 <- st_read(dsn = 'gadm41_NGA_shp', layer = "gadm41_NGA_2", options = "ENCODING=UTF-8")
water.bodies <- poly.adm2[poly.adm2$ENGTYPE_2=='Water body',]$NAME_2
poly.adm2 <- poly.adm2[poly.adm2$ENGTYPE_2=='Local Authority',]
poly.adm2$admin2 <- 1:nrow(poly.adm2)
poly.adm2$admin2.char <- paste0('admin2_',1:nrow(poly.adm2))
n_admin2 <- nrow(poly.adm2)

poly.adm1 <- poly.adm2 %>% group_by(NAME_1) %>% summarise(geometry = st_union(geometry))
poly.adm1$admin1 <- 1:nrow(poly.adm1)
poly.adm1$admin1.char <- paste0('admin1_',1:nrow(poly.adm1))
n_admin1 <- nrow(poly.adm1)

admin.key <- merge(as.data.frame(poly.adm2[,c('NAME_1','NAME_2','admin2','admin2.char')]),
                   as.data.frame(poly.adm1[,c('NAME_1','admin1','admin1.char')]),by='NAME_1')%>%
  dplyr::select(-c(geometry.x,geometry.y))

admin1.mat <- nb2mat(poly2nb(poly.adm1), zero.policy = TRUE)
colnames(admin1.mat) <- rownames(admin1.mat) <- poly.adm1$admin1.char

admin2.mat <- nb2mat(poly2nb(poly.adm2), zero.policy = TRUE)
colnames(admin2.mat) <- rownames(admin2.mat) <- poly.adm2$admin2.char

Q.admin2 <- -admin2.mat
Q.admin2 <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})*Q.admin2
diag(Q.admin2) <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})
diag(Q.admin2)[diag(Q.admin2)==0] <- 1
Q2_inv <- INLA:::inla.ginv(as.matrix(Q.admin2))
Q2_scaled <- INLA::inla.scale.model(Q.admin2, constr=list(A=matrix(1,nrow=1,ncol=n_admin2), e=0))
Q2_scaled_inv <- INLA:::inla.ginv(as.matrix(Q2_scaled))

# load Nigeria population -------

load("~/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Nigeria/worldpop/adm2_weights_u5.rda")

# correct some labeling
weight.adm2.u5[weight.adm2.u5$admin2.name=='Ijebu East',]$admin2.char <- 'admin2_570'
weight.adm2.u5[weight.adm2.u5$admin2.name=='Ijebu North',]$admin2.char <- 'admin2_571'
weight.adm2.u5[weight.adm2.u5$admin2.name=='Ijebu North-East',]$admin2.char <- 'admin2_572'
weight.adm2.u5[weight.adm2.u5$admin2.name=='Ijebu-Ode',]$admin2.char <- 'admin2_569'

# check that all admin3 areas line up correctly
t <- merge(admin.key,weight.adm2.u5)
t[t$NAME_2!=t$admin2.name,]

# round
weight.adm2.u5$admin_pop <- round(weight.adm2.u5$admin_pop)

pop.admin2 <- merge(weight.adm2.u5[,c('admin2.char','admin_pop')],admin.key)
pop.admin2 <- pop.admin2[order(pop.admin2$admin2),]

pop.admin1 <- pop.admin2 %>% group_by(admin1) %>% summarise(admin_pop = sum(admin_pop))

A_2to1 <- matrix(0,n_admin1,n_admin2)
for(area1 in 1:n_admin1){
  which.areas2 <- unique(admin.key[admin.key$admin1==area1,]$admin2)
  A_2to1[area1,which.areas2] <- pop.admin2[pop.admin2$admin2 %in% which.areas2,]$admin_pop/sum(pop.admin2[pop.admin2$admin2 %in% which.areas2,]$admin_pop)
}

# sort population into clusters -------------------
ea_totals <- data.frame(
  NAME_1 = c(
    "Federal Capital Territory", "Benue", "Kogi", "Kwara", "Nasarawa", "Niger", "Plateau",
    "Adamawa", "Bauchi", "Borno", "Gombe", "Taraba", "Yobe",
    "Jigawa", "Kaduna", "Kano", "Katsina", "Kebbi", "Sokoto", "Zamfara",
    "Abia", "Anambra", "Ebonyi", "Enugu", "Imo",
    "Akwa Ibom", "Bayelsa", "Cross River", "Delta", "Edo", "Rivers",
    "Ekiti", "Lagos", "Ogun", "Ondo", "Osun", "Oyo"
  ),
EAs = c(
    3590, 22856, 15846, 16271, 9219, 23445, 15879,
    12808, 19885, 24086, 9494, 10600, 14923,
    21193, 21792, 36359, 33316, 16641, 12779, 17032,
    11569, 21907, 13888, 13997, 19573,
    17113, 9007, 16322, 18209, 12793, 24861,
    11561, 25424, 14493, 19255, 25907, 31106
  )
)

ea_totals <- merge(ea_totals, unique.array(admin.key[,c('admin1','admin1.char','NAME_1')]))
# how many EAs for each admin2
pop.admin2 <- pop.admin2 %>% group_by(admin1) %>% mutate(num_EAs = round(ea_totals$EAs[admin1]*admin_pop/sum(admin_pop)))

sim_pop <- as.data.table(pop.admin2)
sim_pop_long <- sim_pop[rep(seq_len(n_admin2),times = sim_pop$admin_pop)]
# assign each unit to a cluster
sim_pop_long <- sim_pop_long[, cluster := 5000*admin2 + ceiling(runif(.N) * num_EAs)]
sim_pop_long[,admin_pop := NULL]
sim_pop_long[,num_EAs := NULL]


# objects for sampling -------------
# info on each cluster that will be used to determine sampling
cluster_frame <- sim_pop_long %>% group_by(admin1,cluster) %>% summarise(N = n())
# if a cluster is selected, how many individuals will be selected?
cluster_frame$n <- rpois(nrow(cluster_frame),8) #calibrated from data

cluster_alloc <- clusters_states <- data.frame(
  NAME_1 = c(
    "Federal Capital Territory","Benue", "Kogi", "Kwara", "Nasarawa", "Niger", "Plateau",
    "Adamawa", "Bauchi", "Borno", "Gombe", "Taraba", "Yobe",
    "Jigawa", "Kaduna", "Kano", "Katsina", "Kebbi", "Sokoto", "Zamfara",
    "Abia", "Anambra", "Ebonyi", "Enugu", "Imo",
    "Akwa Ibom", "Bayelsa", "Cross River", "Delta", "Edo", "Rivers",
    "Ekiti", "Lagos", "Ogun", "Ondo", "Osun", "Oyo"
  ),
 EAs = c(
    35, 38, 36, 35, 35, 38, 35,
    35, 39, 38, 35, 35, 35,
    39, 43, 53, 42, 37, 37, 36,
    36, 39, 36, 36, 39,
    37, 35, 35, 38, 35, 41,
    35, 53, 37, 36, 36, 42
  )
)

cluster_alloc <- merge(cluster_alloc, unique.array(admin.key[,c('admin1','admin1.char','NAME_1')]))

cluster_frame <- cluster_frame %>% mutate(wt = pop.admin1$admin_pop[admin1]/(cluster_alloc$EAs[admin1]*n))

# helper functions and objects ----------
options(survey.lonely.psu = "adjust")
hyperpc.iid <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

## admin1
admin1.HT.withNA <- function(which.area) {
  admin1 <- NULL
  tmp <- subset(my.svydesign, (admin1 == as.character(which.area)))
  
  if (dim(tmp)[1] == 0) {
    return(c(which.area,rep(NA, 2)))
  } else {
    lm.ob <- survey::svymean(value ~ 1, design = tmp)
    return(c(which.area, lm.ob[1], vcov(lm.ob)))
  }
}

## admin2
admin2.HT.withNA <- function(which.area) {
  admin2 <- NULL
  tmp <- subset(my.svydesign, (admin2 == as.character(which.area)))
  
  if (dim(tmp)[1] == 0) {
    return(c(which.area,rep(NA, 2)))
  } else {
    lm.ob <- survey::svymean(value ~ 1, design = tmp)
    return(c(which.area, lm.ob[1], vcov(lm.ob)))
  }
}


# simulate data at population level -----

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
load(file = 'Nigeria_haz_admin1_weighted_estimates.rda')
load(file = 'Nigeria_haz_admin2_weighted_estimates.rda')

within_admin1_sd <- merge(admin2.dir,admin.key) %>% group_by(admin1) %>% summarise(sd = sd(mean,na.rm=T), n = sum(!is.na(mean)))
admin2_means <- sapply(1:n_admin2, function(k){
  # make area mean the weighted estimate, if it exists
  out <- admin2.dir$mean[k]
  # if it does not, draw from normal with mean equal to admin1 and variance equal to between area variation within that admin1
  if(is.na(out)){
    admin1_tmp <- admin.key[admin.key$admin2==k,]$admin1
    out <- rnorm(1,admin1.dir$mean[admin1_tmp], within_admin1_sd$sd[admin1_tmp])
  }
  return(out)
})

sd_e <- sqrt(2*var(admin2.dir$mean,na.rm=T))
#sd_e <- 0.1

sim_pop_long[, value := admin2_means[admin2] + rnorm(.N,0,sd_e)]

# what are the admin1 means
admin1_means <- sim_pop_long[order(admin1),mean(value),by=admin1]$V1

# # check the theoretical and empirical population averages are equal
# plot(admin2_means, sim_pop_long[order(admin2),mean(value),by=admin2]$V1)
# abline(0,1)
# 
# plot(A_2to1 %*% admin2_means, admin1_means)
# abline(0,1)

# draw samples and get direct estimates -------
nsim <- 1000

direct.adm2.estimates <- direct.adm1.estimates <- list()
for(k in 1:nsim){
  cat(k,'\n')
  
  # which clusters will we sample?
  cluster_sample_ids <- unlist(sapply(1:n_admin1,function(i){
    sample(cluster_frame[cluster_frame$admin1==i,]$cluster,
           cluster_alloc[cluster_alloc$admin1==i,]$EAs,
           prob = cluster_frame[cluster_frame$admin1==i,]$N)
  }))
  cluster_sample <- cluster_frame[cluster_frame$cluster %in% cluster_sample_ids,c('cluster','n','wt')]
  
  sim_sample <- sim_pop_long[
    cluster_sample,                   # right table
    on = "cluster",                   # join by cluster
    allow.cartesian = TRUE            # if some clusters repeat
  ][, .SD[sample(.N, min(.N, unique(n)))], by = cluster] # for each cluster, sample n
  
  sim_sample[,n:=NULL]
  
  # get direct estimates --------
  options(survey.lonely.psu = "adjust")
  sim_sample$admin1.char <- factor(sim_sample$admin1.char)
  my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"), 
                                    strata = ~admin1.char, nest = T, weights = ~wt, data = sim_sample)

  ## admin2
  x <- mapply(admin2.HT.withNA, which.area = 1:n_admin2)
  admin2.dir <- data.frame(t(x))
  colnames(admin2.dir) <- c('admin2','mean','variance')
  
  direct.adm2.estimates[[k]] <- admin2.dir
  
  ## admin1
  x <- mapply(admin1.HT.withNA, which.area = 1:n_admin1)
  admin1.dir <- data.frame(t(x))
  colnames(admin1.dir) <- c('admin1','mean','variance')
  
  direct.adm1.estimates[[k]] <- admin1.dir
  
}

# get FH estimates from direct --------
# true variance of design-based estimates
adm2.dir.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct.adm2.estimates, function(x){x$mean[i]})),na.rm = T)})
adm1.dir.variance <- sapply(1:n_admin1, function(i){var(unlist(lapply(direct.adm1.estimates, function(x){x$mean[i]})),na.rm = T)})

results.adm1 <- results.adm2 <- list() # record estimates (direct and fay herriot)
for(k in 1:100){
  cat(k,'\n')
  
  sim.admin2.dir <- direct.adm2.estimates[[k]]
  sim.admin2.dir$var.truth <- adm2.dir.variance
  sim.admin2.dir$mean.sim <- sapply(1:n_admin2,function(j){
    rnorm(1,admin2_means[j],sqrt(adm2.dir.variance[j]))
  })
  sim.admin2.dir <- sim.admin2.dir %>% mutate(variance = ifelse(variance < 1e-5,NA,variance))
  sim.admin2.dir <- sim.admin2.dir %>% mutate(mean = ifelse(is.na(variance),NA,mean),
                                              var.truth = ifelse(is.na(variance),NA,var.truth),
                                              mean.sim = ifelse(is.na(variance),NA,mean.sim))
  
  
  sim.admin1.dir <- direct.adm1.estimates[[k]]
  sim.admin1.dir$var.truth <- adm1.dir.variance
  sim.admin1.dir$mean.sim <- sapply(1:n_admin1,function(j){
    rnorm(1,admin1_means[j],sqrt(adm1.dir.variance[j]))
  })
  
  # # model is specified exactly correctly ----
  fit1a <- inla(mean.sim ~ f(admin1,model='iid',hyper = hyperpc.iid),
                family = 'gaussian',
                data = sim.admin1.dir,
                scale = 1/sim.admin1.dir$var.truth,
                control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                control.compute = list(config = TRUE))
  samples1a <- get_inla_samples_local(fit1a,1000)
  
  eta.samples1a <- samples1a[,n_admin1+1] + samples1a[,1:n_admin1]
  
  fit1 <- inla(mean.sim ~ f(admin2,model='iid',hyper = hyperpc.iid),
               family = 'gaussian',
               data = sim.admin2.dir,
               scale = 1/sim.admin2.dir$var.truth,
               control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
               control.compute = list(config = TRUE))
  samples1 <- get_inla_samples_local(fit1,1000)
  
  eta.samples1 <- samples1[,n_admin2+1] + samples1[,1:n_admin2]
  
  ## model is specified with correct variance, but assuming normality ----
  fit2a <- inla(mean ~ f(admin1,model='iid',hyper = hyperpc.iid),
                family = 'gaussian',
                data = sim.admin1.dir,
                scale = 1/sim.admin1.dir$var.truth,
                control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                control.compute = list(config = TRUE))
  samples2a <- get_inla_samples_local(fit2a,1000)
  
  eta.samples2a <- samples2a[,n_admin1+1] + samples2a[,1:n_admin1]
  
  fit2 <- inla(mean ~ f(admin2,model='iid',hyper = hyperpc.iid),
               family = 'gaussian',
               data = sim.admin2.dir,
               scale = 1/sim.admin2.dir$var.truth,
               control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
               control.compute = list(config = TRUE))
  samples2 <- get_inla_samples_local(fit2,1000)
  
  eta.samples2 <- samples2[,n_admin2+1] + samples2[,1:n_admin2]
  
  ## real way of doing FH (estimated variance, assuming normality) ----
  fit3a <- inla(mean ~ f(admin1,model='iid',hyper = hyperpc.iid),
                family = 'gaussian',
                data = sim.admin1.dir,
                scale = 1/sim.admin1.dir$variance,
                control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                control.compute = list(config = TRUE))
  samples3a <- get_inla_samples_local(fit3a,1000)
  
  eta.samples3a <- samples3a[,n_admin1+1] + samples3a[,1:n_admin1]
  
  fit3 <- inla(mean ~ f(admin2,model='iid',hyper = hyperpc.iid),
               family = 'gaussian',
               data = sim.admin2.dir,
               scale = 1/sim.admin2.dir$variance,
               control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
               control.compute = list(config = TRUE))
  samples3 <- get_inla_samples_local(fit3,1000)
  
  eta.samples3 <- samples3[,n_admin2+1] + samples3[,1:n_admin2]
  
  ## record estimates -----
  
  results.adm2[[k]] <- merge(sim.admin2.dir[,c('admin2','mean','mean.sim','variance','var.truth')],
                             data.frame(admin2 = 1:n_admin2,
                                        # correctly specified
                                        mean.fh1  = apply(eta.samples1,2,mean),
                                        var.fh1 = apply(eta.samples1,2,var),
                                        lower.fh1 = apply(eta.samples1,2,quantile,prob=0.05),
                                        upper.fh1 = apply(eta.samples1,2,quantile,prob=0.95),
                                        # uses correct variance, but normality is assumed
                                        mean.fh2  = apply(eta.samples2,2,mean),
                                        var.fh2 = apply(eta.samples2,2,var),
                                        lower.fh2 = apply(eta.samples2,2,quantile,prob=0.05),
                                        upper.fh2 = apply(eta.samples2,2,quantile,prob=0.95),
                                        # assumes and normality and uses estimated variance
                                        mean.fh3  = apply(eta.samples3,2,mean),
                                        var.fh3 = apply(eta.samples3,2,var),
                                        lower.fh3 = apply(eta.samples3,2,quantile,prob=0.05),
                                        upper.fh3 = apply(eta.samples3,2,quantile,prob=0.95),
                                        sim = k,
                                        true_mean = admin2_means),all=T)
  
  results.adm1[[k]] <- merge(sim.admin1.dir[,c('admin1','mean','mean.sim','variance','var.truth')],
                             data.frame(admin1 = 1:n_admin1,
                                        # correctly specified
                                        mean.fh1  = apply(eta.samples1a,2,mean),
                                        var.fh1 = apply(eta.samples1a,2,var),
                                        lower.fh1 = apply(eta.samples1a,2,quantile,prob=0.05),
                                        upper.fh1 = apply(eta.samples1a,2,quantile,prob=0.95),
                                        # uses correct variance, but normality is assumed
                                        mean.fh2  = apply(eta.samples2a,2,mean),
                                        var.fh2 = apply(eta.samples2a,2,var),
                                        lower.fh2 = apply(eta.samples2a,2,quantile,prob=0.05),
                                        upper.fh2 = apply(eta.samples2a,2,quantile,prob=0.95),
                                        # assumes and normality and uses estimated variance
                                        mean.fh3  = apply(eta.samples3a,2,mean),
                                        var.fh3 = apply(eta.samples3a,2,var),
                                        lower.fh3 = apply(eta.samples3a,2,quantile,prob=0.05),
                                        upper.fh3 = apply(eta.samples3a,2,quantile,prob=0.95),
                                        sim = k,
                                        true_mean = admin1_means),all=T)
  
}

## SAVE/LOAD RESULTS ###################
setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
# 
# save(sim_pop_long,file='sim4_popdata.rda')
# save(cluster_frame,file='sim4_clusterframe.rda')
# save(cluster_alloc,file='sim4_clusteralloc.rda')
# save(direct.adm2.estimates,file='sim4_admin2.direct.rda')
# save(direct.adm1.estimates,file='sim4_admin1.direct.rda')
# save(results.adm2,file = 'sim4_admin2.fh.rda')
# save(results.adm1,file = 'sim4_admin1.fh.rda')

load(file='sim3_popdata.rda')
load(file='sim3_clusterframe.rda')
load(file='sim3_clusteralloc.rda')
load(file='sim3_admin2.direct.rda')
load(file='sim3_admin1.direct.rda')
load(file = 'sim3_admin2.fh.rda')
load(file = 'sim3_admin1.fh.rda')

#### POST PROC -----

res.adm1.mat <- do.call(rbind,results.adm1)
res.adm1.mat <- res.adm1.mat %>% mutate(lower.direct = mean - 1.645*sqrt(variance),
                                        upper.direct = mean + 1.645*sqrt(variance),
                                        ciwidth0 = upper.direct - lower.direct,
                                        ciwidth1 = upper.fh1 - lower.fh1,
                                        ciwidth2 = upper.fh2 - lower.fh2,
                                        ciwidth3 = upper.fh3 - lower.fh3,
                                        cov0 = lower.direct<true_mean & upper.direct>true_mean,
                                        cov1 = lower.fh1<true_mean & upper.fh1>true_mean,
                                        cov2 = lower.fh2<true_mean & upper.fh2>true_mean,
                                        cov3 = lower.fh3<true_mean & upper.fh3>true_mean,
                                        i.score0 = ciwidth0 + 20*(lower.direct-true_mean)*(true_mean<lower.direct) + 20*(true_mean-upper.direct)*(true_mean>upper.direct),
                                        i.score1 = ciwidth1 + 20*(lower.fh1-true_mean)*(true_mean<lower.fh1) + 20*(true_mean-upper.fh1)*(true_mean>upper.fh1),
                                        i.score2 = ciwidth2 + 20*(lower.fh2-true_mean)*(true_mean<lower.fh2) + 20*(true_mean-upper.fh2)*(true_mean>upper.fh2),
                                        i.score3 = ciwidth3 + 20*(lower.fh3-true_mean)*(true_mean<lower.fh3) + 20*(true_mean-upper.fh3)*(true_mean>upper.fh3))

res.adm2.mat <- do.call(rbind,results.adm2)
res.adm2.mat <- res.adm2.mat %>% mutate(lower.direct = mean - 1.645*sqrt(variance),
                                        upper.direct = mean + 1.645*sqrt(variance),
                                        ciwidth0 = upper.direct - lower.direct,
                                        ciwidth1 = upper.fh1 - lower.fh1,
                                        ciwidth2 = upper.fh2 - lower.fh2,
                                        ciwidth3 = upper.fh3 - lower.fh3,
                                        cov0 = lower.direct<true_mean & upper.direct>true_mean,
                                        cov1 = lower.fh1<true_mean & upper.fh1>true_mean,
                                        cov2 = lower.fh2<true_mean & upper.fh2>true_mean,
                                        cov3 = lower.fh3<true_mean & upper.fh3>true_mean,
                                        i.score0 = ciwidth0 + 20*(lower.direct-true_mean)*(true_mean<lower.direct) + 20*(true_mean-upper.direct)*(true_mean>upper.direct),
                                        i.score1 = ciwidth1 + 20*(lower.fh1-true_mean)*(true_mean<lower.fh1) + 20*(true_mean-upper.fh1)*(true_mean>upper.fh1),
                                        i.score2 = ciwidth2 + 20*(lower.fh2-true_mean)*(true_mean<lower.fh2) + 20*(true_mean-upper.fh2)*(true_mean>upper.fh2),
                                        i.score3 = ciwidth3 + 20*(lower.fh3-true_mean)*(true_mean<lower.fh3) + 20*(true_mean-upper.fh3)*(true_mean>upper.fh3))

coverage.adm1 <- res.adm1.mat %>% group_by(admin1,true_mean) %>% reframe(area_coverage0 = mean(cov0,na.rm=T),
                                                                         area_coverage1 = mean(cov1),
                                                                         area_coverage2 = mean(cov2),
                                                                         area_coverage3 = mean(cov3))


coverage.adm2 <- res.adm2.mat %>% group_by(admin2,true_mean) %>% reframe(area_coverage0 = mean(cov0,na.rm=T),
                                                                         area_coverage1 = mean(cov1),
                                                                         area_coverage2 = mean(cov2),
                                                                         area_coverage3 = mean(cov3))

coverage.adm2.insample <- res.adm2.mat %>% filter(!is.na(mean)) %>% 
  group_by(admin2,true_mean) %>% reframe(area_coverage0 = mean(cov0),
                                         area_coverage1 = mean(cov1),
                                         area_coverage2 = mean(cov2),
                                         area_coverage3 = mean(cov3))

# plots ################

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
pdf('sim3_admin1_plots.pdf')
{
  
  ## ADMIN1
  # coverage
  y_range <- range(coverage.adm1$area_coverage0,coverage.adm1$area_coverage1,coverage.adm1$area_coverage2,coverage.adm1$area_coverage3)
  p0 <- coverage.adm1 %>% ggplot() + geom_point(aes(true_mean,area_coverage0)) + 
    xlab('') + ylab('') + ggtitle('Weighted estimates') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm1$area_coverage0),3)),
               hjust=-0.1, vjust = 1.1)
  p1 <- coverage.adm1 %>% ggplot() + geom_point(aes(true_mean,area_coverage1)) + 
    xlab('') + ylab('') + ggtitle('Correct sampling model') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm1$area_coverage1),3)),
               hjust=-0.1, vjust = 1.1)
  p2 <- coverage.adm1 %>% ggplot() + geom_point(aes(true_mean,area_coverage2)) + 
    xlab('') + ylab('') +ggtitle('Correct variance') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm1$area_coverage2),3)),
               hjust=-0.1, vjust = 1.1)
  p3 <- coverage.adm1 %>% ggplot() + geom_point(aes(true_mean,area_coverage3)) + 
    xlab('') + ylab('') +ggtitle('Naive') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm1$area_coverage3),3)),
               hjust=-0.1, vjust = 1.1)
  
  p <- ggarrange(plotlist = list(p0,p1,p2,p3))
  print(annotate_figure(p, top = text_grob("Coverage of 90% intervals for admin1 estimates", face = "bold", size = 14)))
  
  ## interval widths
  
  tmp <- res.adm1.mat %>% group_by(admin1) %>% reframe(weighted = mean(ciwidth0),
                                                                                correct_sampling = mean(ciwidth1),
                                                                                correct_variance = mean(ciwidth2),
                                                                                naive = mean(ciwidth3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) 
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,correct_sampling)) + geom_abline(intercept = 0, slope=1)
  p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)
  
  p <- ggarrange(plotlist = list(p1,p2,p3),nrow=3)
  print(annotate_figure(p, top = text_grob("Average width of 90% intervals for admin1 estimates", face = "bold", size = 14)))
  
  ## ais
  
  tmp <- res.adm1.mat %>% group_by(admin1) %>% reframe(weighted = mean(i.score0),
                                                       correct_sampling = mean(i.score1),
                                                       correct_variance = mean(i.score2),
                                                       naive = mean(i.score3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) 
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,correct_sampling)) + geom_abline(intercept = 0, slope=1)
  p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)
  
  p <- ggarrange(plotlist = list(p1,p2,p3),nrow=3)
  print(annotate_figure(p, top = text_grob("AIS for 90% intervals for admin1 estimates", face = "bold", size = 14)))
  
  # estimated v true variance
  g3 <- res.adm1.mat %>% group_by(admin1,var.truth) %>% reframe(var_direct = mean(variance)) %>%
    ggplot() + 
    geom_point(aes(var_direct,var.truth)) + 
    geom_abline(slope=1,intercept = 0)  + 
    xlab('Average of estimated variance') + ylab('True variance') +
    ggtitle('Design-based variance of the admin1 weighted estimates')
  
  print(g3)
}
dev.off()

## ADMIN2

pdf('sim3_admin2_plots.pdf')
{
  # coverage
  y_range <- range(coverage.adm2$area_coverage0,coverage.adm2$area_coverage1,coverage.adm2$area_coverage2,coverage.adm2$area_coverage3,na.rm=T)
  p0 <- coverage.adm2 %>% ggplot() + geom_point(aes(true_mean,area_coverage0)) + 
    xlab('') + ylab('') + ggtitle('Weighted estimates') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2$area_coverage0,na.rm=T),3)),
               hjust=-0.1, vjust = 1.1)
  p1 <- coverage.adm2.insample %>% ggplot() + geom_point(aes(true_mean,area_coverage1)) + 
    xlab('') + ylab('') + ggtitle('Correct sampling model') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2$area_coverage1),3)),
               hjust=-0.1, vjust = 1.1)
  p2 <- coverage.adm2.insample %>% ggplot() + geom_point(aes(true_mean,area_coverage2)) + 
    xlab('') + ylab('') +ggtitle('Correct variance') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2$area_coverage2),3)),
               hjust=-0.1, vjust = 1.1)
  p3 <- coverage.adm2.insample %>% ggplot() + geom_point(aes(true_mean,area_coverage3)) + 
    xlab('') + ylab('') +ggtitle('Naive') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2$area_coverage3),3)),
               hjust=-0.1, vjust = 1.1)
  
  p <- ggarrange(plotlist = list(p0,p1,p2,p3))
  print(annotate_figure(p, top = text_grob("Coverage of 90% intervals for in-sample admin2 estimates", face = "bold", size = 14)))
  
  p1 <- coverage.adm2 %>% ggplot() + geom_point(aes(true_mean,area_coverage1)) + 
    xlab('') + ylab('') + ggtitle('Correct sampling model') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2$area_coverage1),3)),
               hjust=-0.1, vjust = 1.1)
  p2 <- coverage.adm2 %>% ggplot() + geom_point(aes(true_mean,area_coverage2)) + 
    xlab('') + ylab('') +ggtitle('Correct variance') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2$area_coverage2),3)),
               hjust=-0.1, vjust = 1.1)
  p3 <- coverage.adm2 %>% ggplot() + geom_point(aes(true_mean,area_coverage3)) + 
    xlab('') + ylab('') +ggtitle('Naive') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2$area_coverage3),3)),
               hjust=-0.1, vjust = 1.1)
  
  p <- ggarrange(plotlist = list(p1,p2,p3),nrow=3)
  print(annotate_figure(p, top = text_grob("Coverage of 90% intervals for all admin2 estimates", face = "bold", size = 14)))
  
  
  ## interval widths
  
  # g1 <- res.adm2.mat %>% 
  #   dplyr::select(admin2,ciwidth0:ciwidth3) %>% 
  #   pivot_longer(ciwidth0:ciwidth3,names_to = 'model',values_to = 'width') %>%
  #   ggplot() + geom_boxplot(aes(model,width)) + ylab('') + xlab('') + 
  #   ggtitle('Width of 90% intervals for admin2 estimates') + theme_bw()
  # 
  # print(g1)
  
  tmp <- res.adm2.mat %>% filter(!is.na(mean)) %>% group_by(admin2) %>% reframe(weighted = mean(ciwidth0),
                                                       correct_sampling = mean(ciwidth1),
                                                       correct_variance = mean(ciwidth2),
                                                       naive = mean(ciwidth3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) +
    ggtitle('in-sample estimates')
 # p2 <- tmp %>% ggplot() + geom_point(aes(naive,correct_sampling)) + geom_abline(intercept = 0, slope=1)
 # p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)
  
  # p <- ggarrange(plotlist = list(p1,p2,p3),nrow=3)
  # print(annotate_figure(p, top = text_grob("Average width of 90% intervals for in-sample admin2 estimates", face = "bold", size = 14)))
  # 
  tmp <- res.adm2.mat %>% group_by(admin2) %>% reframe(correct_sampling = mean(ciwidth1),
                                                       correct_variance = mean(ciwidth2),
                                                       naive = mean(ciwidth3))
  
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,correct_sampling)) + geom_abline(intercept = 0, slope=1)+
    ggtitle('all estimates')
  p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)+
    ggtitle('all estimates')
  
  p <- ggarrange(plotlist = list(p1,p2,p3),nrow=3)
  print(annotate_figure(p, top = text_grob("Average width of 90% intervals for admin2 estimates", face = "bold", size = 14)))
  
  # AIS
  tmp <- res.adm2.mat %>% filter(!is.na(mean)) %>%
    group_by(admin2,true_mean) %>% summarise(weighted = mean(i.score0),
                                             correct_sampling = mean(i.score1),
                                             correct_variance = mean(i.score2),
                                             naive = mean(i.score3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) +
    ggtitle('in-sample estimates')
  # p2 <- tmp %>% ggplot() + geom_point(aes(naive,correct_sampling)) + geom_abline(intercept = 0, slope=1)
  # p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)
  # 
  # p <- ggarrange(plotlist = list(p1,p2,p3),nrow=3)
  # print(annotate_figure(p, top = text_grob("AIS of 90% intervals for in-sample admin2 estimates", face = "bold", size = 14)))
  # 
  tmp <- res.adm2.mat %>% 
    group_by(admin2,true_mean) %>% summarise(correct_sampling = mean(i.score1),
                                             correct_variance = mean(i.score2),
                                             naive = mean(i.score3))
  
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,correct_sampling)) + geom_abline(intercept = 0, slope=1) +
    ggtitle('all estimates')
  p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)+
    ggtitle('all estimates')
  
  p <- ggarrange(plotlist = list(p1,p2,p3),nrow=3)
  print(annotate_figure(p, top = text_grob("AIS of 90% intervals for admin2 estimates", face = "bold", size = 14)))
  
  
  # estimated v true variance
  
  g3 <- res.adm2.mat %>% group_by(admin2,var.truth) %>% reframe(var_direct = mean(variance)) %>%
    ggplot() + 
    geom_point(aes(var_direct,var.truth)) + 
    geom_abline(slope=1,intercept = 0) + 
    xlab('Average of estimated variance') + ylab('True variance') +
    ggtitle('Design-based variance of the admin2 weighted estimates')
  
  print(g3)
  
}
dev.off()




