# Simulates continuous responses for the population (using weight-for-height in Rwanda as a guide)
# takes two stage stratified random samples (at the admin2 level and cluster level)
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
# using Rwanda data as skeleton

# load geometry ------
setwd('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm')
poly.adm3 <- st_read(dsn = 'gadm41_RWA_shp', layer = "gadm41_RWA_3", options = "ENCODING=UTF-8")
water.bodies <- poly.adm3[poly.adm3$ENGTYPE_3=='Water body',]$NAME_3
poly.adm3 <- poly.adm3[poly.adm3$ENGTYPE_3=='Sector',]
poly.adm3$admin3 <- 1:nrow(poly.adm3)
poly.adm3$admin3.char <- paste0('admin3_',1:nrow(poly.adm3))
n_admin3 <- nrow(poly.adm3)

poly.adm2 <- poly.adm3 %>% group_by(NAME_2) %>% summarise(geometry = st_union(geometry))
poly.adm2$admin2 <- 1:nrow(poly.adm2)
poly.adm2$admin2.char <- paste0('admin2_',1:nrow(poly.adm2))
n_admin2 <- nrow(poly.adm2)

poly.adm1 <- poly.adm3 %>% group_by(NAME_1) %>% summarise(geometry = st_union(geometry))
poly.adm1$admin1 <- 1:nrow(poly.adm1)
poly.adm1$admin1.char <- paste0('admin1_',1:nrow(poly.adm1))
n_admin1 <- nrow(poly.adm1)

admin.key <- merge(as.data.frame(poly.adm2[,c('NAME_2','admin2','admin2.char')]),
                   as.data.frame(poly.adm3[,c('NAME_1','NAME_2','NAME_3','admin3','admin3.char')]),by='NAME_2')

admin.key <- merge(admin.key, as.data.frame(poly.adm1[,c('NAME_1','admin1','admin1.char')])) %>%
  dplyr::select(-c(geometry.x,geometry.y,geometry))
admin.key <-admin.key[order(admin.key$admin3),]

admin1.mat <- nb2mat(poly2nb(poly.adm1), zero.policy = TRUE)
colnames(admin1.mat) <- rownames(admin1.mat) <- poly.adm1$admin1.char

admin2.mat <- nb2mat(poly2nb(poly.adm2), zero.policy = TRUE)
colnames(admin2.mat) <- rownames(admin2.mat) <- poly.adm2$admin2.char

admin3.mat <- nb2mat(poly2nb(poly.adm3), zero.policy = TRUE)
colnames(admin3.mat) <- rownames(admin3.mat) <- poly.adm3$admin3.char

Q.admin3 <- -admin3.mat
Q.admin3 <- sapply(1:nrow(Q.admin3),function(i){sum(I(Q.admin3[i,]!=0))})*Q.admin3
diag(Q.admin3) <- sapply(1:nrow(Q.admin3),function(i){sum(I(Q.admin3[i,]!=0))})
diag(Q.admin3)[diag(Q.admin3)==0] <- 1
Q3_inv <- INLA:::inla.ginv(as.matrix(Q.admin3))
Q3_scaled <- INLA::inla.scale.model(Q.admin3, constr=list(A=matrix(1,nrow=1,ncol=n_admin3), e=0))
Q3_scaled_inv <- INLA:::inla.ginv(as.matrix(Q3_scaled))

# load Rwanda population and modify -------

load(file = "/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/RWA_Admin3_pop_u5.rda")## add some population to areas without any (its fake data anyway)
which(pop.admin3$population==0)
pop.admin3[112,]$population <- 1000
# round
pop.admin3$population <- round(pop.admin3$population)

pop.admin2 <- pop.admin3 %>% group_by(admin2) %>% summarise(population = sum(population))

#  pre <- "https://data.worldpop.org/GIS/AgeSex_structures/"
# f0 <- paste0(pre, "Global_2000_2020/2019/RWA/rwa_f_0_2019.tif")
# m0 <- paste0(pre, "Global_2000_2020/2019/RWA/rwa_m_0_2019.tif")
# f1 <- paste0(pre, "Global_2000_2020/2019/RWA/rwa_f_1_2019.tif")
# m1 <- paste0(pre, "Global_2000_2020/2019/RWA/rwa_m_1_2019.tif")
# 
# pop_f_0 <- raster(f0)
# pop_m_0 <- raster(m0)
# pop_f_1 <- raster(f1)
# pop_m_1 <- raster(m1)
# 
# pop_raster <- pop_f_0 + pop_m_0 + pop_f_1 + pop_m_1
# 
# pop.admin3 <- aggPopulation(tiff = pop_raster,
#                           poly.adm = poly.adm3,
#                           by.adm = "admin3.char")
# 
# pop.admin3 <- pop.admin3 %>% rename(admin3.char = admin1.name)
# pop.admin3 <- rbind(pop.admin3,
#                     data.frame(admin3.char = admin.key$admin3.char[!(admin.key$admin3.char %in%  pop.admin3$admin3.char)],
#                                population = 0))
# 
# pop.admin3 <- merge(pop.admin3,admin.key)
# pop.admin3 <- pop.admin3[order(pop.admin3$admin3),]
# 
# save(pop.admin3,file = "/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/RWA_Admin3_pop_u5.rda")

A_3to2 <- matrix(0,n_admin2,n_admin3)
for(area2 in 1:n_admin2){
  which.areas3 <- unique(admin.key[admin.key$admin2==area2,]$admin3)
  A_3to2[area2,which.areas3] <- pop.admin3[pop.admin3$admin3 %in% which.areas3,]$population/sum(pop.admin3[pop.admin3$admin3 %in% which.areas3,]$population)
}

A_2to1 <- matrix(0,n_admin1,n_admin2)
for(area1 in 1:n_admin1){
  which.areas2 <- unique(admin.key[admin.key$admin1==area1,]$admin2)
  A_2to1[area1,which.areas2] <- pop.admin2[pop.admin2$admin2 %in% which.areas2,]$population/sum(pop.admin2[pop.admin2$admin2 %in% which.areas2,]$population)
}




# sort population into clusters --------
ea_dist <- data.frame(
  District = c(
    "Nyarugenge","Gasabo","Kicukiro","Total",
    "Nyanza","Gisagara","Nyaruguru","Huye","Nyamagabe","Ruhango","Muhanga","Kamonyi","Total",
    "Karongi","Rutsiro","Rubavu","Nyabihu","Ngororero","Rusizi","Nyamasheke","Total",
    "Rulindo","Gakenke","Musanze","Burera","Gicumbi","Total",
    "Rwamagana","Nyagatare","Gatsibo","Kayonza","Kirehe","Ngoma","Bugesera","Total",
    "Rwanda"
  ),
  # Urban_EAs = c(
  #   396,585,473,1454,
  #   36,9,8,64,31,40,49,41,278,
  #   35,9,203,44,16,83,8,398,
  #   11,17,116,10,34,188,
  #   39,59,28,35,17,20,38,236,
  #   2554
  # ),
  # Rural_EAs = c(
  #   122,262,72,456,
  #   432,533,391,486,525,511,361,386,3625,
  #   511,482,375,445,484,543,602,3442,
  #   492,603,405,582,611,2693,
  #   467,635,643,426,613,510,576,3870,
  #   14086
  # ),
  Total_EAs = c(
    518,847,545,1910,
    468,542,399,550,556,551,410,427,3903,
    546,491,578,489,500,626,610,3840,
    503,620,521,592,645,2881,
    506,694,671,461,630,530,614,4106,
    16640
  )#,
  # AvgEASize_Urban = c(
  #   135,171,145,153,
  #   181,138,174,177,159,163,213,235,187,
  #   169,162,169,197,189,160,174,171,
  #   190,147,201,150,166,186,
  #   170,206,210,212,139,168,190,191,
  #   165
  # ),
  # AvgEASize_Rural = c(
  #   142,159,139,151,
  #   159,138,153,138,134,137,175,185,151,
  #   133,145,146,129,157,130,134,139,
  #   133,128,152,124,132,133,
  #   145,149,140,166,123,150,136,143,
  #   142
  # ),
  # AvgEASize_Total = c(
  #   137,168,144,153,
  #   160,143,154,142,135,139,180,190,153,
  #   135,145,154,135,158,134,135,142,
  #   134,129,163,124,134,136,
  #   147,154,143,170,123,151,139,146,
  #   146
  # )
)

ea_dist <- ea_dist %>% filter(!(District %in% c('Total','Rwanda')))
ea_dist <- merge(ea_dist,unique.array(admin.key[,c('admin2','admin2.char','NAME_2')]),by.x = 'District',by.y = 'NAME_2')

# how many EAs for each admin3
pop.admin3 <- pop.admin3 %>% group_by(admin2) %>% mutate(num_EAs = round(ea_dist$Total_EAs[admin2]*population/sum(population)))

sim_pop <- as.data.table(pop.admin3)
sim_pop_long <- sim_pop[rep(seq_len(n_admin3),times = sim_pop$population)]
# assign each unit to a cluster
sim_pop_long <- sim_pop_long[, cluster := 1000*admin3 + round(runif(.N,0,num_EAs))]
sim_pop_long[,population := NULL]
sim_pop_long[,num_EAs := NULL]

## objects for sampling --------
# info on each cluster that will be used to determine sampling
cluster_frame <- sim_pop_long %>% group_by(admin2,cluster) %>% summarise(N = n())
# if a cluster is selected, how many individuals will be selected?
cluster_frame$n <- pmax(0,round(rnorm(nrow(cluster_frame),7.5,3)))

cluster_alloc <- data.frame(
  District = c(
    "Nyarugenge","Gasabo","Kicukiro","Total",
    "Nyanza","Gisagara","Nyaruguru","Huye","Nyamagabe","Ruhango","Muhanga","Kamonyi","Total",
    "Karongi","Rutsiro","Rubavu","Nyabihu","Ngororero","Rusizi","Nyamasheke","Total",
    "Rulindo","Gakenke","Musanze","Burera","Gicumbi","Total",
    "Rwamagana","Nyagatare","Gatsibo","Kayonza","Kirehe","Ngoma","Bugesera","Total"
  ),
  EAs_Total = c(
    20,21,20,61,
    18,16,16,16,16,17,16,16,131,
    16,16,16,16,16,16,16,112,
    16,16,16,16,16,80,
    16,17,16,16,18,16,17,116))

cluster_alloc <- cluster_alloc %>% filter(!(District %in% c('Total','Rwanda')))
cluster_alloc  <- merge(cluster_alloc ,unique.array(admin.key[,c('admin2','admin2.char','NAME_2')]),by.x = 'District',by.y = 'NAME_2')

# double sample size
# cluster_alloc$EAs_Total <- 2*cluster_alloc$EAs_Total

cluster_frame <- cluster_frame %>% mutate(wt = pop.admin2$population[admin2]/(cluster_alloc$EAs_Total[admin2]*n))

# helper functions and objects ----------
hyperpc.bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                     phi = list(prior = "pc", param = c(0.5, 2/3)))

## admin2
admin2.HT.withNA <- function(which.area) {
  admin2 <- NULL
  tmp <- subset(my.svydesign, (admin2 == as.character(which.area)))
  
  if (dim(tmp)[1] == 0) {
    return(rep(NA, 2))
  } else {
    lm.ob <- survey::svymean(y ~ 1, design = tmp)
    return(c(which.area, lm.ob[1], vcov(lm.ob)))
  }
}

## admin3
admin3.HT.withNA <- function(which.area) {
  admin3 <- NULL
  tmp <- subset(my.svydesign, (admin3 == as.character(which.area)))
  
  if (dim(tmp)[1] == 0) {
    return(c(which.area,rep(NA, 2)))
  } else {
    lm.ob <- survey::svymean(y ~ 1, design = tmp)
    return(c(which.area, lm.ob[1], vcov(lm.ob)))
  }
}


# simulate data at population level ------
alpha <- -1
sd_b <- 0.5
phi <- 0.75
sd_e <- c(0.5,1)[2]

# iid
#b <- rnorm(n_admin3,0,sd_b)
#bym2
u <- as.vector(Rfast::rmvnorm(1,rep(0,n_admin3),sigma = Q3_scaled_inv))
v <- rnorm(n_admin3,0,1)
b <- sd_b*(sqrt(1-phi)*v + sqrt(phi)*u)
mu <- alpha + b

sim_pop_long[, y := mu[admin3] + rnorm(.N,  0, sd_e)]

# # what is the mean at the admin3 level?
admin3.means <- sim_pop_long[order(admin3),mean(y),by=admin3]
# plot(mu, admin3.means$V1)
# abline(0,1)

# # what is the mean at the admin2 level?
admin2.means <- sim_pop_long[order(admin2),mean(y),by=admin2]
# plot(A_3to2 %*% mu, admin2.means$V1)
# abline(0,1)

# draw samples and get direct estimates -------
nsim <- 1000

direct.adm2.estimates <- direct.adm3.estimates <- list()
for(k in 1:nsim){
  cat(k,'\n')
  
# which clusters will we sample?
cluster_sample_ids <- unlist(sapply(1:n_admin2,function(i){
  sample(cluster_frame[cluster_frame$admin2==i,]$cluster,
         cluster_alloc[cluster_alloc$admin2==i,]$EAs_Total,
         prob = cluster_frame[cluster_frame$admin2==i,]$N)
}))
cluster_sample <- cluster_frame[cluster_frame$cluster %in% cluster_sample_ids,c('cluster','n','wt')]

sim_sample_list <- apply(cluster_sample,1,function(row){
  sim_pop_long[cluster == as.numeric(row['cluster']),][sample(.N,as.numeric(row['n']))]
})

sim_sample <- do.call(rbind,sim_sample_list)
sim_sample <- merge(sim_sample,cluster_sample[,c('cluster','wt')])

# get direct estimates --------
options(survey.lonely.psu = "adjust")
sim_sample$admin2.char <- factor(sim_sample$admin2.char)
my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"), 
                                  strata = ~admin2.char, nest = T, weights = ~wt, data = sim_sample)


## admin2
x <- mapply(admin2.HT.withNA, which.area = 1:n_admin2)
admin2.dir <- data.frame(t(x))
colnames(admin2.dir) <- c('admin2','mean','variance')

direct.adm2.estimates[[k]] <- admin2.dir

## admin3
x <- mapply(admin3.HT.withNA, which.area = 1:n_admin3)
admin3.dir <- data.frame(t(x))
colnames(admin3.dir) <- c('admin3','mean','variance')

direct.adm3.estimates[[k]] <- admin3.dir

}

# get FH estimates from direct --------
# true variance of design-based estimates
adm2.dir.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct.adm2.estimates, function(x){x$mean[i]})),na.rm = T)})
adm3.dir.variance <- sapply(1:n_admin3, function(i){var(unlist(lapply(direct.adm3.estimates, function(x){x$mean[i]})),na.rm = T)})

results.adm3 <- results.adm2 <- list() # record estimates (direct and fay herriot)
for(k in 1:100){
  cat(k,'\n')
  
  sim.admin3.dir.stable <- direct.adm3.estimates[[k]]
  sim.admin3.dir.stable <- sim.admin3.dir.stable[!is.na(sim.admin3.dir.stable$mean),]
  sim.admin3.dir.stable <- sim.admin3.dir.stable[sim.admin3.dir.stable$variance>1e-5,]
  sim.admin3.dir.stable$var.truth <- adm3.dir.variance[sim.admin3.dir.stable$admin3]
  sim.admin3.dir.stable$mean.sim <- rnorm(nrow(sim.admin3.dir.stable),
                                          admin3.means$V1[sim.admin3.dir.stable$admin3],
                                          sqrt(sim.admin3.dir.stable$var.truth))
  
  sim.admin2.dir <- direct.adm2.estimates[[k]]
  sim.admin2.dir$var.truth <- adm2.dir.variance[sim.admin2.dir$admin2]
  sim.admin2.dir$mean.sim <- rnorm(nrow(sim.admin2.dir),
                                   admin2.means$V1[sim.admin2.dir$admin2],
                                   sqrt(sim.admin2.dir$var.truth))
  
  # model is specified exactly correctly ----
  fit1a <- inla(mean.sim ~ f(admin2,model='bym2',hyper = hyperpc.bym2, graph = admin2.mat,
                             scale.model = T, adjust.for.con.comp = T),
                family = 'gaussian',
                data = sim.admin2.dir,
                scale = 1/sim.admin2.dir$var.truth,
                control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                control.compute = list(config = TRUE))
  samples1a <- get_inla_samples_local(fit1a,1000)
  
  eta.samples1a <- samples1a[,2*n_admin2+1] + samples1a[,1:n_admin2]
  
  fit1 <- inla(mean.sim ~ f(admin3,model='bym2',hyper = hyperpc.bym2, graph = admin3.mat,
                            scale.model = T, adjust.for.con.comp = T),
               family = 'gaussian',
               data = sim.admin3.dir.stable,
               scale = 1/sim.admin3.dir.stable$var.truth,
               control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
               control.compute = list(config = TRUE))
  samples1 <- get_inla_samples_local(fit1,1000)
  
  eta.samples1 <- samples1[,2*n_admin3+1] + samples1[,1:n_admin3]
  
  ## model is specified with correct variance, but assuming normality ----
  fit2a <- inla(mean ~ f(admin2,model='bym2',hyper = hyperpc.bym2, graph = admin2.mat,
                         scale.model = T, adjust.for.con.comp = T),
                family = 'gaussian',
                data = sim.admin2.dir,
                scale = 1/sim.admin2.dir$var.truth,
                control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                control.compute = list(config = TRUE))
  samples2a <- get_inla_samples_local(fit2a,1000)
  
  eta.samples2a <- samples2a[,2*n_admin2+1] + samples2a[,1:n_admin2]
  
  fit2 <- inla(mean ~ f(admin3,model='bym2',hyper = hyperpc.bym2, graph = admin3.mat,
                        scale.model = T, adjust.for.con.comp = T),
               family = 'gaussian',
               data = sim.admin3.dir.stable,
               scale = 1/sim.admin3.dir.stable$var.truth,
               control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
               control.compute = list(config = TRUE))
  samples2 <- get_inla_samples_local(fit2,1000)
  
  eta.samples2 <- samples2[,2*n_admin3+1] + samples2[,1:n_admin3]
  
  ## real way of doing FH (estimated variance, assuming normality) ----
  fit3a <- inla(mean ~ f(admin2,model='bym2',hyper = hyperpc.bym2, graph = admin2.mat,
                         scale.model = T, adjust.for.con.comp = T),
                family = 'gaussian',
                data = sim.admin2.dir,
                scale = 1/sim.admin2.dir$variance,
                control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                control.compute = list(config = TRUE))
  samples3a <- get_inla_samples_local(fit3a,1000)
  
  eta.samples3a <- samples3a[,2*n_admin2+1] + samples3a[,1:n_admin2]
  
  fit3 <- inla(mean ~ f(admin3,model='bym2',hyper = hyperpc.bym2, graph = admin3.mat,
                        scale.model = T, adjust.for.con.comp = T),
               family = 'gaussian',
               data = sim.admin3.dir.stable,
               scale = 1/sim.admin3.dir.stable$variance,
               control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
               control.compute = list(config = TRUE))
  samples3 <- get_inla_samples_local(fit3,1000)
  
  eta.samples3 <- samples3[,2*n_admin3+1] + samples3[,1:n_admin3]
  
  ## record estimates -----
  
  results.adm3[[k]] <- merge(sim.admin3.dir.stable[,c('admin3','mean','mean.sim','variance','var.truth')],
                             data.frame(admin3 = 1:n_admin3,
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
                                        true_mean = admin3.means$V1),all=T)
  
  results.adm2[[k]] <- merge(sim.admin2.dir[,c('admin2','mean','mean.sim','variance','var.truth')],
                             data.frame(admin2 = 1:n_admin2,
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
                                        true_mean = admin2.means$V1),all=T)
  
}

## Diagnostics ------

# how many valid estimates does each admin3 have
colSums(do.call(rbind,lapply(direct.adm3.estimates,function(x){(!is.na(x$mean))})))

colSums(do.call(rbind,lapply(direct.adm3.estimates,function(x){(!is.na(x$variance) & x$variance > 1e-10)})))

# look at distribution of direct estimates for each area
for(i in 1:n_admin2){
  hist(unlist(lapply(direct.adm2.estimates, function(x){x$mean[i]})),prob=T)
  lines(seq(-2,2,0.01),dnorm(seq(-2,2,0.01),admin2.means$V1[i],sqrt(adm2.dir.variance[i])))
}
for(i in 1:n_admin3){
  hist(unlist(lapply(direct.adm3.estimates, function(x){x$mean[i]})),prob=T)
  lines(seq(-2,2,0.01),dnorm(seq(-2,2,0.01),admin3.means$V1[i],sqrt(adm3.dir.variance[i])))
}

# distribution of direct estimates for each area, standardized with true mean and variance
for(i in 1:n_admin2){
  hist(unlist(lapply(direct.adm2.estimates, function(x){(x$mean[i] - admin2.means$V1[i])/sqrt(adm2.dir.variance[i])})),prob=T)
}

# distribution of direct estimates for each area, standardized with true mean and estimated variance
for(i in 1:n_admin2){
  hist(unlist(lapply(direct.adm2.estimates, function(x){(x$mean[i] - admin2.means$V1[i])/sqrt(x$var[i])})),prob=T)
}

# distribution of direct estimates, standardized with true mean and estimated variance
hist((res.adm2.mat$mean - res.adm2.mat$true_mean)/sqrt(res.adm2.mat$variance),prob=T)
# distribution of direct estimates, standardized with true mean and variance
hist((res.adm2.mat$mean - res.adm2.mat$true_mean)/sqrt(res.adm2.mat$var.truth),prob=T)
lines(seq(-20,20,0.1),dnorm(seq(-20,20,0.1)))








## POST PROC ##################

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")

save(sim_pop_long,file='sim8_popdata.rda')
save(cluster_frame,file='sim8_clusterframe.rda')
save(cluster_alloc,file='sim8_clusteralloc.rda')
save(direct.adm2.estimates,file='sim8_admin2.direct.rda')
save(direct.adm3.estimates,file='sim8_admin3.direct.rda')
save(results.adm2,file = 'sim8_admin2.fh.rda')
save(results.adm3,file = 'sim8_admin3.fh.rda')

# load(file='sim8_popdata.rda')
# load(file='sim8_clusterframe.rda')
# load(file='sim8_clusteralloc.rda')
# load(file='sim8_admin2.direct.rda')
# load(file='sim8_admin3.direct.rda')
# load(file = 'sim8_admin2.fh.rda')
# load(file = 'sim8_admin3.fh.rda')

res.mat <- do.call(rbind,results.adm3)
res.adm2.mat <- do.call(rbind,results.adm2)


# rmse
rmse <- res.mat %>% group_by(admin3,true_mean) %>% reframe(rmse1 = sqrt(mean((mean.fh1 - true_mean)^2)),
                                                           rmse2 = sqrt(mean((mean.fh2 - true_mean)^2)),
                                                           rmse3 = sqrt(mean((mean.fh3 - true_mean)^2)))

rmse.adm2 <- res.adm2.mat %>% group_by(admin2,true_mean) %>% reframe(rmse1 = sqrt(mean((mean.fh1 - true_mean)^2)),
                                                                     rmse2 = sqrt(mean((mean.fh2 - true_mean)^2)),
                                                                     rmse3 = sqrt(mean((mean.fh3 - true_mean)^2)))


# coverage (area-specific and overall)
coverage <- res.mat %>% dplyr::mutate(cov1 = lower.fh1<true_mean & upper.fh1>true_mean,
                                      cov2 = lower.fh2<true_mean & upper.fh2>true_mean,
                                      cov3 = lower.fh3<true_mean & upper.fh3>true_mean) %>% 
  group_by(admin3,true_mean) %>% reframe(area_coverage1 = mean(cov1),
                                         area_coverage2 = mean(cov2),
                                         area_coverage3 = mean(cov3))


coverage.adm2 <- res.adm2.mat %>%
  mutate(cov1 = lower.fh1<true_mean & upper.fh1>true_mean,
         cov2 = lower.fh2<true_mean & upper.fh2>true_mean,
         cov3 = lower.fh3<true_mean & upper.fh3>true_mean) %>%
  group_by(admin2,true_mean) %>% summarise(area_coverage1 = mean(cov1),
                                           area_coverage2 = mean(cov2),
                                           area_coverage3 = mean(cov3))


colMeans(coverage)[c(3,4,5)]
colMeans(coverage.adm2)[c(3,4,5)]


## PLOTS ####################

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")

pdf('sim8_admin3_plots.pdf')
{
  par(mfrow=c(3,1))
  plot(coverage$true_mean,coverage$area_coverage1,
       xlab = 'Admin 3 (in ascending order of true rate)', ylab='',
       main = 'Coverage of 90% credible interval \n (correct sampling model)')
  abline(h=0.9,col='red')
  legend("topleft", legend = paste0("Coverage = ",100*round(mean(coverage$area_coverage1),4)), bty = "o")
  
  plot(coverage$true_mean,coverage$area_coverage2,
       xlab = 'Admin 3 (in ascending order of true rate)', ylab='',
       main = 'Coverage of 90% credible interval \n (correct V + assuming normality)')
  abline(h=0.9,col='red')
  legend("topleft", legend = paste0("Coverage = ",100*round(mean(coverage$area_coverage2),4)), bty = "o")
  
  plot(coverage$true_mean,coverage$area_coverage3,
       xlab = 'Admin 3 (in ascending order of true rate)', ylab='',
       main = 'Coverage of 90% credible interval \n (estimating V + assuming normality)')
  abline(h=0.9,col='red')
  legend("topleft", legend = paste0("Coverage = ",100*round(mean(coverage$area_coverage3),4)), bty = "o")
  
  plot(coverage$area_coverage2,coverage$area_coverage1,
       main = 'Area coverage of 90% credible interval (Admin3)',
       xlab = 'correct V + assuming normality',
       ylab = 'correct sampling model')
  abline(0,1)
  plot(coverage$area_coverage3,coverage$area_coverage1,
       main = 'Area coverage of 90% credible interval (Admin3)',
       xlab = 'estimating V + assuming normality',
       ylab = 'correct sampling model')
  abline(0,1)
  plot(coverage$area_coverage3,coverage$area_coverage2,
       main = 'Area coverage of 90% credible interval (Admin3)',
       xlab = 'estimating V + assuming normality',
       ylab = 'correct V + assuming normality')
  abline(0,1)
  
  #rmse
  plot(rmse$true_mean,rmse$rmse1,
       xlab = 'Admin 3 (in ascending order of true rate)', ylab='',
       main = 'RMSE (correct sampling model)')
  legend("topleft", legend = paste0("Average = ",round(mean(rmse$rmse1),3)), bty = "o")
  
  plot(rmse$true_mean,rmse$rmse2,
       xlab = 'Admin 3 (in ascending order of true rate)', ylab='',
       main = 'RMSE (correct V + assuming normality)')
  legend("topleft", legend = paste0("Average = ",round(mean(rmse$rmse2),3)), bty = "o")
  
  plot(rmse$true_mean,rmse$rmse3,
       xlab = 'Admin 3 (in ascending order of true rate)', ylab='',
       main = 'RMSE (estimating V + assuming normality)')
  legend("topleft", legend = paste0("Average = ",round(mean(rmse$rmse3),3)), bty = "o")
  
  plot(rmse$rmse2,rmse$rmse1,
       main = 'Root Mean Squared Error (Admin3)',
       xlab = 'correct V + assuming normality',
       ylab = 'correct sampling model')
  abline(0,1)
  plot(rmse$rmse3,rmse$rmse1,
       main = 'Root Mean Squared Error (Admin3)',
       xlab = 'estimating V + assuming normality',
       ylab = 'correct sampling model')
  abline(0,1)
  plot(rmse$rmse3,rmse$rmse2,
       main = 'Root Mean Squared Error (Admin3)',
       xlab = 'estimating V + assuming normality',
       ylab = 'correct V + assuming normality')
  abline(0,1)
  
}
dev.off()

pdf('sim8_admin2_plots.pdf')
{
  par(mfrow=c(3,1))
  plot(coverage.adm2$true_mean,coverage.adm2$area_coverage1,
       xlab = 'Admin 2 (in ascending order of true rate)', ylab='',
       main = 'Coverage of 90% credible interval \n (correct sampling model)')
  abline(h=0.9,col='red')
  legend("topleft", legend = paste0("Coverage = ",100*round(mean(coverage.adm2$area_coverage1),4)), bty = "o")
  
  plot(coverage.adm2$true_mean,coverage.adm2$area_coverage2,
       xlab = 'Admin 2 (in ascending order of true rate)', ylab='',
       main = 'Coverage of 90% credible interval \n (correct V + assuming normality)')
  abline(h=0.9,col='red')
  legend("topleft", legend = paste0("Coverage = ",100*round(mean(coverage.adm2$area_coverage2),4)), bty = "o")
  
  plot(coverage.adm2$true_mean,coverage.adm2$area_coverage3,
       xlab = 'Admin 32 (in ascending order of true rate)', ylab='',
       main = 'Coverage of 90% credible interval \n (estimating V + assuming normality)')
  abline(h=0.9,col='red')
  legend("topleft", legend = paste0("Coverage = ",100*round(mean(coverage.adm2$area_coverage3),4)), bty = "o")
  
  plot(coverage.adm2$area_coverage2,coverage.adm2$area_coverage1,
       main = 'Area coverage of 90% credible interval (Admin2)',
       xlab = 'correct V + assuming normality',
       ylab = 'correct sampling model')
  abline(0,1)
  plot(coverage.adm2$area_coverage3,coverage.adm2$area_coverage1,
       main = 'Area coverage of 90% credible interval (Admin2)',
       xlab = 'estimating V + assuming normality',
       ylab = 'correct sampling model')
  abline(0,1)
  plot(coverage.adm2$area_coverage3,coverage.adm2$area_coverage2,
       main = 'Area coverage of 90% credible interval (Admin2)',
       xlab = 'estimating V + assuming normality',
       ylab = 'correct V + assuming normality')
  abline(0,1)
  
  #rmse
  plot(rmse.adm2$true_mean,rmse.adm2$rmse1,
       xlab = 'Admin 2 (in ascending order of true rate)', ylab='',
       main = 'RMSE (correct sampling model)')
  legend("topleft", legend = paste0("Average = ",round(mean(rmse.adm2$rmse1),3)), bty = "o")
  
  plot(rmse.adm2$true_mean,rmse.adm2$rmse2,
       xlab = 'Admin 2 (in ascending order of true rate)', ylab='',
       main = 'RMSE (correct V + assuming normality)')
  legend("topleft", legend = paste0("Average = ",round(mean(rmse.adm2$rmse2),3)), bty = "o")
  
  plot(rmse.adm2$true_mean,rmse.adm2$rmse3,
       xlab = 'Admin 2 (in ascending order of true rate)', ylab='',
       main = 'RMSE (estimating V + assuming normality)')
  legend("topleft", legend = paste0("Average = ",round(mean(rmse.adm2$rmse3),3)), bty = "o")
  
  plot(rmse.adm2$rmse2,rmse.adm2$rmse1,
       main = 'Root Mean Squared Error (Admin2)',
       xlab = 'correct V + assuming normality',
       ylab = 'correct sampling model')
  abline(0,1)
  plot(rmse.adm2$rmse3,rmse.adm2$rmse1,
       main = 'Root Mean Squared Error (Admin2)',
       xlab = 'estimating V + assuming normality',
       ylab = 'correct sampling model')
  abline(0,1)
  plot(rmse.adm2$rmse3,rmse.adm2$rmse2,
       main = 'Root Mean Squared Error (Admin2)',
       xlab = 'estimating V + assuming normality',
       ylab = 'correct V + assuming normality')
  abline(0,1)
  
}
dev.off()









