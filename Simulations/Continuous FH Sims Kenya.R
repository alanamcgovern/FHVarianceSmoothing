library(sf)
library(tidyverse)
library(INLA)
library(ggpubr)
library(locfit)
library(surveyPrev)
library(raster)
library(spdep)
library(exactextractr)
library(data.table)
library(terra)
library(data.table)
library(cmdstanr)

source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")

# load geometry (country specific) ------
setwd('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm')
poly.adm2 <- st_read(dsn = 'gadm41_KEN_shp', layer = "gadm41_KEN_2", options = "ENCODING=UTF-8")
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

# load population (country specific) -------

pop_dens <- rast("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Kenya/Population/ken_u5_2022_1km.tif")

poly.adm2$admin_pop <- round(exact_extract(pop_dens, poly.adm2, 'sum'))

pop.admin2 <- merge(poly.adm2[,c('admin2.char','admin_pop')],admin.key)
pop.admin2 <- pop.admin2[order(pop.admin2$admin2),]

pop.admin1 <- pop.admin2 %>% group_by(admin1) %>% summarise(admin_pop = sum(admin_pop))

A_2to1 <- matrix(0,n_admin1,n_admin2)
for(area1 in 1:n_admin1){
  which.areas2 <- unique(admin.key[admin.key$admin1==area1,]$admin2)
  A_2to1[area1,which.areas2] <- pop.admin2[pop.admin2$admin2 %in% which.areas2,]$admin_pop/sum(pop.admin2[pop.admin2$admin2 %in% which.areas2,]$admin_pop)
}


# sampling info (survey specific) -------------

# total number of EAs in each strata
ea_totals <- data.table(
  NAME_1 = rep(c(
    "Mombasa","Kwale","Kilifi","Tana River","Lamu","Taita Taveta","Garissa","Wajir","Mandera","Marsabit",
    "Isiolo","Meru","Tharaka-Nithi","Embu","Kitui","Machakos","Makueni","Nyandarua","Nyeri","Kirinyaga",
    "Murang'a","Kiambu","Turkana","West Pokot","Samburu","Trans Nzoia","Uasin Gishu","Elgeyo-Marakwet",
    "Nandi","Baringo","Laikipia","Nakuru","Narok","Kajiado","Kericho","Bomet","Kakamega","Vihiga",
    "Bungoma","Busia","Siaya","Kisumu","Homa Bay","Migori","Kisii","Nyamira","Nairobi"
  ),times=2),
  urban = c(rep(1,n_admin1),rep(0,n_admin1)),
  EAs = c(4031,365,1141,226,99,299,441,320,386,215, #urban
          273,492,125,263,224,1413,297,242,560,457,
          406,6012,351,95,138,463,1863,74,
          201,291,476,3431,375,2272,336,115,642,160,
          557,336,280,1222,362,461,558,157,15621,
          #rural
          0,1462,1971,607,250,710,871,977,901,651,
          232,3771,1109,1555,3507,2720,2481,1637,1960,1490,
          2793,2085,1343,1486,559,1770,1492,1089,
          1893,1606,1076,2703,2245,1094,1829,1818,3963,1241,
          3175,1727,2396,1711,2526,2124,2968,1541,0))

ea_totals <- merge(ea_totals, unique.array(admin.key[,c('admin1','admin1.char','NAME_1')]))

# average number of individuals to sample from each cluster 
n_per_EA <- 10

# number of EAs sampled from each strata

cluster_alloc <- data.table(
  NAME_1 = rep(c(
    "Mombasa","Kwale","Kilifi","Tana River","Lamu","Taita Taveta","Garissa","Wajir","Mandera","Marsabit",
    "Isiolo","Meru","Tharaka-Nithi","Embu","Kitui","Machakos","Makueni","Nyandarua","Nyeri","Kirinyaga",
    "Murang'a","Kiambu","Turkana","West Pokot","Samburu","Trans Nzoia","Uasin Gishu","Elgeyo-Marakwet",
    "Nandi","Baringo","Laikipia","Nakuru","Narok","Kajiado","Kericho","Bomet","Kakamega","Vihiga",
    "Bungoma","Busia","Siaya","Kisumu","Homa Bay","Migori","Kisii","Nyamira","Nairobi"
  ),times=2),
  urban = c(rep(1,n_admin1),rep(0,n_admin1)),
  EAs = c(
    #urban
    39,12,17,20,13,14,13,13,14,13,
    17,11,9,11,8,17,10,10,13,13,
    11,27,12,8,11,13,20,7,
    9,11,14,22,11,23,11,7,11,9,
    11,11,9,18,10,12,11,9,58,
    # rural
    0,22,20,13,19,19,21,21,20,20,
    16,27,24,24,28,21,26,25,23,22,
    26,15,23,26,22,23,17,26,
    26,23,20,18,25,15,25,28,27,25,
    26,24,27,19,26,24,26,26,0)
)

cluster_alloc <- merge(cluster_alloc, unique.array(admin.key[,c('admin1','admin1.char','NAME_1')]))


## urban is oversampled in most (90%) areas
merge(ea_totals %>% group_by(NAME_1) %>% summarise(urb_frac_pop = sum(urban*EAs)/sum(EAs)),cluster_alloc %>% 
        group_by(NAME_1) %>% summarise(urb_frac_sample = sum(urban*EAs)/sum(EAs))) %>% 
  ggplot() + geom_point(aes(urb_frac_pop,urb_frac_sample)) + geom_abline(intercept = 0,slope = 1)

# sort population into clusters -----------
cluster_frame <- ea_totals[rep(1:.N, EAs)]
cluster_frame[, cluster := seq_len(.N)]
cluster_frame[,EAs := NULL]
cluster_frame[,strata := admin1 + urban*n_admin1]

# how many EAs for each admin2?
EA_props <- lapply(1:n_admin1,function(i){pop.admin2[pop.admin2$admin1==i,]$admin_pop/sum(pop.admin2[pop.admin2$admin1==i,]$admin_pop)})
ids_list <- lapply(1:n_admin1,function(i){pop.admin2[pop.admin2$admin1==i,]$admin2})

cluster_frame[, admin2 := {
  props <- EA_props[[admin1]]           # get vector of proportions for this state
  ids <- ids_list[[admin1]]
  
  n <- .N                                  # number of EAs in this state
  group_sizes <- round(props * n)          # convert proportions to group sizes
  # adjust if rounding doesn’t sum exactly to n
  diff_n <- n - sum(group_sizes)
  if (diff_n != 0) group_sizes[1:abs(diff_n)] <- group_sizes[1:abs(diff_n)] + sign(diff_n)
  
  idx <- sample(seq_len(n))
  rep(ids, times = group_sizes)[order(idx)]
}, by = admin1]

cluster_frame[,admin2_strata := admin2 + urban*n_admin2]

sim_pop <- as.data.table(pop.admin2)
sim_pop_long <- sim_pop[rep(seq_len(n_admin2),times = sim_pop$admin_pop)]

# assign each unit to a cluster
EA_list <- lapply(1:n_admin2,function(i){cluster_frame[cluster_frame$admin2==i,]$cluster})

sim_pop_long[, cluster := {
  vals <- EA_list[[admin2]]
  sample(vals, .N, replace = TRUE)
}, by = admin2]

sim_pop_long[,admin_pop := NULL]
sim_pop_long[,geometry := NULL]
sim_pop_long <- merge(cluster_frame[,c('cluster','urban','strata','admin2_strata')],sim_pop_long)

cluster_frame <- merge(cluster_frame,sim_pop_long[,.N,cluster],by='cluster')

# objects for sampling ---------
cluster_frame$n <-  rpois(nrow(cluster_frame),n_per_EA) #number that will be sampled from that cluster if sampled

# get population urban rural percent for weights
urban_df <- sim_pop_long[,mean(urban),by=NAME_1]
pop.admin1 <- merge(pop.admin1,merge(urban_df, unique.array(admin.key[,c('admin1','admin1.char','NAME_1')])))
cluster_alloc <- cluster_alloc %>% mutate(strata_total =  pop.admin1$admin_pop[admin1]*(urban*pop.admin1$V1[admin1] + (1-urban)*(1-pop.admin1$V1[admin1])))

cluster_frame <- merge(cluster_frame,cluster_alloc,by=c('admin1','NAME_1','urban','admin1.char')) %>% mutate(wt = strata_total/(EAs*n))
cluster_frame[,EAs:=NULL]
cluster_frame[,strata_total:=NULL]

cluster_alloc$strata <- cluster_alloc$admin1 + cluster_alloc$urban*n_admin1

# helper functions and objects ----------
options(survey.lonely.psu = "adjust")

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

fit_inla_with_samples <- function(outcome,
                                  scale_var,
                                  data,
                                  random_effect = "admin1",
                                  random_model = c("iid", "bym2"),
                                  adj_matrix = NULL,      # only needed for bym2
                                  covariates = NULL,
                                  n_samples = 1000,
                                  hyperpc.iid = list(prec = list(prior = "pc.prec", param = c(1, 0.01))),
                                  hyperpc.bym2 = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                                                      phi = list(prior = "pc", param = c(0.5, 0.5)))) {
  
  # ---- Build formula ----
  covariate_part <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    "1"
  }
  
  if (random_model == "iid") {
    random_term <- paste0("f(", random_effect, ", model='iid', hyper=hyperpc.iid)")
  } else if (random_model == "bym2") {
    random_term <- paste0("f(", random_effect,
                          ", model='bym2', graph=adj_matrix, hyper=hyperpc.bym2, scale.model=T, adjust.for.con.comp=T)")
  }
  
  formula_str <- paste0(outcome, " ~ ", covariate_part, " + ", random_term)
  formula <- as.formula(formula_str)
  
  # ---- Fit INLA model ----
  fit <- inla(
    formula,
    family = "gaussian",
    data = data,
    scale = 1 / data[[scale_var]],
    control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
    control.compute = list(config = TRUE)
  )
  
  # ---- Get posterior samples ----
  samples <- get_inla_samples_local(fit, n_samples)
  
  # ---- Build model matrix for fixed effects ----
  X <- model.matrix(as.formula(paste0("~ ", covariate_part)), data = data)
  n_fixed <- ncol(X)
  
  # ---- Figure out parameter layout ----
  n_random <- length(unique(data[[random_effect]]))
  if(random_model == 'iid'){
    intercept_col <- n_random + 1
  }else if(random_model =='bym2'){
    intercept_col <- 2*n_random + 1
  }
  
  fixed_cols <- if (n_fixed > 1) {
    (intercept_col + 1):(intercept_col + (n_fixed - 1))
  } else NULL
  
  eta_samples <- samples[, intercept_col] + samples[, 1:n_random]
  
  # Add fixed effects contribution
  if (!is.null(fixed_cols)) {
    eta_samples <- eta_samples + samples[, fixed_cols] %*% t(X[, -1, drop = FALSE])
  }
  
  return(eta_samples)
}

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
satt_mod <- cmdstan_model("IID_FH_VarSmooth_Satt.stan")

# load covariates ---------
setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")

# ur_weights_adm1 <- readRDS('Kenya_2022_admin1_u5_ur_weights.rds')
# ur_weights_adm2 <- readRDS('Kenya_2022_admin2_u5_ur_weights.rds')

setwd("/Users/alanamcgovern/Desktop/Research/KEN_Covariates")
load(file='Kenya_admin2_covariates.rda')
load(file='Kenya_admin1_covariates.rda')

cmat_admin1 <- merge(cmat_admin1,sim_pop_long[order(admin1),mean(urban),by=admin1]) %>% rename(urb_frac = V1)
cmat_admin2 <- merge(cmat_admin2,sim_pop_long[order(admin2),mean(urban),by=admin2]) %>% rename(urb_frac = V1)

# mean covariates: density_log, nt_lights_log, tthc_log, precip, temp, elev, urb_frac
# var covariates: pop_var_log, nt_lights_var, tthc_var, precip_var, temp_var, elev_var, area_log

# set simulation parameters ---------
admin_level <- c(1,2)[2]

load(file=paste0('/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/ken_2022_cont_outcomes.rda'))
dir.dat <- dat[!is.na(dat[,c('haz')]),]
dir.dat$value <- dir.dat[,c('haz')]

my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"), 
                                  strata = ~v023, nest = T, weights = ~v005, data = dir.dat)

admin1 <- NULL
v025 <- NULL
x <- lapply(1:n_admin1,function(which.area){
  
  tmp <- subset(my.svydesign, (admin1 == as.character(which.area) & v025==1))
  
  if (dim(tmp)[1] == 0) {
    out <- c(which.area,1,rep(NA, 2))
  } else {
    lm.ob <- survey::svymean(value ~ 1, design = tmp)
    out <- c(which.area,1, lm.ob[1], vcov(lm.ob))
  }
  
  tmp <- subset(my.svydesign, (admin1 == as.character(which.area) & v025==2))
  
  if (dim(tmp)[1] == 0) {
    out <- rbind(out,c(which.area,0,rep(NA, 2)))
  } else {
    lm.ob <- survey::svymean(value ~ 1, design = tmp)
    out <- rbind(out,c(which.area,0, lm.ob[1], vcov(lm.ob)))
  }
  
})

admin1.strata.dir <- as.data.frame(do.call(rbind,x))
colnames(admin1.strata.dir) <- c('admin1','urban','mean','variance')
admin1.strata.dir$strata <- admin1.strata.dir$admin1 + admin1.strata.dir$urban*n_admin1

if(admin_level==2){
  admin2 <- NULL
  v025 <- NULL
  x <- lapply(1:n_admin2,function(which.area){
    
    tmp <- subset(my.svydesign, (admin2 == as.character(which.area) & v025==1))
    
    if (dim(tmp)[1] == 0) {
      out <- c(which.area,1,rep(NA, 2))
    } else {
      lm.ob <- survey::svymean(value ~ 1, design = tmp)
      out <- c(which.area,1, lm.ob[1], vcov(lm.ob))
    }
    
    tmp <- subset(my.svydesign, (admin2 == as.character(which.area) & v025==2))
    
    if (dim(tmp)[1] == 0) {
      out <- rbind(out,c(which.area,0,rep(NA, 2)))
    } else {
      lm.ob <- survey::svymean(value ~ 1, design = tmp)
      out <- rbind(out,c(which.area,0, lm.ob[1], vcov(lm.ob)))
    }
    return(out)
    
  })
  
  admin2.strata.dir <- as.data.frame(do.call(rbind,x))
  colnames(admin2.strata.dir) <- c('admin2','urban','mean','variance')
  admin2.strata.dir$admin2_strata <- admin2.strata.dir$admin2 + admin2.strata.dir$urban*n_admin2
  
  missing_strata <- which(is.na(admin2.strata.dir$mean))
  
  for(j in missing_strata){
    #which larger strata is a part of?
    urb_tmp <- admin2.strata.dir[j,]$urban
    admin1_tmp <- admin.key[admin.key$admin2==admin2.strata.dir[j,]$admin2,]$admin1
    
    # don't need to add rural means for areas that are all urban
    if(admin1_tmp %in% c(28,30) & urb_tmp==0){
      admin2.strata.dir[j,]$mean <- NA
    }else{
      #use the larger strata mean as the mean
      mean_tmp <- admin1.strata.dir[admin1.strata.dir$admin1==admin1_tmp & admin1.strata.dir$urban==urb_tmp,]$mean
      # and the variance of within-strata means as the variance
      sd_tmp <- as.numeric(admin2.strata.dir %>% filter(admin2 %in% admin.key[admin.key$admin1==admin1_tmp,]$admin2, urban == urb_tmp) %>% summarise(sd(mean,na.rm=T)))
      if(is.na(sd_tmp)){
        sd_tmp <- as.numeric(admin2.strata.dir %>% filter(admin2 %in% admin.key[admin.key$admin1==admin1_tmp,]$admin2) %>% summarise(sd(mean,na.rm=T)))
      }
      admin2.strata.dir[j,]$mean <- rnorm(1,mean_tmp,sd_tmp)
    }
  }  
}

alpha <- 0.5  #dir.dat %>% group_by(admin1) %>% summarise(v = log(var(value))) %>% summarise(mean(v))
tau <- 30 #dir.dat %>% group_by(admin1) %>% summarise(v = log(var(value))) %>% summarise(1/var(v))

kappa <- rnorm(n_admin1,1,0.25) # always at admin1 level (or at least that's how we will fit the model)

# simulate data at population level (country/data specific) -----

sig2_R <- exp(rnorm(n_admin1,alpha,1/tau))
sig2 <- c(sig2_R, kappa*sig2_R)

if(admin_level==1){
  admin1.strata.dir <- admin1.strata.dir[order(admin1.strata.dir$strata),]
  sim_pop_long[, value := admin1.strata.dir$mean[strata] + rnorm(.N,0,sqrt(sig2[strata]))]
}
if(admin_level==2){
  admin2.strata.dir <- admin2.strata.dir[order(admin2.strata.dir$admin2_strata),]
  sim_pop_long[, value := admin2.strata.dir$mean[admin2_strata] + rnorm(.N,0,sqrt(sig2[strata]))]
}

# is the indexing correct?
# plot(sim_pop_long[order(strata),mean(value),by=strata]$V1,admin1.strata.dir[admin1.strata.dir$strata %in% cluster_frame$strata,]$mean)
# abline(0,1)
# 
# plot(sim_pop_long[order(admin2_strata),mean(value),by=admin2_strata]$V1,admin2.strata.dir[admin2.strata.dir$admin2_strata %in% cluster_frame$admin2_strata,]$mean)
# abline(0,1)

 admin2_means <- (merge(admin2.strata.dir %>% dplyr::select(admin2,urban,mean) %>% 
             pivot_wider(names_from=urban,values_from=mean),cmat_admin2[,c('admin2','urb_frac')]) %>%
  mutate(mean = ifelse(urb_frac==1,`1`,`0`*(1-urb_frac) + `1`*urb_frac)))$mean

admin1_means <- A_2to1%*%admin2_means

# more checks
# plot(sim_pop_long[order(admin1),mean(value),by=admin1]$V1,admin1_means)
# plot(sim_pop_long[order(admin2),mean(value),by=admin2]$V1,admin2_means)

# sample using 2 stage sampling and get direct estimates -----------
nsim <- 100

direct.adm2 <- direct.adm1 <- sampled_clusters <- list()
for(k in 1:nsim){
  cat(k,'\n')
  
  # which clusters will we sample?
  cluster_sample_ids <- unlist(sapply( unique(cluster_frame$strata),function(i){
    sample(cluster_frame[cluster_frame$strata==i,]$cluster,
           cluster_alloc[cluster_alloc$strata==i,]$EAs,
           prob = cluster_frame[cluster_frame$strata==i,]$N)
  }))
  cluster_sample <- cluster_frame[cluster_frame$cluster %in% cluster_sample_ids,c('cluster','n','wt')]
  
  sim_sample <- sim_pop_long[
    cluster_sample,                   # right table
    on = "cluster"                   # join by cluster
  ][, .SD[sample(.N, min(.N, unique(n)))], by = cluster] # for each cluster, sample n
  
  sim_sample[,n:=NULL]
  
  ## save info (weights and sample sizes) on sampled cluster for Satt approx
  sampled_clusters[[k]] <- sim_sample %>% group_by(admin1,admin2,strata,cluster) %>% reframe(n=n(),wt=unique(wt))
  
  
  # get direct estimates (not accounting for stratification) --------
  options(survey.lonely.psu = "adjust")
  my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"), 
                                    strata = ~strata, nest = T, weights = ~wt, data = sim_sample)
  
  # ## admin2
  x <- mapply(admin2.HT.withNA, which.area = 1:n_admin2)
  admin2.dir <- data.frame(t(x))
  colnames(admin2.dir) <- c('admin2','mean','variance')
  sample_info <- merge(sim_sample[,length(unique(cluster)),by=admin2],sim_sample[,.N,by=admin2])
  colnames(sample_info)[2:3] <- c('n_clusters','n_obs')

  direct.adm2[[k]] <- merge(sample_info,admin2.dir,all=T)

  ## admin1
  x <- mapply(admin1.HT.withNA, which.area = 1:n_admin1)
  admin1.dir <- data.frame(t(x))
  colnames(admin1.dir) <- c('admin1','mean','variance')
  sample_info <- merge(sim_sample[,length(unique(cluster)),by=admin1],sim_sample[,.N,by=admin1])
  colnames(sample_info)[2:3] <- c('n_clusters','n_obs')
  
  direct.adm1[[k]] <- merge(sample_info,admin1.dir,all=T)
  
}

adm2.dir.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct.adm2, function(x){x$mean[i]})),na.rm = T)})
adm1.dir.variance <- sapply(1:n_admin1, function(i){var(unlist(lapply(direct.adm1, function(x){x$mean[i]})),na.rm = T)})

# diagnostics ----------
par(mfrow=c(3,3))
for(area in 1:n_admin1){
  plot(sapply(1:nsim, function(i){var(unlist(lapply(direct.adm1[1:i], function(x){x$mean[area]})),na.rm = T)}),type='l',main = paste('Area ', area))
}

for(area in 1:n_admin1){
  hist(unlist(lapply(direct.adm1, function(x){x$var[area]})))
  abline(v=adm1.dir.variance[area],col='red')
}

for(area in sample(1:n_admin2,18)){
  hist(unlist(lapply(direct.adm2, function(x){x$var[area]})))
  abline(v=adm2.dir.variance[area],col='red')
}

# save direct estimates and other simulation info -------------

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
# 
params_list <- list(admin2.strata = admin2.strata.dir,
                    alpha = alpha,
                    tau = tau,
                    kappa = kappa,
                    sig2 = sig2)
save(params_list,file='simX_params.rda')
save(sim_pop_long,file='simX_popdata.rda')
save(cluster_frame,file='simX_clusterframe.rda')
save(cluster_alloc,file='simX_clusteralloc.rda')
save(sampled_clusters,file='simX_sampled_clusters.rda')
save(direct.adm1,file='simX_admin1.direct.rda')
save(direct.adm2,file='simX_admin2.direct.rda')

# get FH estimates (oracle and no smoothing) from direct --------

results.adm1 <- results.adm2 <- list() # record estimates (direct and fay herriot)
#covariates = mean_covariates
covariates = c('temp','elev','tthc_log','nt_lights_log')
random_model = c("iid","bym2")[2]
for(k in 1:100){
  cat(k,'\n')
  
   sim.admin2.dir <- direct.adm2[[k]]
   sim.admin2.dir$var.truth <- adm2.dir.variance
   sim.admin2.dir <- sim.admin2.dir %>% mutate(variance = ifelse(variance < 1e-5,NA,variance))
   sim.admin2.dir <- sim.admin2.dir %>% mutate(mean = ifelse(is.na(variance),NA,mean),
                                               var.truth = ifelse(is.na(variance),NA,var.truth))
   sim.admin2.dir <- merge(sim.admin2.dir,cmat_admin2)
   sim.admin2.dir <- sim.admin2.dir[order(sim.admin2.dir$admin2),]
  
  sim.admin1.dir <- direct.adm1[[k]]
  sim.admin1.dir$var.truth <- adm1.dir.variance
  sim.admin1.dir <- merge(sim.admin1.dir,cmat_admin1)
  sim.admin1.dir <- sim.admin1.dir[order(sim.admin1.dir$admin1),]

  ## model is specified with correct variance ----
  eta.samples1a <- fit_inla_with_samples(outcome = 'mean',
                                         scale_var = 'var.truth',
                                         data = sim.admin1.dir,
                                         random_effect = "admin1",
                                         covariates = covariates,
                                         adj_matrix = admin1.mat,
                                         random_model = random_model)
  
  eta.samples1 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'var.truth',
                                        data = sim.admin2.dir,
                                        random_effect = "admin2",
                                        covariates = covariates,
                                        adj_matrix = admin2.mat,
                                        random_model = random_model)

  # naive FH (estimated variance, assuming normality) ----
  eta.samples2a <- fit_inla_with_samples(outcome = 'mean',
                                         scale_var = 'variance',
                                         data = sim.admin1.dir,
                                         random_effect = "admin1",
                                         covariates = covariates,
                                         adj_matrix = admin1.mat,
                                         random_model = random_model)
  
  eta.samples2 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = sim.admin2.dir,
                                        random_effect = "admin2",
                                        covariates = covariates,
                                        adj_matrix = admin2.mat,
                                        random_model = random_model)
  
  ## record estimates -----
  
  results.adm2[[k]] <- merge(sim.admin2.dir[,c('admin2','mean','variance','var.truth')],
                             data.frame(admin2 = 1:n_admin2,
                                        # uses correct variance
                                        mean.fh2  = apply(eta.samples2,2,mean),
                                        var.fh2 = apply(eta.samples2,2,var),
                                        lower.fh2 = apply(eta.samples2,2,quantile,prob=0.05),
                                        upper.fh2 = apply(eta.samples2,2,quantile,prob=0.95),
                                        # standard FH
                                        mean.fh1  = apply(eta.samples1,2,mean),
                                        var.fh1 = apply(eta.samples1,2,var),
                                        lower.fh1 = apply(eta.samples1,2,quantile,prob=0.05),
                                        upper.fh1 = apply(eta.samples1,2,quantile,prob=0.95),
                                        sim = k,
                                        true_mean = admin2_means),all=T)
  
  results.adm1[[k]] <- merge(sim.admin1.dir[,c('admin1','mean','variance','var.truth')],
                             data.frame(admin1 = 1:n_admin1,
                                        # uses correct variance
                                        mean.fh2  = apply(eta.samples2a,2,mean),
                                        var.fh2 = apply(eta.samples2a,2,var),
                                        lower.fh2 = apply(eta.samples2a,2,quantile,prob=0.05),
                                        upper.fh2 = apply(eta.samples2a,2,quantile,prob=0.95),
                                        # standard FH
                                        mean.fh1  = apply(eta.samples1a,2,mean),
                                        var.fh1 = apply(eta.samples1a,2,var),
                                        lower.fh1 = apply(eta.samples1a,2,quantile,prob=0.05),
                                        upper.fh1 = apply(eta.samples1a,2,quantile,prob=0.95),
                                        sim = k,
                                        true_mean = admin1_means),all=T)
  
}

# SAVE/LOAD RESULTS ---------

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
# 
save(results.adm1,file = 'sim8_complex_admin1.fh.rda')
save(results.adm2,file = 'sim8_complex_admin2.fh.rda')
#### POST PROC -----

res.adm1.mat <- do.call(rbind,results.adm1)
res.adm1.mat <- res.adm1.mat %>% mutate(lower.direct = mean - 1.645*sqrt(variance),
                                        upper.direct = mean + 1.645*sqrt(variance),
                                        ciwidth0 = upper.direct - lower.direct,
                                        ciwidth4 = upper.fh4 - lower.fh4,
                                        ciwidth2 = upper.fh2 - lower.fh2,
                                        ciwidth3 = upper.fh3 - lower.fh3,
                                        se0 = (true_mean - mean)^2,
                                        se4 = (true_mean - mean.fh4)^2,
                                        se2 = (true_mean - mean.fh2)^2,
                                        se3 = (true_mean - mean.fh3)^2,
                                        cov0 = lower.direct<true_mean & upper.direct>true_mean,
                                        cov4 = lower.fh4<true_mean & upper.fh4>true_mean,
                                        cov2 = lower.fh2<true_mean & upper.fh2>true_mean,
                                        cov3 = lower.fh3<true_mean & upper.fh3>true_mean,
                                        i.score0 = ciwidth0 + 20*(lower.direct-true_mean)*(true_mean<lower.direct) + 20*(true_mean-upper.direct)*(true_mean>upper.direct),
                                        i.score4 = ciwidth4 + 20*(lower.fh4-true_mean)*(true_mean<lower.fh4) + 20*(true_mean-upper.fh4)*(true_mean>upper.fh4),
                                        i.score2 = ciwidth2 + 20*(lower.fh2-true_mean)*(true_mean<lower.fh2) + 20*(true_mean-upper.fh2)*(true_mean>upper.fh2),
                                        i.score3 = ciwidth3 + 20*(lower.fh3-true_mean)*(true_mean<lower.fh3) + 20*(true_mean-upper.fh3)*(true_mean>upper.fh3))

res.adm2.mat <- do.call(rbind,results.adm2)
res.adm2.mat <- res.adm2.mat %>% mutate(lower.direct = mean - 1.645*sqrt(variance),
                                        upper.direct = mean + 1.645*sqrt(variance),
                                        ciwidth0 = upper.direct - lower.direct,
                                         ciwidth4 = upper.fh4 - lower.fh4,
                                        ciwidth2 = upper.fh2 - lower.fh2,
                                        ciwidth3 = upper.fh3 - lower.fh3,
                                        se0 = (true_mean - mean)^2,
                                          se4 = (true_mean - mean.fh4)^2,
                                        se2 = (true_mean - mean.fh2)^2,
                                        se3 = (true_mean - mean.fh3)^2,
                                        cov0 = lower.direct<true_mean & upper.direct>true_mean,
                                         cov4 = lower.fh4<true_mean & upper.fh4>true_mean,
                                        cov2 = lower.fh2<true_mean & upper.fh2>true_mean,
                                        cov3 = lower.fh3<true_mean & upper.fh3>true_mean,
                                        i.score0 = ciwidth0 + 20*(lower.direct-true_mean)*(true_mean<lower.direct) + 20*(true_mean-upper.direct)*(true_mean>upper.direct),
                                         i.score4 = ciwidth4 + 20*(lower.fh4-true_mean)*(true_mean<lower.fh4) + 20*(true_mean-upper.fh4)*(true_mean>upper.fh4),
                                        i.score2 = ciwidth2 + 20*(lower.fh2-true_mean)*(true_mean<lower.fh2) + 20*(true_mean-upper.fh2)*(true_mean>upper.fh2),
                                        i.score3 = ciwidth3 + 20*(lower.fh3-true_mean)*(true_mean<lower.fh3) + 20*(true_mean-upper.fh3)*(true_mean>upper.fh3))

coverage.adm1 <- res.adm1.mat %>% group_by(admin1,true_mean) %>% reframe(area_coverage0 = mean(cov0,na.rm=T),
                                                                         area_coverage4 = mean(cov4),
                                                                         area_coverage2 = mean(cov2),
                                                                         area_coverage3 = mean(cov3))


coverage.adm2 <- res.adm2.mat %>% group_by(admin2,true_mean) %>% reframe(area_coverage0 = mean(cov0,na.rm=T),
                                                                         area_coverage4 = mean(cov4),
                                                                         area_coverage2 = mean(cov2),
                                                                         area_coverage3 = mean(cov3))

coverage.adm2.insample <- res.adm2.mat %>% filter(!is.na(mean)) %>%
  group_by(admin2,true_mean) %>% reframe(area_coverage0 = mean(cov0),
                                         area_coverage4 = mean(cov4),
                                         area_coverage2 = mean(cov2),
                                         area_coverage3 = mean(cov3))

rmse.adm1 <- res.adm1.mat %>% group_by(admin1,true_mean) %>% reframe(rmse0 = sqrt(mean(se0)),
                                                                     rmse4 = sqrt(mean(se4)),
                                                                     rmse2 = sqrt(mean(se2)),
                                                                     rmse3 = sqrt(mean(se3)))


rmse.adm2 <- res.adm2.mat %>% group_by(admin2,true_mean) %>% reframe(rmse4 = sqrt(mean(se4)),
                                                        rmse2 = sqrt(mean(se2)),
                                                        rmse3 = sqrt(mean(se3)))

rmse.adm2.insample <- res.adm2.mat %>% filter(!is.na(mean)) %>%group_by(admin2,true_mean) %>% reframe(rmse0 = sqrt(mean(se0)),
                                                                                                       rmse4 = sqrt(mean(se4)),
                                                                                                       rmse2 = sqrt(mean(se2)),
                                                                                                       rmse3 = sqrt(mean(se3)))
## plots --------------
pdf('sim8_complex_admin1_plots.pdf')
{
  
  ## ADMIN1
  # rmse
  y_range <- range(rmse.adm1$rmse0,#rmse.adm1$rmse4,
                   rmse.adm1$rmse2,rmse.adm1$rmse3)
  p0 <- rmse.adm1 %>% ggplot() + geom_point(aes(true_mean,rmse0)) + 
    xlab('') + ylab('') + ggtitle('Weighted estimates') +
    scale_y_continuous(limits = y_range) + 
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm1$rmse0),3)),
               hjust=-0.1, vjust = 1.1)
  p1 <- rmse.adm1 %>% ggplot() + geom_point(aes(true_mean,rmse4)) +
    xlab('') + ylab('') + ggtitle('Variance smoothing model') +
    scale_y_continuous(limits = y_range) +
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm1$rmse4),3)),
               hjust=-0.1, vjust = 1.1)
  p2 <- rmse.adm1 %>% ggplot() + geom_point(aes(true_mean,rmse2)) + 
    xlab('') + ylab('') + ggtitle('Correct variance') +
    scale_y_continuous(limits = y_range) + 
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm1$rmse2),3)),
               hjust=-0.1, vjust = 1.1)
  p3 <- rmse.adm1 %>% ggplot() + geom_point(aes(true_mean,rmse3)) + 
    xlab('') + ylab('') + ggtitle('Naive') +
    scale_y_continuous(limits = y_range) + 
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm1$rmse3),3)),
               hjust=-0.1, vjust = 1.1)
  p <- ggarrange(plotlist = list(p0,p1,
                                 p2,p3))
  print(annotate_figure(p, top = text_grob("RMSE of admin1 estimates (compared to truth)", face = "bold", size = 14)))
  
  
  # coverage
  y_range <- range(coverage.adm1$area_coverage0,#coverage.adm1$area_coverage4,
                   coverage.adm1$area_coverage2,coverage.adm1$area_coverage3)
  p0 <- coverage.adm1 %>% ggplot() + geom_point(aes(true_mean,area_coverage0)) + 
    xlab('') + ylab('') + ggtitle('Weighted estimates') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm1$area_coverage0),3)),
               hjust=-0.1, vjust = 1.1)
  p1 <- coverage.adm1 %>% ggplot() + geom_point(aes(true_mean,area_coverage4)) +
    xlab('') + ylab('') + ggtitle('Variance smoothing model') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm1$area_coverage4),3)),
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
  
  p <- ggarrange(plotlist = list(p0,p1,
                                 p2,p3))
  print(annotate_figure(p, top = text_grob("Coverage of 90% intervals for admin1 estimates", face = "bold", size = 14)))
  
  ## interval widths
  
  tmp <- res.adm1.mat %>% group_by(admin1) %>% reframe(weighted = mean(ciwidth0),
                                                       var_smooth = mean(ciwidth4),
                                                       correct_variance = mean(ciwidth2),
                                                       naive = mean(ciwidth3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) 
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,var_smooth)) + geom_abline(intercept = 0, slope=1)
  p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)
  
  p <- ggarrange(plotlist = list(p1,p2,
                                 p3),nrow=3)
  print(annotate_figure(p, top = text_grob("Average width of 90% intervals for admin1 estimates", face = "bold", size = 14)))
  
  ## ais
  
  tmp <- res.adm1.mat %>% group_by(admin1) %>% reframe(weighted = mean(i.score0),
                                                       var_smooth = mean(i.score4),
                                                       correct_variance = mean(i.score2),
                                                       naive = mean(i.score3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) 
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,var_smooth)) + geom_abline(intercept = 0, slope=1)
  p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)
  
  p <- ggarrange(plotlist = list(p1,p2,
                                 p3),nrow=3)
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

pdf('sim8_complex_admin2_plots.pdf')
{
  
  y_range <- range(rmse.adm2$rmse4,
                   rmse.adm2$rmse2,rmse.adm2$rmse3)
  p1 <- rmse.adm2 %>% ggplot() + geom_point(aes(true_mean,rmse4)) +
    xlab('') + ylab('') + ggtitle('FH with variance smoothing') +
    scale_y_continuous(limits = y_range) +
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm2$rmse4),3)),
               hjust=-0.1, vjust = 1.1)
  p2 <- rmse.adm2 %>% ggplot() + geom_point(aes(true_mean,rmse2)) + 
    xlab('') + ylab('') + ggtitle('Correct variance') +
    scale_y_continuous(limits = y_range) + 
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm2$rmse2),3)),
               hjust=-0.1, vjust = 1.1)
  p3 <- rmse.adm2 %>% ggplot() + geom_point(aes(true_mean,rmse3)) + 
    xlab('') + ylab('') + ggtitle('Naive FH') +
    scale_y_continuous(limits = y_range) + 
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm2$rmse3),3)),
               hjust=-0.1, vjust = 1.1)
  p <- ggarrange(plotlist = list(p2,p3,p1),nrow=3)
  print(annotate_figure(p, top = text_grob("RMSE of all admin2 estimates (compared to truth)", face = "bold", size = 14)))
  
  # coverage
  y_range <- range(coverage.adm2$area_coverage0,coverage.adm2$area_coverage4,
                   coverage.adm2$area_coverage2,coverage.adm2$area_coverage3,na.rm=T)
  p0 <- coverage.adm2 %>% ggplot() + geom_point(aes(true_mean,area_coverage0)) + 
    xlab('') + ylab('') + ggtitle('Weighted estimates') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2$area_coverage0,na.rm=T),3)),
               hjust=-0.1, vjust = 1.1)
  p1 <- coverage.adm2.insample %>% ggplot() + geom_point(aes(true_mean,area_coverage4)) +
    xlab('') + ylab('') + ggtitle('FH with variance smoothing') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2.insample$area_coverage4),3)),
               hjust=-0.1, vjust = 1.1)
  p2 <- coverage.adm2.insample %>% ggplot() + geom_point(aes(true_mean,area_coverage2)) + 
    xlab('') + ylab('') +ggtitle('Correct variance') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2.insample$area_coverage2),3)),
               hjust=-0.1, vjust = 1.1)
  p3 <- coverage.adm2.insample %>% ggplot() + geom_point(aes(true_mean,area_coverage3)) + 
    xlab('') + ylab('') +ggtitle('Naive FH') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2.insample$area_coverage3),3)),
               hjust=-0.1, vjust = 1.1)
  
  p <- ggarrange(plotlist = list(p0,
                                 p2,p3,p1))
  print(annotate_figure(p, top = text_grob("Coverage of 90% intervals for in-sample admin2 estimates", face = "bold", size = 14)))
  
  p1 <- coverage.adm2 %>% ggplot() + geom_point(aes(true_mean,area_coverage4)) +
    xlab('') + ylab('') + ggtitle('FH with variance smoothing') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2$area_coverage4),3)),
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
  
  p <- ggarrange(plotlist = list(p2,p3,p1),nrow=3)
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
                                                                                var_smooth = mean(ciwidth4),
                                                                                correct_variance = mean(ciwidth2),
                                                                                naive = mean(ciwidth3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) +
    ggtitle('in-sample estimates')
   
  tmp <- res.adm2.mat %>% group_by(admin2) %>% reframe(var_smooth= mean(ciwidth4),
                                                       correct_variance = mean(ciwidth2),
                                                       naive = mean(ciwidth3))
  
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,var_smooth)) + geom_abline(intercept = 0, slope=1)+
    ggtitle('all estimates')
  p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)+
    ggtitle('all estimates')
  
  p <- ggarrange(plotlist = list(p1,p2,p3),nrow=3)
  print(annotate_figure(p, top = text_grob("Average width of 90% intervals for admin2 estimates", face = "bold", size = 14)))
  
  # AIS
  tmp <- res.adm2.mat %>% filter(!is.na(mean)) %>%
    group_by(admin2,true_mean) %>% summarise(weighted = mean(i.score0),
                                             var_smooth = mean(i.score4),
                                             correct_variance = mean(i.score2),
                                             naive = mean(i.score3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) +
    ggtitle('in-sample estimates')
 
  tmp <- res.adm2.mat %>% 
    group_by(admin2,true_mean) %>% summarise(var_smooth = mean(i.score4),
                                             correct_variance = mean(i.score2),
                                             naive = mean(i.score3))
  
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,var_smooth)) + geom_abline(intercept = 0, slope=1) +
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

