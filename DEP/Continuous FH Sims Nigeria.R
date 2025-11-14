# Simulates continuous responses for the population (using height-for-age in Nigeria as a guide)
# takes two stage stratified random samples AND one stage strtatifed and then for each,
# fits FH models 1) under the correct sampling model 2) with correct variance but assuming normality 3) naively (with estimated variance and assuming normality)

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

# load Nigeria population -------

pop_dens <- rast("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Nigeria/Population/nga_u5_2018_100m.tif")

poly.adm2$admin_pop <- round(exact_extract(pop_dens, poly.adm2, 'sum'))

pop.admin2 <- merge(poly.adm2[,c('admin2.char','admin_pop')],admin.key)
pop.admin2 <- pop.admin2[order(pop.admin2$admin2),]

pop.admin1 <- pop.admin2 %>% group_by(admin1) %>% summarise(admin_pop = sum(admin_pop))

A_2to1 <- matrix(0,n_admin1,n_admin2)
for(area1 in 1:n_admin1){
  which.areas2 <- unique(admin.key[admin.key$admin1==area1,]$admin2)
  A_2to1[area1,which.areas2] <- pop.admin2[pop.admin2$admin2 %in% which.areas2,]$admin_pop/sum(pop.admin2[pop.admin2$admin2 %in% which.areas2,]$admin_pop)
}

# load covariates ----------------
load(file = "/Users/alanamcgovern/Desktop/Research/NGA_Covariates/Nigera_admin1_covariates.rda")
load(file = "/Users/alanamcgovern/Desktop/Research/NGA_Covariates/Nigera_admin2_covariates.rda")

if (st_is_longlat(poly.adm2) == FALSE) {
  poly.adm2 <- st_transform(poly.adm2, 4326)
}
centroids <- st_centroid(poly.adm2)
coords <- st_coordinates(centroids)

cmat_admin2$long <- coords[,('X')]
cmat_admin2$lat <- coords[,('Y')]

if (st_is_longlat(poly.adm1) == FALSE) {
  poly.adm2 <- st_transform(poly.adm1, 4326)
}
centroids <- st_centroid(poly.adm1)
coords <- st_coordinates(centroids)

cmat_admin1$long <- coords[,('X')]
cmat_admin1$lat <- coords[,('Y')]

# rescale some covariates
cmat_admin2 <- cmat_admin2 %>% mutate(elev = (elev - mean(elev))/sd(elev),
                                      precip = (precip - mean(precip))/sd(precip),
                                      temp = (temp - mean(temp))/sd(temp),
                                      long = (long-mean(long)/sd(long)),
                                      lat = (lat-mean(lat)/sd(lat)))

cmat_admin1 <- cmat_admin1 %>% mutate(elev = (elev - mean(elev))/sd(elev),
                                      precip = (precip - mean(precip))/sd(precip),
                                      temp = (temp - mean(temp))/sd(temp),
                                      long = (long-mean(long)/sd(long)),
                                      lat = (lat-mean(lat)/sd(lat)))

mean_covariates <- c('nt_lights_log','tthc_log','precip','temp','elev','density_log','long','lat')
var_covariates <- c('nt_lights_var','tthc_var','precip_var','temp_var','elev_var','area','pop_var') ## TAKE THE LOG OF THESE

# sort population into clusters -------------------

ea_totals <- data.table(
  NAME_1 = rep(c("Adamawa","Bauchi","Borno","Gombe","Taraba","Yobe",
            "Jigawa","Kaduna","Kano","Katsina","Kebbi","Sokoto","Zamfara",
            "Benue","Federal Capital Territory","Kogi","Kwara","Nasarawa","Niger","Plateau",
            "Abia","Anambra","Ebonyi","Enugu","Imo",
            "Akwa Ibom","Bayelsa","Cross River","Delta","Edo","Rivers",
            "Ekiti","Lagos","Ogun","Ondo","Osun","Oyo"),times=2),
  urban = c(rep(1,n_admin1),rep(0,n_admin1)),
  EAs = c(3040, 4539, 14760, 2035, 3799, 3851, ## urban
            2576, 9625, 17799, 8436, 3715, 1806, 3128,
            5179, 6682, 3412, 5345, 5469, 6466, 3456,
            5130, 7522, 3169, 8483, 5053,
            7717, 2931, 4738, 10214, 7259, 14120,
            4181, 32280, 21572, 6610, 10093, 16086,
          11743,19234,4229,7084,10686,7658, ## rural
            15685,18884,25380,20025,12225,13321,9773,
            15210,3901,7106,7986,7093,14305,13659,
            5243,5337,6416,8415,10520,
            8378,1638,3264,2930,1813,3500,
            1122,1121,4167,3847,2156,4452))

ea_totals <- merge(ea_totals, unique.array(admin.key[,c('admin1','admin1.char','NAME_1')]))

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
sim_pop_long <- merge(cluster_frame[,c('cluster','urban','strata')],sim_pop_long)

cluster_frame <- merge(cluster_frame,sim_pop_long[,.N,cluster],by='cluster')

# objects for sampling -------------
# if a cluster is selected, how many individuals will be selected?
cluster_frame$n <- rpois(nrow(cluster_frame),8) #calibrated from data

cluster_alloc <- data.table(
 NAME_1 = rep(c(
    "Adamawa","Bauchi","Borno","Gombe","Taraba","Yobe",
    "Jigawa","Kaduna","Kano","Katsina","Kebbi","Sokoto","Zamfara",
    "Benue","Federal Capital Territory","Kogi","Kwara","Nasarawa","Niger","Plateau",
    "Abia","Anambra","Ebonyi","Enugu","Imo",
    "Akwa Ibom","Bayelsa","Cross River","Delta","Edo","Rivers",
    "Ekiti","Lagos","Ogun","Ondo","Osun","Oyo"
  ),times=2),
urban = c(rep(1,n_admin1),rep(0,n_admin1)),
EAs = c(
    7,  7, 29,  7, 10, 11,
    5, 16, 19, 11,  8,  5,  8,
    12, 23, 13, 16, 17, 15,  9,
    19, 23, 13, 25, 13,
    20, 19, 21, 32, 30, 34,
    25, 49, 39, 24, 32, 35,
  29, 34,  9, 24, 25, 22,
    31, 27, 26, 28, 28, 32, 27,
    29, 13, 23, 21, 19, 25, 30,
    16, 15, 21, 12, 26,
    19, 10, 14,  8,  6,  8,
    6,  2,  7, 13,  6,  8)
)

cluster_alloc <- merge(cluster_alloc, unique.array(admin.key[,c('admin1','admin1.char','NAME_1')]))

# get population urban rural percent for weights
urban_df <- sim_pop_long[,mean(urban),by=NAME_1]
pop.admin1 <- merge(pop.admin1,merge(urban_df, unique.array(admin.key[,c('admin1','admin1.char','NAME_1')])))
cluster_alloc <- cluster_alloc %>% mutate(strata_total =  pop.admin1$admin_pop[admin1]*(urban*pop.admin1$V1[admin1] + (1-urban)*(1-pop.admin1$V1[admin1])))

cluster_frame <- merge(cluster_frame,cluster_alloc,by=c('admin1','NAME_1','urban','admin1.char')) %>% mutate(wt_complex = strata_total/(EAs*n))
cluster_frame[,EAs:=NULL]
cluster_frame[,strata_total:=NULL]

cluster_alloc$strata <- cluster_alloc$admin1 + cluster_alloc$urban*n_admin1

# srs_sample_sizes <- sim_pop %>% group_by(admin1) %>% summarise(N=sum(admin_pop))
# sample_prop <- 0.00035
# # assume sample size is roughly proportional to population size (add some noise)
# srs_sample_sizes$n <- round(rnorm(nrow(srs_sample_sizes),sample_prop*srs_sample_sizes$N,10))
# srs_sample_sizes[srs_sample_sizes$n < 2,]$n <- 2
# srs_sample_sizes$wt_srs <- srs_sample_sizes$N/srs_sample_sizes$n

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

# simulate data at population level (just using admin1 means for now) -----

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
 load(file = 'Nigeria_haz_admin1_weighted_estimates.rda')
 load(file = 'Nigeria_haz_admin2_weighted_estimates.rda')
#load(file = 'Nigeria_whz_admin1_weighted_estimates.rda')
#load(file = 'Nigeria_whz_admin2_weighted_estimates.rda')
 
admin1_means <- admin1.dir[order(admin1.dir$admin1),]$mean

within_admin1_var <- 2 # from data
rho <- 1 # proportion which is individual level (as opposed to cluster level) -- in data it is about 95%

sd_cluster <- sqrt((1-rho)*within_admin1_var)
cluster_frame$eps <- rnorm(nrow(cluster_frame),0,sd_cluster)

sd_e <- sqrt(rho*within_admin1_var)

sim_pop_long <- sim_pop_long[cluster_frame[,c('cluster','eps')], on = 'cluster']
sim_pop_long[, value := admin1_means[admin1] + eps + rnorm(.N,0,sd_e)]

# sample using 2 stage sampling and get direct estimates -----------
nsim <- 500

direct.adm2.complex <- direct.adm1.complex <- sampled_clusters <- list()
for(k in 1:nsim){
  cat(k,'\n')
  
  # which clusters will we sample?
  cluster_sample_ids <- unlist(sapply(1:(2*n_admin1),function(i){
    sample(cluster_frame[cluster_frame$strata==i,]$cluster,
           cluster_alloc[cluster_alloc$strata==i,]$EAs,
           prob = cluster_frame[cluster_frame$strata==i,]$N)
  }))
  cluster_sample <- cluster_frame[cluster_frame$cluster %in% cluster_sample_ids,c('cluster','n','wt_complex')]
  
  sim_sample <- sim_pop_long[
    cluster_sample,                   # right table
    on = "cluster"                   # join by cluster
  ][, .SD[sample(.N, min(.N, unique(n)))], by = cluster] # for each cluster, sample n
  
  sim_sample[,n:=NULL]
  
  ## save info (weights and sample sizes) on sampled cluster for Satt approx
  sampled_clusters[[k]] <- sim_sample %>% group_by(admin1,cluster) %>% reframe(n=n(),wt=unique(wt_complex))
  
  
  # get direct estimates (not accounting for stratification) --------
  options(survey.lonely.psu = "adjust")
  sim_sample$admin1.char <- factor(sim_sample$admin1.char)
  my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"), 
                                    strata = ~admin1.char, nest = T, weights = ~wt_complex, data = sim_sample)
  
  # ## admin2
  # x <- mapply(admin2.HT.withNA, which.area = 1:n_admin2)
  # admin2.dir <- data.frame(t(x))
  # colnames(admin2.dir) <- c('admin2','mean','variance')
  # sample_info <- merge(sim_sample[,length(unique(cluster)),by=admin2],sim_sample[,.N,by=admin2])
  # colnames(sample_info)[2:3] <- c('n_clusters','n_obs')
  # 
  # direct.adm2.complex[[k]] <- merge(sample_info,admin2.dir,all=T)
  
  ## admin1
  x <- mapply(admin1.HT.withNA, which.area = 1:n_admin1)
  admin1.dir <- data.frame(t(x))
  colnames(admin1.dir) <- c('admin1','mean','variance')
  sample_info <- merge(sim_sample[,length(unique(cluster)),by=admin1],sim_sample[,.N,by=admin1])
  colnames(sample_info)[2:3] <- c('n_clusters','n_obs')

  direct.adm1.complex[[k]] <- merge(sample_info,admin1.dir,all=T)
  
}

#adm2.dir.complex.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct.adm2.complex, function(x){x$mean[i]})),na.rm = T)})
adm1.dir.complex.variance <- sapply(1:n_admin1, function(i){var(unlist(lapply(direct.adm1.complex, function(x){x$mean[i]})),na.rm = T)})

# diagnostics ----

## have variances converged?
par(mfrow=c(3,3))
for(area in 1:n_admin1){
  plot(sapply(1:nsim, function(i){var(unlist(lapply(direct.adm1.complex[1:i], function(x){x$mean[area]})),na.rm = T)}),type='l',main = paste('Area ', area))
}

for(area in 1:n_admin1){
  hist(unlist(lapply(direct.adm1.complex, function(x){x$var[area]})))
  abline(v=adm1.dir.complex.variance[area],col='red')
}


par(mfrow=c(3,3))
for(area in sample(1:n_admin2,9)){
  plot(sapply(1:nsim, function(i){var(unlist(lapply(direct.adm2.complex[1:i], function(x){x$mean[area]})),na.rm = T)}),type='l',main = paste('Area ', area))
}

# get FH estimates from direct --------
# CHOOSE SRS OR COMPLEX SAMPLING
#adm2.dir.variance <- adm2.dir.complex.variance
adm1.dir.variance <- adm1.dir.complex.variance
#direct.adm2.estimates <- direct.adm2.complex
direct.adm1.estimates <- direct.adm1.complex

results.adm1 <- results.adm2 <- list() # record estimates (direct and fay herriot)
#covariates = mean_covariates
covariates = NULL
random_model = c("iid","bym2")[1]
for(k in 1:100){
  cat(k,'\n')
  
  # sim.admin2.dir <- direct.adm2.estimates[[k]]
  # sim.admin2.dir$var.truth <- adm2.dir.variance
  # sim.admin2.dir$mean.sim <- sapply(1:n_admin2,function(j){
  #   rnorm(1,admin2_means[j],sqrt(adm2.dir.variance[j]))
  # })
  # sim.admin2.dir <- sim.admin2.dir %>% mutate(variance = ifelse(variance < 1e-5,NA,variance))
  # sim.admin2.dir <- sim.admin2.dir %>% mutate(mean = ifelse(is.na(variance),NA,mean),
  #                                             var.truth = ifelse(is.na(variance),NA,var.truth),
  #                                             mean.sim = ifelse(is.na(variance),NA,mean.sim))
  # sim.admin2.dir <- merge(sim.admin2.dir,cmat_admin2,by='admin2')
  
  sim.admin1.dir <- direct.adm1.estimates[[k]]
  sim.admin1.dir$var.truth <- adm1.dir.variance
  # sim.admin1.dir$mean.sim <- sapply(1:n_admin1,function(j){
  #   rnorm(1,admin1_means[j],sqrt(adm1.dir.variance[j]))
  # })
  # sim.admin1.dir <- merge(sim.admin1.dir,cmat_admin1,by='admin1')
  # 
  # # model is specified exactly correctly ----
  
  # eta.samples1a <- fit_inla_with_samples(outcome = 'mean.sim',
  #                                        scale_var = 'var.truth',
  #                                        data = sim.admin1.dir,
  #                                        random_effect = "admin1",
  #                                        covariates = covariates,
  #                                        adj_matrix = admin1.mat,
  #                                        random_model = random_model)
  # 
  # eta.samples1 <- fit_inla_with_samples(outcome = 'mean.sim',
  #                                       scale_var = 'var.truth',
  #                                       data = sim.admin2.dir,
  #                                       random_effect = "admin2",
  #                                       covariates = covariates,
  #                                       adj_matrix = admin2.mat,
  #                                       random_model = random_model)
  
  ## model is specified with correct variance ----
  eta.samples2a <- fit_inla_with_samples(outcome = 'mean',
                                         scale_var = 'var.truth',
                                         data = sim.admin1.dir,
                                         random_effect = "admin1",
                                         covariates = covariates,
                                         adj_matrix = admin1.mat,
                                         random_model = random_model)
  
  # eta.samples2 <- fit_inla_with_samples(outcome = 'mean',
  #                                       scale_var = 'var.truth',
  #                                       data = sim.admin2.dir,
  #                                       random_effect = "admin2",
  #                                       covariates = covariates,
  #                                       adj_matrix = admin2.mat,
  #                                       random_model = random_model)
  # 
  ## naive FH (estimated variance, assuming normality) ----
  eta.samples3a <- fit_inla_with_samples(outcome = 'mean',
                                         scale_var = 'variance',
                                         data = sim.admin1.dir,
                                         random_effect = "admin1",
                                         covariates = covariates,
                                         adj_matrix = admin1.mat,
                                         random_model = random_model)
  
  # eta.samples3 <- fit_inla_with_samples(outcome = 'mean',
  #                                       scale_var = 'variance',
  #                                       data = sim.admin2.dir,
  #                                       random_effect = "admin2",
  #                                       covariates = covariates,
  #                                       adj_matrix = admin2.mat,
  #                                       random_model = random_model)
  
  ## naive FH with variance smoothing
  Cons <- v_hat_scaled <- df <- rep(NA,n_admin1)
  for(area in 1:n_admin1){
    tmp <- sampled_clusters[[k]] %>% filter(admin1==area)
    N = nrow(tmp)
    
    omega <- tmp$n*tmp$wt ## Nx1
    D = diag(omega^2)
    M = diag(1,N) - (rep(1,N) %*% t(omega))/as.numeric(t(rep(1,N)) %*% t(t(omega)))
    S = diag(1/tmp$n) # structure of the covariance matrix
    q <- eigen(S^(1/2)%*%t(M)%*%D%*%M%*%S^(1/2))$values
    
    # outputs
    Cons[area] <- c(sum(tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2) #multiplier for variance term in weighted mean distribution
    
    # recale to be in terms of chi square approximation
    v_hat_scaled[area] <- sum(q)/(sum(q^2)*N/(sum(tmp$wt*tmp$n)^2*(N-1)))*admin1.dir$variance[area]
    
    df[area] <- sum(q)^2/sum(q^2)
  }
  
  data_list = list(m=n_admin1,
                   m_data=n_admin1,
                   data_areas = 1:n_admin1,
                   y=sim.admin1.dir$mean,
                   v_hat_scaled = v_hat_scaled,
                   # df for chi square approx
                   df = df,
                   # constant for likelihood
                   Cons = Cons,
                   # covariates in mean model
                   p_mean = 1,
                   X = matrix(1,n_admin1,1),
                   # p_mean = length(mean_covariates) +1,
                   # X = cbind(1,admin1.dir[,mean_covariates]),
                   # covariates in log-variance model (none for now)
                   p_var = 1,
                   Z = matrix(1,n_admin1,1)) 
  
  fit <- satt_mod$sample(
    data = data_list,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh=0
  )
  
  eta.samples4a <- as.data.frame(fit$draws(format='df',variables=c('theta'))) %>% dplyr::select(-c('.chain','.draw','.iteration'))
  
  
  # 
  ## record estimates -----
  
  # results.adm2[[k]] <- merge(sim.admin2.dir[,c('admin2','mean','mean.sim','variance','var.truth')],
  #                            data.frame(admin2 = 1:n_admin2,
  #                                       # correctly specified
  #                                       mean.fh1  = apply(eta.samples1,2,mean),
  #                                       var.fh1 = apply(eta.samples1,2,var),
  #                                       lower.fh1 = apply(eta.samples1,2,quantile,prob=0.05),
  #                                       upper.fh1 = apply(eta.samples1,2,quantile,prob=0.95),
  #                                       # uses correct variance, but normality is assumed
  #                                       mean.fh2  = apply(eta.samples2,2,mean),
  #                                       var.fh2 = apply(eta.samples2,2,var),
  #                                       lower.fh2 = apply(eta.samples2,2,quantile,prob=0.05),
  #                                       upper.fh2 = apply(eta.samples2,2,quantile,prob=0.95),
  #                                       # assumes and normality and uses estimated variance
  #                                       mean.fh3  = apply(eta.samples3,2,mean),
  #                                       var.fh3 = apply(eta.samples3,2,var),
  #                                       lower.fh3 = apply(eta.samples3,2,quantile,prob=0.05),
  #                                       upper.fh3 = apply(eta.samples3,2,quantile,prob=0.95),
  #                                       sim = k,
  #                                       true_mean = admin2_means),all=T)
  
  results.adm1[[k]] <- merge(sim.admin1.dir[,c('admin1','mean',#'mean.sim',
                                               'variance','var.truth')],
                             data.frame(admin1 = 1:n_admin1,
                                        # correctly specified
                                        # mean.fh1  = apply(eta.samples1a,2,mean),
                                        # var.fh1 = apply(eta.samples1a,2,var),
                                        # lower.fh1 = apply(eta.samples1a,2,quantile,prob=0.05),
                                        # upper.fh1 = apply(eta.samples1a,2,quantile,prob=0.95),
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
                                        # variance smoothing with Satterwhaite
                                        mean.fh4  = apply(eta.samples4a,2,mean),
                                        var.fh4 = apply(eta.samples4a,2,var),
                                        lower.fh4 = apply(eta.samples4a,2,quantile,prob=0.05),
                                        upper.fh4 = apply(eta.samples4a,2,quantile,prob=0.95),
                                        sim = k,
                                        true_mean = admin1_means),all=T)
  
}

# SAVE/LOAD RESULTS ---------

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
# 
 save(sim_pop_long,file='sim7_popdata.rda')
  save(admin1_means,file='sim7_admin1_means.rda')
  save(cluster_frame,file='sim7_clusterframe.rda')
  save(cluster_alloc,file='sim7_clusteralloc.rda')
save(sampled_clusters,file='sim7_clusteralloc.rda')
# save(direct.adm2.complex,file='sim3-6_complex_admin2.direct.rda')
save(direct.adm1.complex,file='sim7_complex_admin1.direct.rda')
#save(results.adm2,file = 'sim6_complex_admin2.fh.rda')
save(results.adm1,file = 'sim7_complex_admin1.fh.rda')

# load(file='sim3_admin2_means.rda')
#load('sim3_complex_admin1.fh.rda')
##load('sim3_complex_admin2.fh.rda')
# load(file = 'sim6_complex_admin1.fh.rda')
# load(file = 'sim6_complex_admin2.fh.rda')
# 
# load(file='sim3-6_complex_admin1.direct.rda')

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

# res.adm2.mat <- do.call(rbind,results.adm2)
# res.adm2.mat <- res.adm2.mat %>% mutate(lower.direct = mean - 1.645*sqrt(variance),
#                                         upper.direct = mean + 1.645*sqrt(variance),
#                                         ciwidth0 = upper.direct - lower.direct,
#                                         ciwidth1 = upper.fh1 - lower.fh1,
#                                         ciwidth2 = upper.fh2 - lower.fh2,
#                                         ciwidth3 = upper.fh3 - lower.fh3,
#                                         se0 = (true_mean - mean)^2,
#                                         se1 = (true_mean - mean.fh1)^2,
#                                         se2 = (true_mean - mean.fh2)^2,
#                                         se3 = (true_mean - mean.fh3)^2,
#                                         cov0 = lower.direct<true_mean & upper.direct>true_mean,
#                                         cov1 = lower.fh1<true_mean & upper.fh1>true_mean,
#                                         cov2 = lower.fh2<true_mean & upper.fh2>true_mean,
#                                         cov3 = lower.fh3<true_mean & upper.fh3>true_mean,
#                                         i.score0 = ciwidth0 + 20*(lower.direct-true_mean)*(true_mean<lower.direct) + 20*(true_mean-upper.direct)*(true_mean>upper.direct),
#                                         i.score1 = ciwidth1 + 20*(lower.fh1-true_mean)*(true_mean<lower.fh1) + 20*(true_mean-upper.fh1)*(true_mean>upper.fh1),
#                                         i.score2 = ciwidth2 + 20*(lower.fh2-true_mean)*(true_mean<lower.fh2) + 20*(true_mean-upper.fh2)*(true_mean>upper.fh2),
#                                         i.score3 = ciwidth3 + 20*(lower.fh3-true_mean)*(true_mean<lower.fh3) + 20*(true_mean-upper.fh3)*(true_mean>upper.fh3))

coverage.adm1 <- res.adm1.mat %>% group_by(admin1,true_mean) %>% reframe(area_coverage0 = mean(cov0,na.rm=T),
                                                                         area_coverage4 = mean(cov4),
                                                                         area_coverage2 = mean(cov2),
                                                                         area_coverage3 = mean(cov3))


# coverage.adm2 <- res.adm2.mat %>% group_by(admin2,true_mean) %>% reframe(area_coverage0 = mean(cov0,na.rm=T),
#                                                                          area_coverage1 = mean(cov1),
#                                                                          area_coverage2 = mean(cov2),
#                                                                          area_coverage3 = mean(cov3))
# 
# coverage.adm2.insample <- res.adm2.mat %>% filter(!is.na(mean)) %>% 
#   group_by(admin2,true_mean) %>% reframe(area_coverage0 = mean(cov0),
#                                          area_coverage1 = mean(cov1),
#                                          area_coverage2 = mean(cov2),
#                                          area_coverage3 = mean(cov3))

rmse.adm1 <- res.adm1.mat %>% group_by(admin1,true_mean) %>% reframe(rmse0 = sqrt(mean(se0)),
                                                                     rmse4 = sqrt(mean(se4)),
                                                                    rmse2 = sqrt(mean(se2)),
                                                                    rmse3 = sqrt(mean(se3)))


# rmse.adm2 <- res.adm2.mat %>% group_by(admin2,true_mean) %>% reframe(rmse1 = sqrt(mean(se1)),
#                                                         rmse2 = sqrt(mean(se2)),
#                                                         rmse3 = sqrt(mean(se3)))
# 
# rmse.adm2.insample <- res.adm2.mat %>% filter(!is.na(mean)) %>%group_by(admin2,true_mean) %>% reframe(rmse0 = sqrt(mean(se0)),
#                                                                                                        rmse1 = sqrt(mean(se1)),
#                                                                                                        rmse2 = sqrt(mean(se2)),
#                                                                                                        rmse3 = sqrt(mean(se3)))


## PLOTS ------

pdf('sim7_complex_admin1_plots.pdf')
{
  
  ## ADMIN1
  # rmse
  y_range <- range(rmse.adm1$rmse0,rmse.adm1$rmse4,rmse.adm1$rmse2,rmse.adm1$rmse3)
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
  p <- ggarrange(plotlist = list(p0,p1,p2,p3))
  print(annotate_figure(p, top = text_grob("RMSE of admin1 estimates (compared to truth)", face = "bold", size = 14)))
  

  # coverage
  y_range <- range(coverage.adm1$area_coverage0,coverage.adm1$area_coverage4,coverage.adm1$area_coverage2,coverage.adm1$area_coverage3)
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
  
  p <- ggarrange(plotlist = list(p0,p1,p2,p3))
  print(annotate_figure(p, top = text_grob("Coverage of 90% intervals for admin1 estimates", face = "bold", size = 14)))
  
  ## interval widths
  
  tmp <- res.adm1.mat %>% group_by(admin1) %>% reframe(weighted = mean(ciwidth0),
                                                       var_smooth = mean(ciwidth4),
                                                       correct_variance = mean(ciwidth2),
                                                       naive = mean(ciwidth3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) 
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,var_smooth)) + geom_abline(intercept = 0, slope=1)
  p3 <- tmp %>% ggplot() + geom_point(aes(naive,correct_variance)) + geom_abline(intercept = 0, slope=1)
  
  p <- ggarrange(plotlist = list(p1,p2,p3),nrow=3)
  print(annotate_figure(p, top = text_grob("Average width of 90% intervals for admin1 estimates", face = "bold", size = 14)))
  
  ## ais
  
  tmp <- res.adm1.mat %>% group_by(admin1) %>% reframe(weighted = mean(i.score0),
                                                       var_smooth = mean(i.score4),
                                                       correct_variance = mean(i.score2),
                                                       naive = mean(i.score3))
  
  p1 <- tmp %>% ggplot() + geom_point(aes(naive,weighted)) + geom_abline(intercept = 0, slope=1) 
  p2 <- tmp %>% ggplot() + geom_point(aes(naive,var_smooth)) + geom_abline(intercept = 0, slope=1)
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

#pdf('sim6_complex_admin2_plots.pdf')
{

  y_range <- range(rmse.adm2$rmse1,rmse.adm2$rmse2,rmse.adm2$rmse3)
  p1 <- rmse.adm2 %>% ggplot() + geom_point(aes(true_mean,rmse1)) + 
    xlab('') + ylab('') + ggtitle('Correct sampling model') +
    scale_y_continuous(limits = y_range) + 
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm2$rmse1),3)),
               hjust=-0.1, vjust = 1.1)
  p2 <- rmse.adm2 %>% ggplot() + geom_point(aes(true_mean,rmse2)) + 
    xlab('') + ylab('') + ggtitle('Correct variance') +
    scale_y_continuous(limits = y_range) + 
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm2$rmse2),3)),
               hjust=-0.1, vjust = 1.1)
  p3 <- rmse.adm2 %>% ggplot() + geom_point(aes(true_mean,rmse3)) + 
    xlab('') + ylab('') + ggtitle('Naive') +
    scale_y_continuous(limits = y_range) + 
    geom_label(y=Inf,x=-Inf,label=paste('RMSE = ',round(mean(rmse.adm2$rmse3),3)),
               hjust=-0.1, vjust = 1.1)
  p <- ggarrange(plotlist = list(p1,p2,p3))
  print(annotate_figure(p, top = text_grob("RMSE of all admin2 estimates (compared to truth)", face = "bold", size = 14)))
  
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
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2.insample$area_coverage1),3)),
               hjust=-0.1, vjust = 1.1)
  p2 <- coverage.adm2.insample %>% ggplot() + geom_point(aes(true_mean,area_coverage2)) + 
    xlab('') + ylab('') +ggtitle('Correct variance') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2.insample$area_coverage2),3)),
               hjust=-0.1, vjust = 1.1)
  p3 <- coverage.adm2.insample %>% ggplot() + geom_point(aes(true_mean,area_coverage3)) + 
    xlab('') + ylab('') +ggtitle('Naive') +
    scale_y_continuous(limits = y_range) + geom_hline(yintercept = 0.9,col='red') +
    geom_label(y=Inf,x=-Inf,label=paste('Coverage = ',round(mean(coverage.adm2.insample$area_coverage3),3)),
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
#dev.off()

## testing ---------------

## is sampling variance distributed chi-square with df = n_clusters - 1 ?
## didn't save that info for admin2, but for admin1 looks good

test_mat <- do.call(rbind,direct.adm1.complex)
test_mat$var_truth <- rep(adm1.dir.complex.variance,1000)

test_mat <- merge(test_mat,unique.array(admin.key[,c('admin1','NAME_1')]))
test_mat <- merge(test_mat,cluster_alloc,by='NAME_1')


test_mat$test_stat <- (test_mat$EAs-1)*test_mat$variance/test_mat$var_truth


hist(test_mat[NAME_1=='Kaduna',]$test_stat,prob=T)
lines(1:200,dchisq(1:200,42))

res.adm2.mat %>% group_by(admin2,true_mean) %>% reframe(val0 = sqrt(mean(se0,,na.rm=T)),
                                                        val1 = sqrt(mean(se1)),
                                                        val2 = sqrt(mean(se2)),
                                                        val3 = sqrt(mean(se3)))



