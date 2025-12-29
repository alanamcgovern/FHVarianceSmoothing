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
library(tidyverse)
library(sampling)

source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")
setwd('/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing')

# load geometry (country specific) ------
setwd('Kenya_Example')
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

admin1.nbs <- poly2nb(poly.adm1)
nodes1 <- NULL
for(i in 1:n_admin1){
  nodes1 <- rbind(nodes1,data.frame(node1=i,node2=as.numeric(admin1.nbs[[i]])))
}
nodes1 <- nodes1[nodes1$node1<nodes1$node2,]

Q.admin1 <- -admin1.mat
Q.admin1 <- sapply(1:nrow(Q.admin1),function(i){sum(I(Q.admin1[i,]!=0))})*Q.admin1
diag(Q.admin1) <- sapply(1:nrow(Q.admin1),function(i){sum(I(Q.admin1[i,]!=0))})
diag(Q.admin1)[diag(Q.admin1)==0] <- 1
Q1_inv <- INLA:::inla.ginv(as.matrix(Q.admin1))
Q1_scaled <- INLA::inla.scale.model(Q.admin1, constr=list(A=matrix(1,nrow=1,ncol=n_admin1), e=0))
Q1_scaled_inv <- INLA:::inla.ginv(as.matrix(Q1_scaled))


admin2.mat <- nb2mat(poly2nb(poly.adm2), zero.policy = TRUE)
colnames(admin2.mat) <- rownames(admin2.mat) <- poly.adm2$admin2.char

admin2.nbs <- poly2nb(poly.adm2)
nodes2 <- NULL
for(i in 1:n_admin2){
  nodes2 <- rbind(nodes2,data.frame(node1=i,node2=as.numeric(admin2.nbs[[i]])))
}
nodes2 <- nodes2[nodes2$node1<nodes2$node2,]

Q.admin2 <- -admin2.mat
Q.admin2 <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})*Q.admin2
diag(Q.admin2) <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})
diag(Q.admin2)[diag(Q.admin2)==0] <- 1
Q2_inv <- INLA:::inla.ginv(as.matrix(Q.admin2))
Q2_scaled <- INLA::inla.scale.model(Q.admin2, constr=list(A=matrix(1,nrow=1,ncol=n_admin2), e=0))
Q2_scaled_inv <- INLA:::inla.ginv(as.matrix(Q2_scaled))


# load population (country specific) -------

pop_dens <- rast("ken_u5_2022_1km.tif")

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
cluster_alloc$strata <- cluster_alloc$admin1 + cluster_alloc$urban*n_admin1

## oversample urban even more
   cluster_alloc[cluster_alloc$urban==1,]$EAs <- 2*cluster_alloc[cluster_alloc$urban==1,]$EAs
   cluster_alloc[cluster_alloc$urban==0,]$EAs <- round(0.5*cluster_alloc[cluster_alloc$urban==0,]$EAs)


## urban is oversampled in most (90%) areas
# merge(ea_totals %>% group_by(NAME_1) %>% summarise(urb_frac_pop = sum(urban*EAs)/sum(EAs)),cluster_alloc %>% 
#         group_by(NAME_1) %>% summarise(urb_frac_sample = sum(urban*EAs)/sum(EAs))) %>% 
#   ggplot() + geom_point(aes(urb_frac_pop,urb_frac_sample)) + geom_abline(intercept = 0,slope = 1)

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
  
  n <- .N                                  # number of EAs in this admin1
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
sim_pop_long <- merge(cluster_frame,sim_pop_long,by=c('NAME_1','admin1','admin1.char','admin2','cluster'))

# add population totals to cluster frame
cluster_frame <- merge(cluster_frame,sim_pop_long[,.N,cluster],by='cluster')

# objects for sampling ---------
# average number of individuals to sample from each cluster 
n_per_EA <- 10

cluster_frame$n <-  rnbinom(nrow(cluster_frame),size=10,mu = n_per_EA) + 1 #number that will be sampled from that cluster if sampled
#cluster_frame$n <-  rnbinom(nrow(cluster_frame), size=10, mu = n_per_EA)
#cluster_frame$n <-  rpois(nrow(cluster_frame),n_per_EA)
cluster_frame[cluster_frame$n<1,]$n <- 1
cluster_frame$n <- pmin(cluster_frame$n,0.75*round(cluster_frame$N))

# get population urban rural percent for weights
cluster_alloc <- left_join(cluster_alloc, sim_pop_long[,.N,by=strata],by='strata') %>% rename(strata_total=N)

cluster_frame <- merge(cluster_frame,cluster_alloc,by=c('admin1','NAME_1','urban','admin1.char','strata')) %>% 
  mutate(wt = strata_total/(EAs*n)/1e3, # individual sampling weight for all in cluster
         pik = EAs*N/strata_total) # inclusion probability of cluster
cluster_frame[,EAs:=NULL]
cluster_frame[,strata_total:=NULL]

# helper function ----------

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

# generate/load covariates ---------
# use same simulated covariates across all settings

# cmat_admin2 <- data.frame(admin2 = 1:n_admin2,
#                           X1 =  rnorm(n_admin2,0,1),
#                           X2 =  rnorm(n_admin2,0,1),
#                           X3 =  rnorm(n_admin2,0,1),
#                           X4 =  rnorm(n_admin2,0,1),
#                           X5 =  rnorm(n_admin2,0,1),
#                           X6 =  rnorm(n_admin2,0,1))
# 
# save(cmat_admin2,file='/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Simulations/simulated_covariates.rda')

load(file='/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Simulations/simulated_covariates.rda')
cmat_admin2$urb_frac <- sim_pop_long[order(admin2),mean(urban),by=admin2]$V1

mean_covariates <- c('X1','X2','X3')
var_covariates <- c('X4','X5','X6')

mean_model_matrix <- cbind(rep(1,n_admin2),cmat_admin2[,mean_covariates])
var_model_matrix <- cbind(rep(1,n_admin2),cmat_admin2[,var_covariates])

# simulation parameters and generate surface ------
sig_u <- c(sqrt(0.1),sqrt(0.05))
phi <- c(0.75,0.75)

W1 <- sig_u[1]^2*((1-phi[1])*diag(1,n_admin2) + phi[1]*Q2_scaled_inv)
W2 <- sig_u[2]^2*((1-phi[2])*diag(1,n_admin2) + phi[2]*Q2_scaled_inv)

beta <- c(-1, -0.1, -0.1, 0.2)
gamma <- c(0.5, -0.15, -0.1, 0.25)

# delta_mean <- 0
# kappa_var <- 1

# setting 1/2/4
delta_mean <- 2
kappa_var <- 1

# setting 3
# delta_mean <- rnorm(n_admin2,1,1)
# kappa_var <- rnorm(n_admin2,2,0.5)

# admin2_rural_mean <- beta %*% t(mean_model_matrix) + as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),W1))
# admin2_strata_mean <- c(admin2_rural_mean, admin2_rural_mean + delta_mean)
# admin2_mean <- admin2_strata_mean[1:n_admin2]*(1-cmat_admin2$urb_frac) + admin2_strata_mean[n_admin2 + 1:n_admin2]*cmat_admin2$urb_frac
# 
# sig2_R <- exp(gamma %*% t(var_model_matrix) + as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),W2)))
# sig2 <- c(sig2_R, kappa_var*sig2_R)
# 
# # setting 1-3
# sim_pop_long[, value := admin2_strata_mean[admin2_strata] + rnorm(.N,0,sqrt(sig2[admin2_strata]))]

# setting 4
 df_t <- 5
 s = sqrt(sig2*(df_t-2)/df_t)
# sim_pop_long[, value := admin2_strata_mean[admin2_strata] + s[admin2_strata]*rt(.N,df_t)]

## check
# plot(admin2_mean,sim_pop_long[order(admin2),mean(value),by=admin2]$V1)
# abline(0,1)

# take sample and get direct estimates --------
direct <- sampled_clusters <- true_vals <- list()
for(k in 1:100){
  
  cat(k,'\n')
  
  admin2_rural_mean <- beta %*% t(mean_model_matrix) + as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),W1))
  admin2_strata_mean <- c(admin2_rural_mean, admin2_rural_mean + delta_mean)
  admin2_mean <- admin2_strata_mean[1:n_admin2]*(1-cmat_admin2$urb_frac) + admin2_strata_mean[n_admin2 + 1:n_admin2]*cmat_admin2$urb_frac

  sig2_R <- exp(gamma %*% t(var_model_matrix) + as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),W2)))
  sig2 <- c(sig2_R, kappa_var*sig2_R)

 # sim_pop_long[, value := admin2_strata_mean[admin2_strata] + rnorm(.N,0,sqrt(sig2[admin2_strata]))]
  
   sim_pop_long[, value := admin2_strata_mean[admin2_strata] + s[admin2_strata]*rt(.N,df_t)]

  # PPS sampling of clusters
  cluster_sample_ids <- unlist(sapply(unique(cluster_frame$strata),function(i){
    cluster_frame[cluster_frame$strata==i,]$cluster[UPsystematic(cluster_frame[cluster_frame$strata==i,]$pik)>=1]
  }))
  cluster_sample <- cluster_frame[cluster_frame$cluster %in% cluster_sample_ids,c('cluster','n','wt')]
  
  sim_sample <- sim_pop_long[
    cluster_sample,                   # right table
    on = "cluster"                   # join by cluster
  ][, .SD[sample(.N, unique(n))], by = cluster] # for each cluster, sample n
  
  sim_sample[,n:=NULL]
  
  ## save info (weights and sample sizes) on sampled cluster for Satt approx
  sampled_clusters[[k]] <- sim_sample %>% group_by(admin1,admin2,urban,strata,cluster) %>% reframe(n=n(),wt=unique(wt))

  true_vals[[k]] <- list(sig2 = sig2,
                         admin2_mean = admin2_mean,
                         admin2_strata_mean = admin2_strata_mean)

  my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"),strata = ~strata,
                                    weights = ~wt, data = sim_sample)
  
  x <- mapply(admin2.HT.withNA, which.area = 1:n_admin2)
  admin.dir <- data.frame(t(x))
  colnames(admin.dir) <- c('admin2','mean','variance')
  direct[[k]] <- admin.dir
}


# save all the relevant objects -------
setting <- 14

setwd(paste0('/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Simulations'))

folder <- paste0('Sim',setting)
if(!dir.exists(folder)){
  dir.create(folder, recursive = TRUE)
}
setwd(folder)

sim_pop_long[,value:=NULL]
save(sim_pop_long,file='popdata.rda')
save(cmat_admin2,file='cmat_admin2.rda')
save(cluster_alloc,file='cluster_alloc.rda')
save(cluster_frame,file='cluster_frame.rda')
save(sampled_clusters,file='sampled_clusters.rda')
save(direct,file='direct.rda')
save(true_vals,file='true_vals.rda')

objects <- list(admin.key = admin.key,
                nodes1 = nodes1,
                nodes2 = nodes2,
                car_scale1 = Q1_scaled[1,1]/Q.admin1[1,1],
                car_scale2 = Q2_scaled[1,1]/Q.admin2[1,1],
                A_2to1 = A_2to1)

save(objects,file='objects.rda')

params <- list(#sig2 = sig2,
               beta = beta,
               gamma = gamma,
               #admin2_mean = admin2_mean,
               #admin2_strata_mean = admin2_strata_mean,
               sig_u = sig_u,
               phi = phi)

save(params,file='params.rda')


