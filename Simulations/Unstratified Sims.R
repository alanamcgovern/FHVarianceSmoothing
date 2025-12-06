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


# sampling info NOT STRATIFIED (survey specific) -------------

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

# NOT STRATIFIED SAMPLING!!!
ea_totals <- ea_totals %>% pivot_wider(id_cols = NAME_1, names_from = urban, values_from = EAs) %>% mutate(EAs = `1`+`0`)
ea_totals <- merge(ea_totals[,c('NAME_1','EAs')], unique.array(admin.key[,c('admin1','admin1.char','NAME_1')]))
ea_totals <- as.data.table(ea_totals)

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

cluster_alloc <- cluster_alloc %>% pivot_wider(id_cols = NAME_1, names_from = urban, values_from = EAs) %>% mutate(EAs = `1`+`0`)
cluster_alloc <- merge(cluster_alloc[,c('NAME_1','EAs')], unique.array(admin.key[,c('admin1','admin1.char','NAME_1')]))
cluster_alloc <- as.data.table(cluster_alloc)

## urban is oversampled in most (90%) areas
# merge(ea_totals %>% group_by(NAME_1) %>% summarise(urb_frac_pop = sum(urban*EAs)/sum(EAs)),cluster_alloc %>% 
#         group_by(NAME_1) %>% summarise(urb_frac_sample = sum(urban*EAs)/sum(EAs))) %>% 
#   ggplot() + geom_point(aes(urb_frac_pop,urb_frac_sample)) + geom_abline(intercept = 0,slope = 1)

# sort population into clusters -----------
cluster_frame <- ea_totals[rep(1:.N, EAs)]
cluster_frame[, cluster := seq_len(.N)]
cluster_frame[,EAs := NULL]
cluster_frame[,strata := admin1]

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

cluster_frame[,admin2_strata := admin2]

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
sim_pop_long <- merge(cluster_frame[,c('cluster','strata','admin2_strata')],sim_pop_long)

cluster_frame <- merge(cluster_frame,sim_pop_long[,.N,cluster],by='cluster')

# objects for sampling ---------
cluster_frame$n <-  rpois(nrow(cluster_frame),n_per_EA) #number that will be sampled from that cluster if sampled
cluster_frame[cluster_frame$n<1,]$n <- 1
#cluster_frame$n <-  n_per_EA

# add population size per admin1
cluster_alloc <- cluster_alloc %>% mutate(strata_total =  pop.admin1$admin_pop[admin1])

cluster_frame <- merge(cluster_frame,cluster_alloc,by=c('admin1','NAME_1','admin1.char')) %>% mutate(wt = strata_total/(EAs*n))
cluster_frame[,EAs:=NULL]
cluster_frame[,strata_total:=NULL]

cluster_alloc$strata <- cluster_alloc$admin1

# functions/models -----
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

setwd('/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing')
mod0 <- cmdstan_model("Stan/NonStrat_Smooth.stan")

hyperpc.iid = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
hyperpc.bym2 = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                    phi = list(prior = "pc", param = c(0.5, 0.5)))

# load covariates ---------
load(file='Kenya_Example/KEN_Covariates/Kenya_admin2_covariates.rda')

# mean covariates: density_log, nt_lights_log, tthc_log, precip, temp, elev
# var covariates: pop_var_log, nt_lights_var, tthc_var, precip_var, temp_var, elev_var, area_log

# simulation parameters and generate surface ------
sig_u <- c(1,0.2)
phi <- c(0.7,0.25)

W1 <- sig_u[1]^2*((1-phi[1])*diag(1,n_admin2) + phi[1]*Q2_scaled_inv)
W2 <- sig_u[2]^2*((1-phi[2])*diag(1,n_admin2) + phi[2]*Q2_scaled_inv)

beta <- c(-1, -0.1, -0.25, 0.1)
gamma <- c(0.5, 0.05, 0.01, -0.01, 0.05)

mean_covariates <- c('nt_lights_log','elev','temp')
var_covariates <- c('nt_lights_var_log','elev_var_log','temp_var_log','area_log')

###

mean_model_matrix <- cbind(rep(1,n_admin2),cmat_admin2[,mean_covariates])
var_model_matrix <- cbind(rep(1,n_admin2),cmat_admin2[,var_covariates])

admin2_mean <- beta %*% t(mean_model_matrix) + as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),W1))
sig2 <- exp(gamma %*% t(var_model_matrix) + as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),W2)))

sim_pop_long[, value := admin2_mean[admin2] + rnorm(.N,0,sqrt(sig2[admin2]))]

plot(sim_pop_long[order(admin2),mean(value),by=admin2]$V1,admin2_mean)
plot(sim_pop_long[order(admin2),var(value),by=admin2]$V1,sig2)

# take sample and get direct estimates --------
direct <- sampled_clusters <- list()
for(k in 1:100){
  
  cat(k,'\n')
  
  # which clusters will we sample?
  cluster_sample_ids <- unlist(sapply(unique(cluster_frame$strata),function(i){
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
  sampled_clusters[[k]] <- sim_sample %>% group_by(admin1,admin2,cluster) %>% reframe(n=n(),wt=unique(wt))
  
  my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"),strata = ~admin1,
                                    weights = ~wt, data = sim_sample)
  
  x <- mapply(admin2.HT.withNA, which.area = 1:n_admin2)
  admin.dir <- data.frame(t(x))
  colnames(admin.dir) <- c('admin2','mean','variance')
  direct[[k]] <- admin.dir
}

# get FH estimates ------
results <- sig2_results <- NULL
for(k in 1:100){
  cat(k,'\n')
  
  dir.dat <- direct[[k]]
  dir.dat <- merge(dir.dat,cmat_admin2)
  
  change_id <- which(dir.dat$variance < 1e-10)
  if(length(change_id)>0){
    dir.dat[change_id,]$variance <- NA
    dir.dat[change_id,]$mean <- NA
  }
  
  # no variance smoothing -----------
  formula_str <- paste0("mean ~ ", paste(mean_covariates, collapse = " + "), " + f(admin2,model='iid',hyper=hyperpc.iid)")
  formula <- as.formula(formula_str)
  
  fit1 <- inla(formula,
               family = "gaussian",
               data = dir.dat,
               quantiles = c(0.05,0.95),
               scale = 1 / dir.dat$variance,
               control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
               control.compute = list(config = TRUE))
  
  # get input for models --------
  data_areas <- which(!is.na(dir.dat$mean))
  Cons_exact <- Cons_naive <- df_exact <- df_naive <- v_hat_scaled_exact <- v_hat_scaled_naive <- rep(NA,n_admin2)
  for(area in data_areas){
    tmp <- sampled_clusters[[k]] %>% filter(admin1 == admin.key[admin.key$admin2 == area,]$admin1)
    tmp[tmp$admin2 != area,]$wt <- 0
    N_D <- sum(tmp$wt>0)
    N <- nrow(tmp)
    
    omega <- tmp$n*tmp$wt
    M <- diag(1,N) - (rep(1,N)%*%t(omega)/sum(omega))
    D <- diag(omega^2)
    S <- diag(1/tmp$n)
    out <- eigen(sqrt(S)%*%t(M)%*%D%*%M%*%sqrt(S))
    q <- out$values
    
    Cons_naive[area] <- 1/N_D
    df_naive[area] <- N_D - 1
    v_hat_scaled_naive[area] <- dir.dat$variance[area]*(N_D^2*(N-1))/N
    
    Cons_exact[area] <- sum(tmp$wt^2*tmp$n)/(sum(omega)^2)
    df_exact[area] <- sum(q)^2/sum(q^2)
    v_hat_scaled_exact[area] <- dir.dat$variance[area]*(sum(omega)^2*(N-1))/N*sum(q)/sum(q^2)
   
  }
  
  # Naive smoothing ------------
  
  data_list2 = list(m=n_admin2,
                    m_data=length(data_areas),
                    data_areas = data_areas,
                    y=dir.dat$mean[data_areas],
                    v_hat_scaled = v_hat_scaled_naive[data_areas],
                    Cons=Cons_naive[data_areas], df=df_naive[data_areas],
                    p_mean = ncol(mean_model_matrix),
                    X = mean_model_matrix,
                    p_var = ncol(var_model_matrix),
                    Z = var_model_matrix,
                    bym2_mean = 1,
                    bym2_var = 1,
                    N_edges = nrow(nodes2),
                    node1 = nodes2$node1,
                    node2 = nodes2$node2,
                    car_scale = Q2_scaled[1,1]/Q.admin2[1,1]) 
  
  fit2 <- mod0$sample(
    data = data_list2,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    show_messages = F,
    show_exceptions = F,
    refresh=0)
  
  theta_fit2 <- as.matrix(unclass(fit2$draws(variables = c("theta"), format = "matrix")))
  sig2_fit2 <- exp(as.matrix(unclass(fit2$draws(variables = c("log_sig2"), format = "matrix"))))
  
  # "Exact" smoothing ----------
  
  data_list3 = list(m=n_admin2,
                    m_data=length(data_areas),
                    data_areas = data_areas,
                    y=dir.dat$mean[data_areas],
                    v_hat_scaled = v_hat_scaled_exact[data_areas],
                    Cons=Cons_exact[data_areas], df=df_exact[data_areas],
                    p_mean = ncol(mean_model_matrix),
                    X = mean_model_matrix,
                    p_var = ncol(var_model_matrix),
                    Z = var_model_matrix,
                    bym2_mean = 1,
                    bym2_var = 1,
                    N_edges = nrow(nodes2),
                    node1 = nodes2$node1,
                    node2 = nodes2$node2,
                    car_scale = Q2_scaled[1,1]/Q.admin2[1,1]) 
  
  fit3 <- mod0$sample(
    data = data_list3,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    show_messages = F,
    show_exceptions = F,
    refresh=0)
  
  theta_fit3 <- as.matrix(unclass(fit3$draws(variables = c("theta"), format = "matrix")))
  sig2_fit3 <- exp(as.matrix(unclass(fit3$draws(variables = c("log_sig2"), format = "matrix"))))

  # Save ------------------------
  results[[k]] <- data.frame(sim=k,
                             area = 1:n_admin2,
                             model = rep(c('std','naive','exact'),each=n_admin2),
                             mean = c(fit1$summary.fitted.values$mean,
                                      apply(theta_fit2,2,mean),
                                      apply(theta_fit3,2,mean)),
                             upper = c(fit1$summary.fitted.values$`0.95quant`,
                                       apply(theta_fit2,2,quantile,probs=0.95),
                                       apply(theta_fit3,2,quantile,probs=0.95)),
                             lower =  c(fit1$summary.fitted.values$`0.05quant`,
                                        apply(theta_fit2,2,quantile,probs=0.05),
                                        apply(theta_fit3,2,quantile,probs=0.05)))
  
  sig2_results[[k]] <- data.frame(sim=k,
                                  area = 1:n_admin2,
                                  model = rep(c('naive','exact'),each=n_admin2),
                                  mean = c(apply(sig2_fit2,2,mean),
                                           apply(sig2_fit3,2,mean)),
                                  upper = c(apply(sig2_fit2,2,quantile,probs=0.95),
                                            apply(sig2_fit3,2,quantile,probs=0.95)),
                                  lower =  c(apply(sig2_fit2,2,quantile,probs=0.05),
                                             apply(sig2_fit3,2,quantile,probs=0.05)))
  
}

# look at results -------
results_mat <- do.call(rbind,results)
results_mat$true_mean <- admin2_mean[results_mat$area]

results_mat <- results_mat %>% mutate(error = mean - admin2_mean[area],
                                      cov = (upper > admin2_mean[area]) &  (lower < admin2_mean[area]),
                                      width = upper - lower)
#results_mat %>% group_by(area,model) %>% reframe(rmse = sqrt(mean(error^2)),coverage = mean(cov),ci_width = mean(width))

results_mat %>% group_by(true_mean,area,model) %>% summarise(coverage = mean(cov)) %>%
  ggplot() + geom_line(aes(x=true_mean,y=coverage,col=model))

results_mat %>% group_by(model) %>% reframe(rmse = sqrt(mean(error^2)),coverage = mean(cov),ci_width = mean(width))

sig2_results_mat <- do.call(rbind,sig2_results)
t <- sig2_results_mat %>% filter(model=='exact') %>% group_by(area) %>% summarise(val=mean(mean))
plot(sig2, t$val)
abline(0,1)

t <- sig2_results_mat %>% filter(model=='naive') %>% group_by(area) %>% summarise(val=mean(mean))
p <- sim_pop_long %>% group_by(admin2,cluster) %>% summarise(val = mean(value)) %>% group_by(admin2) %>% summarise(val = var(val))

# testing ------
plot(sig2,colMeans(sig2_fit3))


t <- sim_pop_long %>% group_by(admin2) %>% summarise(val=var(value))
hist(t$val)

var(sim_pop_long$value)


t <- sim_pop_long %>% group_by(admin2,cluster) %>% reframe(y=mean(value),n=n())
t <- t %>% group_by(admin2) %>% summarise(val = var(y*sqrt(n)))
plot(sig2,t$val)
abline(0,1)
