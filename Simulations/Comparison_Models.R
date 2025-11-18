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

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing")

# load simulation objects -------------
load(file='Simulations/sim1_sampled_clusters.rda')
load(file='Simulations/sim1_admin1.direct.rda')
load(file='Simulations/sim1_admin2.direct.rda')
load(file='Simulations/sim1_popdata.rda')
load(file='Simulations/sim1_params.rda')

# load geometry ------
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


# load covariates ----------
load(file='KEN_Covariates/Kenya_admin2_covariates.rda')
load(file='KEN_Covariates/Kenya_admin1_covariates.rda')

cmat_admin1 <- merge(cmat_admin1,sim_pop_long[order(admin1),mean(urban),by=admin1]) %>% rename(urb_frac = V1)
cmat_admin2 <- merge(cmat_admin2,sim_pop_long[order(admin2),mean(urban),by=admin2]) %>% rename(urb_frac = V1)

# helper functions and objects ----------
adm2.dir.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct.adm2, function(x){x$mean[i]})),na.rm = T)})
adm1.dir.variance <- sapply(1:n_admin1, function(i){var(unlist(lapply(direct.adm1, function(x){x$mean[i]})),na.rm = T)})

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

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing")
mod0 <- cmdstan_model("Stan/BYM2_FH_NaiveSmooth.stan")

# get FH estimates (comparison models) from direct --------

results.adm1 <- results.adm2 <- list() # record estimates 
#covariates = mean_covariates
covariates = c('urb_frac','temp','elev','tthc_log','nt_lights_log')
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
  sim.admin2.dir <- as.data.frame(sim.admin2.dir)
  data_areas <- sim.admin2.dir[!is.na(sim.admin2.dir$mean),]$admin2
  
  sim.admin1.dir <- direct.adm1[[k]]
  sim.admin1.dir$var.truth <- adm1.dir.variance
  sim.admin1.dir <- merge(sim.admin1.dir,cmat_admin1)
  sim.admin1.dir <- sim.admin1.dir[order(sim.admin1.dir$admin1),]
  sim.admin1.dir <- as.data.frame(sim.admin1.dir)
  
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
  
  # FH no smoothing ----
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
  
  # FH with naive smoothing -----------
  
  Cons <- v_hat_scaled <- df<- rep(NA,n_admin1)
  for(area in 1:n_admin1){
    tmp <- sampled_clusters[[k]] %>% filter(admin1==area) %>% group_by(cluster,urban) %>% reframe(n=n())
    N = c(sum(tmp$urban),sum(1-tmp$urban))
    
    # outputs
    Cons[area] <- 1/(sum(tmp$n)) #multiplier for variance term in weighted mean distribution
    
    # recale to be in terms of chi square approximation
    v_hat_scaled[area] <- sim.admin1.dir$variance[area]*sum(N)^2/mean(N/(N-1))
    
    df[area] <- sum(N) - length(unique(tmp$urban))
  }
  
  data_list = list(m=n_admin1,
                   m_data = n_admin1,
                   data_areas = 1:n_admin1,
                   y=sim.admin1.dir$mean,
                   v_hat_scaled = v_hat_scaled,
                   Cons=Cons,
                   df=df,
                   # covariates in mean model
                   p_mean = length(covariates) + 1,
                   X = cbind(rep(1,n_admin1),sim.admin1.dir[,covariates]),
                   # covariates in log-variance model (none for now)
                   p_var = 1,
                   Z = matrix(1,n_admin1,1),
                   N_edges = nrow(nodes1),
                   node1 = nodes1$node1,
                   node2 = nodes1$node2,
                   car_scale = Q1_scaled[1,1]/Q.admin1[1,1]) 
  
  
  # Fit the model
  fit <- mod0$sample(
    data = data_list,
    chains = 4,
    parallel_chains = 4,
    refresh = 0,
    iter_warmup = 1000,
    iter_sampling = 1000
  )
  
  eta.samples3a <- as.matrix(unclass(fit$draws(variables = c("theta"), format = "matrix")))
  
  Cons <- v_hat_scaled <- df<- rep(NA,length(data_areas))
  for(area in data_areas){
    id <- which(area==data_areas)
    
    tmp <- sampled_clusters[[k]] %>% filter(admin1 == admin.key[admin.key$admin2==area,]$admin1) %>% 
      group_by(cluster,urban,admin2) %>% reframe(n=n())
    N = c(sum(tmp$urban),sum(1-tmp$urban))
    
    tmp_D <- sampled_clusters[[k]] %>% filter(admin2 == area) %>% group_by(cluster,urban) %>% reframe(n=n())
    N_D = c(sum(tmp_D$urban),sum(1-tmp_D$urban))
    
    # outputs 
    Cons[id] <- 1/(sum(tmp_D$n)) #multiplier for variance term in weighted mean distribution
    
    # rescale to be in terms of chi square approximation
    v_hat_scaled[id] <- sim.admin2.dir[sim.admin2.dir$admin2==area,]$variance*sum(N_D)^2/mean(N/(N-1))
    
    df[id] <- max(1,sum(N_D) - length(unique(tmp_D$urban)))
  }
  
  data_list = list(m=n_admin2,
                   m_data = length(data_areas),
                   data_areas = data_areas,
                   y=sim.admin2.dir[!is.na(sim.admin2.dir$mean),]$mean,
                   v_hat_scaled = v_hat_scaled,
                   df=df,
                   # constant for likelihood
                   Cons = Cons,
                   # covariates in mean model
                   p_mean = length(covariates) + 1,
                   X = cbind(rep(1,n_admin2),sim.admin2.dir[,covariates]),
                   # covariates in log-variance model (none for now)
                   p_var = 1,
                   Z = matrix(1,n_admin2,1),
                   N_edges = nrow(nodes2),
                   node1 = nodes2$node1,
                   node2 = nodes2$node2,
                   car_scale = Q2_scaled[1,1]/Q.admin2[1,1]) 
  
  # Fit the model
  fit <- mod0$sample(
    data = data_list,
    chains = 4,
    parallel_chains = 4,
    refresh = 0,
    iter_warmup = 1000,
    iter_sampling = 1000
  )
  
  eta.samples3 <- as.matrix(unclass(fit$draws(variables = c("theta"), format = "matrix")))
  
  ## record estimates -----
  
  results.adm2[[k]] <- merge(sim.admin2.dir[,c('admin2','mean','variance','var.truth')],
                             data.frame(admin2 = 1:n_admin2,
                                        # standard FH
                                        mean.fh.std  = apply(eta.samples2,2,mean),
                                        var.fh.std = apply(eta.samples2,2,var),
                                        lower.fh.std = apply(eta.samples2,2,quantile,prob=0.05),
                                        upper.fh.std = apply(eta.samples2,2,quantile,prob=0.95),
                                        # uses correct variance
                                        mean.fh.oracle  = apply(eta.samples1,2,mean),
                                        var.fh.oracle = apply(eta.samples1,2,var),
                                        lower.fh.oracle = apply(eta.samples1,2,quantile,prob=0.05),
                                        upper.fh.oracle = apply(eta.samples1,2,quantile,prob=0.95),
                                        # naive smoothing
                                        mean.fh.naive  = apply(eta.samples3,2,mean),
                                        var.fh.naive = apply(eta.samples3,2,var),
                                        lower.fh.naive = apply(eta.samples3,2,quantile,prob=0.05),
                                        upper.fh.naive = apply(eta.samples3,2,quantile,prob=0.95),
                                        sim = k,
                                        true_mean = params_list$admin2_means),all=T)
  
  results.adm1[[k]] <- merge(sim.admin1.dir[,c('admin1','mean','variance','var.truth')],
                             data.frame(admin1 = 1:n_admin1,
                                        # standard FH
                                        mean.fh.std  = apply(eta.samples2a,2,mean),
                                        var.fh.std = apply(eta.samples2a,2,var),
                                        lower.fh.std = apply(eta.samples2a,2,quantile,prob=0.05),
                                        upper.fh.std = apply(eta.samples2a,2,quantile,prob=0.95),
                                        # uses correct variance
                                        mean.fh.oracle  = apply(eta.samples1a,2,mean),
                                        var.fh.oracle = apply(eta.samples1a,2,var),
                                        lower.fh.oracle = apply(eta.samples1a,2,quantile,prob=0.05),
                                        upper.fh.oracle = apply(eta.samples1a,2,quantile,prob=0.95),
                                        # naive smoothing
                                        mean.fh.naive  = apply(eta.samples3a,2,mean),
                                        var.fh.naive = apply(eta.samples3a,2,var),
                                        lower.fh.naive = apply(eta.samples3a,2,quantile,prob=0.05),
                                        upper.fh.naive = apply(eta.samples3a,2,quantile,prob=0.95),
                                        sim = k,
                                        true_mean = params_list$admin1_means),all=T)
  
}

# SAVE RESULTS ---------
setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Simulations")
save(results.adm1,file = 'Simulations/sim1_admin1.fh.comparisons.rda')
save(results.adm2,file = 'Simulations/sim1_admin2.fh.comparisons.rda')

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

