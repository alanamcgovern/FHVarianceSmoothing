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
setting <- 2

load(file=paste0('Simulations/Sim',setting,'/direct.rda'))
load(file=paste0('Simulations/Sim',setting,'/params.rda'))
load(file=paste0('Simulations/Sim',setting,'/cmat_admin2.rda'))

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

admin2.mat <- nb2mat(poly2nb(poly.adm2), zero.policy = TRUE)
colnames(admin2.mat) <- rownames(admin2.mat) <- poly.adm2$admin2.char

# helper functions and objects ----------
adm2.dir.variance <- sapply(1:n_admin2, function(i){var(unlist(lapply(direct, function(x){x$mean[i]})),na.rm = T)})

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
                                  hyperpc.iid = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01))),
                                  hyperpc.bym2 = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)),
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

# get FH estimates (comparison models) from direct --------

results.adm2 <- list() # record estimates 
mean_covariates = c('urb_frac','X1','X2','X3')
random_model = c("iid","bym2")[2]
for(k in 1:100){
  cat(k,'\n')
  
  sim.admin2.dir <- direct[[k]]
  sim.admin2.dir$var.truth <- adm2.dir.variance
  sim.admin2.dir <- sim.admin2.dir %>% mutate(variance = ifelse(variance < 1e-5,NA,variance))
  sim.admin2.dir <- sim.admin2.dir %>% mutate(mean = ifelse(is.na(variance),NA,mean),
                                              var.truth = ifelse(is.na(variance),NA,var.truth))
  sim.admin2.dir <- merge(sim.admin2.dir,cmat_admin2)
  sim.admin2.dir <- sim.admin2.dir[order(sim.admin2.dir$admin2),]
  sim.admin2.dir <- as.data.frame(sim.admin2.dir)
  data_areas <- sim.admin2.dir[!is.na(sim.admin2.dir$mean),]$admin2
  
  ## model is specified with correct variance ----
  eta.samples1 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'var.truth',
                                        data = sim.admin2.dir,
                                        random_effect = "admin2",
                                        covariates = mean_covariates,
                                        adj_matrix = admin2.mat,
                                        random_model = random_model)
  
  # FH no smoothing ----
   eta.samples2 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = sim.admin2.dir,
                                        random_effect = "admin2",
                                        covariates = mean_covariates,
                                        adj_matrix = admin2.mat,
                                        random_model = random_model)
  
  ## record estimates -----
  
  results.adm2[[k]] <- rbind(data.frame(mean  = apply(eta.samples2,2,mean),
                                  variance = apply(eta.samples2,2,var),
                                  lower = apply(eta.samples2,2,quantile,prob=0.05),
                                  upper = apply(eta.samples2,2,quantile,prob=0.95),
                                  area = 1:n_admin2,
                                  sim=k,
                                  model='std'),
                             data.frame(mean  = apply(eta.samples1,2,mean),
                                        variance = apply(eta.samples1,2,var),
                                        lower = apply(eta.samples1,2,quantile,prob=0.05),
                                        upper = apply(eta.samples1,2,quantile,prob=0.95),
                                        area = 1:n_admin2,
                                        sim=k,
                                        model='oracle'))

}

# SAVE RESULTS ---------
setwd(paste0("/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Simulations/Sim",setting))
results.adm2 <- do.call(rbind,results.adm2)

save(results.adm2,file = paste0('admin2.fh.comparisons.rda'))

