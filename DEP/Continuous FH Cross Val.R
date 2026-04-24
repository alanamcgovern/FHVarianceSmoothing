library(rdhs)
library(haven)
library(surveyPrev)
library(sf)
library(spdep)
library(tidyverse)
library(INLA)
library(igraph)
library(ggpubr)
source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")

# get polygons and cluster info -----
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

# load covariates ----------------
load(file = "/Users/alanamcgovern/Desktop/Research/NGA_Covariates/Nigera_admin1_covariates.rda")
load(file = "/Users/alanamcgovern/Desktop/Research/NGA_Covariates/Nigera_admin2_covariates.rda")

# rescale some covariates
cmat_admin2 <- cmat_admin2 %>% mutate(elev = (elev - mean(elev))/sd(elev),
                                      precip = (precip - mean(precip))/sd(precip),
                                      temp = (temp - mean(temp))/sd(temp))

cmat_admin1 <- cmat_admin1 %>% mutate(elev = (elev - mean(elev))/sd(elev),
                                      precip = (precip - mean(precip))/sd(precip),
                                      temp = (temp - mean(temp))/sd(temp))

mean_covariates <- c('nt_lights_log','tthc_log','precip','temp','elev','density_log')
var_covariates <- c('nt_lights_var','tthc_var','precip_var','temp_var','elev_var','area','pop_var') ## TAKE THE LOG OF THESE

# helper functions and objects ----------
hyperpc.iid <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

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

# specify outcome and load direct estimates -----
var_t <- 'haz'

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
load(file = paste0('Nigeria_',var_t,'_admin1_weighted_estimates.rda'))
load(file = paste0('Nigeria_',var_t,'_admin2_weighted_estimates.rda'))

admin1.dir <- merge(admin1.dir,cmat_admin1)
admin2.dir <- merge(admin2.dir,cmat_admin2)

admin2.dir.stable <- admin2.dir
admin2.dir.stable$variance <- ifelse(admin2.dir.stable$variance<1e-10,NA,admin2.dir.stable$variance)
admin2.dir.stable$mean <- ifelse(is.na(admin2.dir.stable$variance),NA,admin2.dir.stable$mean)

# CV admin1 ----------------------
area_seq <- admin1.dir$admin1

start.time <- Sys.time()
admin1_out <- lapply(area_seq,function(area_t){
  dir_tmp <- admin1.dir
  dir_tmp[dir_tmp$admin1==area_t,]$mean <- dir_tmp[dir_tmp$admin1==area_t,]$variance <- NA
  
  # just iid
  eta.samples1 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = dir_tmp,
                                        random_effect = "admin1",
                                        #covariates = mean_covariates,
                                        #adj_matrix = admin1.mat,
                                        random_model = 'iid')
  
  # iid + covariates
  eta.samples2 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = dir_tmp,
                                        random_effect = "admin1",
                                        covariates = mean_covariates,
                                        #adj_matrix = admin1.mat,
                                        random_model = 'iid')
  
  
  # bym2
  eta.samples3 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = dir_tmp,
                                        random_effect = "admin1",
                                        #covariates = mean_covariates,
                                        adj_matrix = admin1.mat,
                                        random_model = 'bym2')
  
  # bym2 + covariates
  eta.samples4 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = dir_tmp,
                                        random_effect = "admin1",
                                        covariates = mean_covariates,
                                        adj_matrix = admin1.mat,
                                        random_model = 'bym2')
  
  out <- data.frame(mean = c(mean(eta.samples1[,area_t]),mean(eta.samples2[,area_t]),mean(eta.samples3[,area_t]),mean(eta.samples4[,area_t])),
                    lower90 = c(quantile(eta.samples1[,area_t],c(0.05)),quantile(eta.samples2[,area_t],c(0.05)),quantile(eta.samples3[,area_t],c(0.05)),quantile(eta.samples4[,area_t],c(0.05))),
                    upper90 = c(quantile(eta.samples1[,area_t],c(0.95)),quantile(eta.samples2[,area_t],c(0.95)),quantile(eta.samples3[,area_t],c(0.95)),quantile(eta.samples4[,area_t],c(0.95))),
                    model = c('iid','iid+cov','bym2','bym2+cov'))
  out$direct <- admin1.dir$mean[area_t]
  
  return(out)
  
})
Sys.time() - start.time

cv_admin1_results <- do.call(rbind,admin1_out)
save(cv_admin1_results, file=paste0('Nigeria_',var_t,'_admin1_cv.rda'))
load(paste0('Nigeria_',var_t,'_admin1_cv.rda'))
cv_admin1_results %>% group_by(model) %>% reframe(rmpse = sqrt(mean(mean-direct)^2), 
                                                  width = mean(upper90-lower90),
                                                  cov = mean(upper90>direct & lower90<direct))

# CV admin2 ---------------------
area_seq <- admin2.dir[!is.na(admin2.dir$mean),]$admin2
admin2_out <- list()
#admin2_out <- lapply(area_seq,function(area_t)
for(area_t in area_seq){
  cat(area_t,'\n')
  dir_tmp <- admin2.dir.stable
  dir_tmp[dir_tmp$admin2==area_t,]$mean <- dir_tmp[dir_tmp$admin2==area_t,]$variance <- NA
  
  # just iid
  eta.samples1 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = dir_tmp,
                                        random_effect = "admin2",
                                        #covariates = mean_covariates,
                                        #adj_matrix = admin2.mat,
                                        random_model = 'iid')
  
  # iid + covariates
  eta.samples2 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = dir_tmp,
                                        random_effect = "admin2",
                                        covariates = mean_covariates,
                                        #adj_matrix = admin2.mat,
                                        random_model = 'iid')
  
  
  # bym2
  eta.samples3 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = dir_tmp,
                                        random_effect = "admin2",
                                        #covariates = mean_covariates,
                                        adj_matrix = admin2.mat,
                                        random_model = 'bym2')
  
  # bym2 + covariates
  eta.samples4 <- fit_inla_with_samples(outcome = 'mean',
                                        scale_var = 'variance',
                                        data = dir_tmp,
                                        random_effect = "admin2",
                                        covariates = mean_covariates,
                                        adj_matrix = admin2.mat,
                                        random_model = 'bym2')
  
  out <- data.frame(mean = c(mean(eta.samples1[,area_t]),mean(eta.samples2[,area_t]),mean(eta.samples3[,area_t]),mean(eta.samples4[,area_t])),
                    lower90 = c(quantile(eta.samples1[,area_t],c(0.05)),quantile(eta.samples2[,area_t],c(0.05)),quantile(eta.samples3[,area_t],c(0.05)),quantile(eta.samples4[,area_t],c(0.05))),
                    upper90 = c(quantile(eta.samples1[,area_t],c(0.95)),quantile(eta.samples2[,area_t],c(0.95)),quantile(eta.samples3[,area_t],c(0.95)),quantile(eta.samples4[,area_t],c(0.95))),
                    model = c('iid','iid+cov','bym2','bym2+cov'))
  out$direct <- admin2.dir$mean[area_t]
  
  admin2_out[[area_t]] <- out
  
}#)

cv_admin2_results <- do.call(rbind,admin2_out)
save(cv_admin2_results, file=paste0('Nigeria_',var_t,'_admin2_cv.rda'))

load(file=paste0('Nigeria_',var_t,'_admin2_cv.rda'))

cv_admin2_results %>% group_by(model) %>% reframe(rmpse = sqrt(mean(mean-direct)^2), 
                                                  width = mean(upper90-lower90),
                                                  cov = mean(upper90>direct & lower90<direct))


