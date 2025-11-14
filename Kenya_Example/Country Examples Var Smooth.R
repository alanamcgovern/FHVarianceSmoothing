library(rdhs)
library(haven)
library(surveyPrev)
library(sf)
library(spdep)
library(tidyverse)
library(INLA)
#library(igraph)
library(ggpubr)
library(cmdstanr)
library(expm)
source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")

## I like Kenya the best from what I see here (least shrinkage of the mean)
country <- c('Nigeria','Angola','Ethiopia','Kenya','Namibia','Uganda')[4]

## choose country --------
if(country=='Nigeria'){ # very high variance estimates (about 5x higher) and extreme shrinkage of mean
  country.abbrev <- 'NGA'
  dhs.abbrev <- 'NG'
  survey_year <- 2024
}else if(country=='Angola'){ # same as Nigeria but less dramatic
  country.abbrev <- 'AGO'
  dhs.abbrev <- 'AO'
  survey_year <- 2015
}else if(country=='Ethiopia'){ # similar to Angola
  country.abbrev <- 'ETH'
  dhs.abbrev <- 'ET'
  survey_year <- 2019
}else if(country=='Kenya'){ # similar to Angola and Ethiopia
  country.abbrev <- 'KEN'
  dhs.abbrev <- 'KE'
  survey_year <- 2022
}else if(country=='Namibia'){ # extreme shrinkage of mean (similar to Nigeria)
  country.abbrev <- 'NAM'
  dhs.abbrev <- 'NM'
  survey_year <- 2013
}else if(country=='Uganda'){ # extreme shrinkage of mean (similar to Nigeria)
  country.abbrev <- 'UGA'
  dhs.abbrev <- 'UG'
  survey_year <- 2016
}else{
  message('Add country info')
}


# get polygons (country specific) -----
setwd('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm')
poly.adm2 <- st_read(dsn = paste0('gadm41_',country.abbrev,'_shp'), layer = paste0("gadm41_",country.abbrev,"_2"), options = "ENCODING=UTF-8")
if(country=='Uganda')
  poly.adm2 <- poly.adm2[poly.adm2$ENGTYPE_2!='Water body',]

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


# read in survey data (country specific) -------

load(file=paste0('/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/',country.abbrev,'_',survey_year,'_cont_outcomes.rda'))

# geo <- getDHSgeo(country=country,year=survey_year)
# cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, by.adm1 = "NAME_1",by.adm2 = "NAME_2")
# cluster.info$data <- cluster.info$data %>%
#   rename(NAME_1 = admin1.name, NAME_2 = admin2.name) %>%
#   dplyr::select(-c(LONGNUM,LATNUM,geometry,admin2.name.full))
# cluster.info$data <- merge(cluster.info$data,admin.key,by=c('NAME_1', 'NAME_2'))
# 
# set_rdhs_config(email = "amcgov@uw.edu",
#                 project = "Spatial Modeling for Subnational Administrative Level 2 Small-Area Estimation - Under 5 Mortality Rate")
# 
# survey_codes <- dhs_datasets(countryIds = dhs.abbrev) %>%
#   dplyr::filter(SurveyId == paste0(dhs.abbrev,survey_year,'DHS') & FileFormat=='Stata dataset (.dta)')
# 
# path.children <- get_datasets(survey_codes[survey_codes$FileType=="Children's Recode",]$FileName, clear_cache = T)
# #path.birth <- get_datasets(survey_codes[survey_codes$FileType=="Births Recode",]$FileName, clear_cache = T)
# #paths.household <- get_datasets(survey_codes[survey_codes$FileType=="Household Recode",]$FileName, clear_cache = T)
# 
# raw.dat.tmp <- readRDS(paste0(path.children))
# raw.dat.tmp <- raw.dat.tmp %>%
#   #dplyr::select(find_continuous_vars(raw.dat.tmp)) %>%
#   rename(cluster = v001)
# #get_variable_labels(raw.dat.tmp)
# 
# dat <- merge(raw.dat.tmp,cluster.info$data,by='cluster') %>%
#   mutate(haz = if_else(hw70>9000,NA,hw70/100), # height for age
#          waz = if_else(hw71>9000,NA,hw71/100), # weight for age
#          whz = if_else(hw72>9000,NA,hw72/100), # weight for height
#          hemo = if_else(hw53>900,NA,hw53), # hemoglobin
#          hemo_adj = if_else(hw56>900,NA,hw56)) # hemoglobin adjusted
# 
# save(dat,file=paste0('/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/',country.abbrev,'_',survey_year,'_cont_outcomes.rda'))
# 

# load covariates (country specific) ----------------
setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")

ur_weights_adm1 <- readRDS('Kenya_2022_admin1_u5_ur_weights.rds')
ur_weights_adm2 <- readRDS('Kenya_2022_admin2_u5_ur_weights.rds')

setwd("/Users/alanamcgovern/Desktop/Research/KEN_Covariates")
load(file='Kenya_admin2_covariates.rda')
load(file='Kenya_admin1_covariates.rda')

cmat_admin1 <- merge(cmat_admin1,ur_weights_adm1)
cmat_admin2 <- merge(cmat_admin2,ur_weights_adm2)

# mean covariates: density_log, nt_lights_log, tthc_log, precip, temp, elev, urb_frac
# var covariates: pop_var_log, nt_lights_var, tthc_var, precip_var, temp_var, elev_var, area_log

# specify outcome and get direct estimates ----------
var_t <- 'haz'

dir.dat <- dat[!is.na(dat[,c(var_t)]),]
dir.dat$value <- dir.dat[,c(var_t)]
options(survey.lonely.psu = "adjust")
my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"), 
                                  strata = ~v023, nest = T, weights = ~v005, data = dir.dat)

## admin1
admin1.HT.withNA <- function(which.area) {
  admin1 <- NULL
  tmp <- subset(my.svydesign, (admin1 == as.character(which.area) & v023==1))
  
  if (dim(tmp)[1] == 0) {
    return(c(which.area,rep(NA, 2)))
  } else {
    lm.ob <- survey::svymean(value ~ 1, design = tmp)
    return(c(which.area, lm.ob[1], vcov(lm.ob)))
  }
}

x <- mapply(admin1.HT.withNA, which.area = 1:n_admin1)
admin1.dir <- data.frame(t(x))
colnames(admin1.dir) <- c('admin1','mean','variance')
admin1.dir$admin1.char <- paste0('admin1_',admin1.dir$admin1)
admin1.dir <- merge(admin1.dir,cmat_admin1)
admin1.dir <- admin1.dir[order(admin1.dir$admin1),]

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

x <- mapply(admin2.HT.withNA, which.area = 1:n_admin2)
admin2.dir <- data.frame(t(x))
colnames(admin2.dir) <- c('admin2','mean','variance')
admin2.dir <- merge(admin2.dir,admin.key)
admin2.dir <- merge(admin2.dir,cmat_admin2)
admin2.dir <- admin2.dir[order(admin2.dir$admin2),]


admin2.dir.stable <- admin2.dir
admin2.dir.stable$variance <- ifelse(admin2.dir.stable$variance < 1e-10,NA,admin2.dir.stable$variance)
admin2.dir.stable$mean <- ifelse(is.na(admin2.dir.stable$variance),NA,admin2.dir.stable$mean)

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
save(admin1.dir,file = paste0(country,'_',var_t,'_admin1_weighted_estimates.rda'))
save(admin2.dir,file = paste0(country,'_',var_t,'_admin2_weighted_estimates.rda'))


# get FH estimates (fixed variance) ----------
 
hyperpc.iid = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
hyperpc.bym2 = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                    phi = list(prior = "pc", param = c(0.5, 0.5)))
 
 fh1 <- inla(mean ~ f(admin1,model='iid',hyper=hyperpc.iid) + urb_frac,
             family = "gaussian",
             data = admin1.dir,
             scale = 1/admin1.dir$variance,
             control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
             control.compute = list(config = TRUE))
 # samples <- get_inla_samples_local(fh1, 1000)
 # fh1.est <- data.frame(admin1 = 1:n_admin1,
 #                       mean = colMeans(samples[,1:n_admin1] + samples[,n_admin1+1]),
 #                       lower90 = apply(samples[,1:n_admin1] + samples[,n_admin1+1],2,quantile,0.05),
 #                       upper90 = apply(samples[,1:n_admin1] + samples[,n_admin1+1],2,quantile,0.95))
 summary(fh1)
 
 
 ## tthc, temp, elev
 fh1.cov <- inla(mean ~ f(admin1,model='iid',hyper=hyperpc.iid) + urb_frac + 
                        density_log+ nt_lights_log+ tthc_log+ precip+ temp+ elev,
                      family = "gaussian",
                      data = admin1.dir,
                 quantiles = c(0.05,0.95),
                      scale = 1/admin1.dir$variance,
                      control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                      control.compute = list(config = TRUE))
 summary(fh1.cov)
 
 ## tthc, temp, elev
 fh1.bym2.cov <- inla(mean ~ f(admin1,model='bym2',hyper=hyperpc.bym2,graph = admin1.mat,scale.model = T) + urb_frac + 
                        density_log+ nt_lights_log+ tthc_log+ precip+ temp+ elev,
                  family = "gaussian",
                  data = admin1.dir,
                  quantiles = c(0.05,0.95),
                  scale = 1/admin1.dir$variance,
                  control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                  control.compute = list(config = TRUE))
 summary(fh1.bym2.cov)
 
 var1 <- inla(log(variance) ~ 1,
              data = admin1.dir, family = "gaussian")
 summary(var1)
 
 var1.cov <- inla(log(variance) ~ pop_var_log+ nt_lights_var_log+ tthc_var_log+ precip_var_log+ temp_var_log+ elev_var_log+ area_log,
                  quantiles = c(0.05,0.95),
                  data = admin1.dir, family = "gaussian")
 summary(var1.cov)
 
 # tthc
 var1.bym2 <- inla(log(variance) ~ f(admin1,model='bym2',hyper=hyperpc.bym2,graph = admin1.mat,scale.model = T),
                   data = admin2.dir.stable, family = "gaussian")
 summary(var1.bym2)
 
 #area
 var1.bym2.cov <- inla(log(variance) ~ f(admin1,model='bym2',hyper=hyperpc.bym2,graph = admin1.mat,scale.model = T) +
                         pop_var_log+ nt_lights_var_log+ tthc_var_log+ precip_var_log+ temp_var_log+ elev_var_log+ area_log,
                       quantiles = c(0.05,0.95),
                   data = admin2.dir.stable, family = "gaussian")
 summary(var1.bym2.cov)
 
 
 fh2 <- inla(mean ~ f(admin2,model='iid',hyper=hyperpc.iid) + urb_frac,
             family = "gaussian",
             data = admin2.dir.stable,
             scale = 1/admin2.dir.stable$variance,
             control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
             control.compute = list(config = TRUE))
 # samples <- get_inla_samples_local(fh2, 1000)
 # fh2.est <- data.frame(admin2 = 1:n_admin2,
 #                       mean = colMeans(samples[,1:n_admin2] + samples[,n_admin2+1]),
 #                       lower90 = apply(samples[,1:n_admin2] + samples[,n_admin2+1],2,quantile,0.05),
 #                       upper90 = apply(samples[,1:n_admin2] + samples[,n_admin2+1],2,quantile,0.95))
 summary(fh2)
 
 # density, nt_lights, tthc, temp, elev
 fh2.cov <- inla(mean ~ f(admin2,model='iid',hyper=hyperpc.iid) + urb_frac + 
                        density_log+ nt_lights_log+ tthc_log+ precip+ temp+ elev,
                      family = "gaussian",
                 quantiles = c(0.05,0.95),
                      data = admin2.dir.stable,
                      scale = 1/admin2.dir.stable$variance,
                      control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                      control.compute = list(config = TRUE))
 summary(fh2.cov)
 
 # nt_lights
 fh2.bym2.cov <- inla(mean ~ f(admin2,model='bym2',hyper=hyperpc.bym2,graph = admin2.mat,scale.model = T)+ urb_frac + 
                    density_log+ nt_lights_log+ tthc_log+ 
                      precip+ temp+ elev,
             family = "gaussian",
             quantiles = c(0.05,0.95),
             data = admin2.dir.stable,
             scale = 1/admin2.dir.stable$variance,
             control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
             control.compute = list(config = TRUE))
 summary(fh2.bym2.cov)
 
 var2 <- inla(log(variance) ~ 1,
      data = admin2.dir.stable, family = "gaussian")
 summary(var2)
 
 #area
 var2.cov <- inla(log(variance) ~ pop_var_log+ nt_lights_var_log+ tthc_var_log+ precip_var_log+ temp_var_log+ elev_var_log+ area_log,
                  quantiles = c(0.05,0.95),
              data = admin2.dir.stable, family = "gaussian")
 summary(var2.cov)
 
 #area
 var2.bym2.cov <- inla(log(variance) ~ f(admin2,model='bym2',hyper=hyperpc.bym2,graph = admin2.mat,scale.model = T) +
                         pop_var_log+ nt_lights_var_log+ tthc_var_log+ precip_var_log+ temp_var_log+ elev_var_log+ area_log,
                       quantiles = c(0.05,0.95),
                   data = admin2.dir.stable, family = "gaussian")
 summary(var2.bym2.cov)
 
 
# get FH estimates (smoothed variance w/ N-1 approx) ----------
 
 setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
 mod0 <- cmdstan_model("IID_FH_VarSmooth.stan")
 
 ## ADMIN1 ------------
 data_list = list(m=n_admin1,
                  m_data = n_admin1,
                  data_areas = 1:n_admin1,
                  y=admin1.dir$mean,
                  v=admin1.dir$variance,
                  N=(dir.dat %>% group_by(admin1) %>% summarise(val = length(unique(cluster))))$val, # sampled in strata
                  N_D=(dir.dat %>% group_by(admin1) %>% summarise(val = length(unique(cluster))))$val,
                  # covariates in mean model
                  p_mean = 1,
                  X = matrix(1,n_admin1,1),
                  # covariates in log-variance model (none for now)
                  p_var = 1,
                  Z = matrix(1,n_admin1,1)) 
 
 
 # Fit the model
 fit <- mod0$sample(
   data = data_list,
   chains = 4,
   parallel_chains = 4,
   iter_warmup = 1000,
   iter_sampling = 1000
 )
 
 
 draws <- fit$draws(c("theta"))
 fh1.est.smooth <- posterior::summarise_draws(draws)
 
 ## ADMIN2 (one strata) ---------
 data_areas <- admin2.dir.stable[!is.na(admin2.dir.stable$mean),]$admin2
 
 N_p <- sapply(data_areas, function(area){
   length(unique(dir.dat[dir.dat$admin1==admin.key[admin.key$admin2==area,]$admin1,]$cluster))
 })
 
 N_u <- sapply(data_areas, function(area){
   length(unique(dir.dat[dir.dat$admin2==area,]$cluster))
 })
 
 data_list = list(m=n_admin2,
                  m_data = length(data_areas),
                  data_areas = data_areas,
                  y=admin2.dir.stable[!is.na(admin2.dir.stable$mean),]$mean,
                  v=admin2.dir.stable[!is.na(admin2.dir.stable$mean),]$variance,
                  N=N_p, # sampled in strata
                  N_D=N_u,
                  # covariates in mean model
                  p_mean = 1,
                  X = matrix(1,n_admin2,1),
                  # covariates in log-variance model (none for now)
                  p_var = 1,
                  Z = matrix(1,n_admin2,1)) 
 
 
 # Fit the model
 fit <- mod0$sample(
   data = data_list,
   chains = 4,
   parallel_chains = 4,
   iter_warmup = 1000,
   iter_sampling = 1000
 )
 
 draws <- fit$draws(c("theta"))
 fh2.est.smooth <- posterior::summarise_draws(draws)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 # get FH estimates (smoothed variance w/ satt approx) ----------
 setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
 mod <- cmdstan_model("IID_FH_VarSmooth_Satt.stan")
 
 ## ADMIN1 -----
 
 Cons <- v_hat_scaled<- sum_q <- sum_q2 <- rep(NA,n_admin1)
 for(area in 1:n_admin1){
   tmp <- dir.dat[dir.dat$admin1==area,] %>% group_by(cluster) %>% reframe(n=n(),wt=unique(v005/1e7))
   N = nrow(tmp)
   
   omega <- tmp$n*tmp$wt ## Nx1
   D = diag(omega^2)
   M = diag(1,N) - (rep(1,N) %*% t(omega))/as.numeric(t(rep(1,N)) %*% t(t(omega)))
   S = diag(1/tmp$n) # structure of the covariance matrix
   q <- eigen(sqrtm(S)%*%t(M)%*%D%*%M%*%sqrtm(S))$values
   
   sum_q[area] <- sum(q)
   sum_q2[area] <- sum(q^2)
   
   # outputs
   Cons[area] <- c(sum(tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2) #multiplier for variance term in weighted mean distribution
   
   # recale to be in terms of chi square approximation
   v_hat_scaled[area] <- (N-1)*sum(tmp$wt*tmp$n)^2/N*admin1.dir$variance[area]

 }
 
 data_list = list(m=n_admin1,
                  m_data = n_admin1,
                  data_areas = 1:n_admin1,
                  y=admin1.dir$mean,
                  v_hat_scaled = v_hat_scaled,
                  # constant for likelihood
                  Cons = Cons,
                  # quantities for chi-square approx
                  sum_q = sum_q,
                  sum_q2 = sum_q2,
                  # covariates in mean model
                  p_mean = 1,
                  X = matrix(1,n_admin1,1),
                  # covariates in log-variance model (none for now)
                  p_var = 1,
                  Z = matrix(1,n_admin1,1)) 
 
 # Fit the model
 fit <- mod$sample(
   data = data_list,
   chains = 4,
   parallel_chains = 4,
   iter_warmup = 1000,
   iter_sampling = 1000,
   refresh=0
 )
 
 
 draws <- fit$draws(c("theta"))
 fh1.est.smooth.satt <- posterior::summarise_draws(draws)
 
 plot(fh1$summary.fitted.values$mean,fh1.est.smooth.satt$mean)
 abline(0,1)
 
 sig2.draws.satt <- exp(fit$draws(c("log_sig2")))
 plot((admin1.dir$variance),posterior::summarise_draws(sig2.draws.satt)$mean*data_list$Cons)
 abline(0,1)
 
 plot(fh1.est$upper90 - fh1.est$lower90,fh1.est.smooth.satt$q95 - fh1.est.smooth.satt$q5)
 abline(0,1)
 
 ## ADMIN 2 ----
 
 data_areas <- admin2.dir.stable[!is.na(admin2.dir.stable$mean),]$admin2
 
 Cons <- v_hat_scaled <- sum_q <- sum_q2 <- rep(NA,length(data_areas))
 for(i in 1:length(data_areas)){
   
   # number of clusters in corresponding admin1 area
   N_p <- length(unique(dir.dat[dir.dat$admin1==admin.key[admin.key$admin2==data_areas[i],]$admin1,]$cluster))
   
   tmp <- dir.dat[dir.dat$admin2==data_areas[i],] %>% group_by(cluster) %>% reframe(n=n(),wt=unique(v005/1e7))
   # number of clusters in the admin2 area
   N_u = nrow(tmp)
   
   omega <- tmp$n*tmp$wt ## Nx1
   D = diag(omega^2)
   M = diag(1,N_u) - (rep(1,N_u) %*% t(omega))/as.numeric(t(rep(1,N_u)) %*% t(t(omega)))
   S = diag(1/tmp$n) # structure of the covariance matrix
   q <- eigen(sqrtm(S)%*%t(M)%*%D%*%M%*%sqrtm(S))$values
   
   # outputs
   Cons[i] <- c(sum(tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2) #multiplier for variance term in weighted mean distribution
   
   sum_q[i] <- sum(q)
   sum_q2[i] <- sum(q^2)
   
   # rescale to be in terms of chi square approximation
   v_hat_scaled[i] <- sum(tmp$wt*tmp$n)^2*(N_p-1)/N_p*admin2.dir.stable[admin2.dir.stable$admin2==data_areas[i],]$variance

 }
 
 data_list = list(m=n_admin2,
                  m_data = length(data_areas),
                  data_areas = data_areas,
                  y=admin2.dir.stable[!is.na(admin2.dir.stable$mean),]$mean,
                  v_hat_scaled = v_hat_scaled,
                  # for chi square approx
                  sum_q = sum_q,
                  sum_q2 = sum_q2,
                  # constant for likelihood
                  Cons = Cons,
                  # covariates in mean model
                  p_mean = 1,
                  X = matrix(1,n_admin2,1),
                  # covariates in log-variance model (none for now)
                  p_var = 1,
                  Z = matrix(1,n_admin2,1)) 
 
 # Fit the model
 fit <- mod$sample(
   data = data_list,
   chains = 4,
   parallel_chains = 4,
   iter_warmup = 1000,
   iter_sampling = 1000,
   refresh=0
 )
 
 draws <- fit$draws(c("theta"))
 fh2.est.smooth.satt <- posterior::summarise_draws(draws)
 
 plot(fh2.est$mean,fh2.est.smooth.satt$mean)
 abline(0,1)
 
 sig2.draws.satt <- exp(fit$draws(c("log_sig2_data")))
 plot((admin2.dir.stable[data_areas,]$variance),posterior::summarise_draws(sig2.draws.satt)$mean*data_list$Cons)
 abline(0,1)
 
 plot(fh2.est$upper90 - fh2.est$lower90,fh2.est.smooth.satt$q95 - fh2.est.smooth.satt$q5)
 abline(0,1)
 
 
 # get FH estimates (satt smoothed, accounting for stratification) ---------
 setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
 mod <- cmdstan_model("IID_FH_Satt_Strat.stan")
 
 ## ADMIN1 -----
 
 dir.dat$urban <- ifelse(dir.dat$v025==2,0,dir.dat$v025)
 
 Cons <- v_hat_scaled <-rep(NA,n_admin1)
 q <- q_id <-  nu_main <- nu_urban <- start_q <- NULL
  for(area in 1:n_admin1){
   tmp <- dir.dat[dir.dat$admin1==area,] %>% group_by(urban, cluster) %>% reframe(n=n(),wt=unique(v005/1e7)) %>% arrange(desc(urban))
   
   N <- c(sum(tmp$urban), sum(1-tmp$urban))
   
   omega <- tmp$wt*tmp$n
   
   M <- diag(omega) %*% (diag(1,sum(N)) - (matrix(1,sum(N),ncol=1) %*% t(omega))/sum(omega))
   
   # if both strata have samples
   if(sum(N==0)==0){
     B_comps <- lapply(1:2,function(i){
       if(N[i]>1){
         (N[i]/(N[i]-1))*(diag(1,N[i]) - matrix(1,N[i],1)%*%matrix(1,1,N[i])/N[i])
       }else if(N[i]==1){
         matrix(1,1,1)
       }else{NULL} # shouldnt get here
     })
     
     B <- bdiag(B_comps)
   }else{
     which.strat <- which(N!=0) # which strata has clusters?
     B <- (N[which.strat]/(N[which.strat]-1))*(diag(1,N[which.strat]) - matrix(1,N[which.strat],1)%*%matrix(1,1,N[which.strat])/N[which.strat])
   }
   
   C <- t(M)%*%B%*%M 

   S = diag(1/tmp$n) # structure of the covariance matrix
   out <- eigen(sqrtm(S)%*%C%*%sqrtm(S))
   
   q_tmp <- out$values
   U_tmp <- out$vectors
   keep_id <- which(abs(q_tmp)>1e-10) # only keep non-zero eigenvalues and corresponding eigenvectors
   
   start_q <- c(start_q, length(q)+1)
   q <- c(q,q_tmp[keep_id])
   q_id <- c(q_id,rep(area,length(keep_id)))
   U <- U_tmp[,keep_id]
   
   # outputs
   Cons[area] <- sum(tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2 #multiplier for variance term in weighted mean distribution
   
   # recale to be in terms of chi square approximation
   v_hat_scaled[area] <- sum(tmp$wt*tmp$n)^2*admin1.dir$variance[area]
   
   nu_main <- c(nu_main, as.vector(t(U)%*%sqrtm(S)%*%C%*%matrix(1,nrow=sum(N))))
   nu_urban <- c(nu_urban, as.vector(t(U)%*%sqrtm(S)%*%C%*%matrix(c(rep(1,N[1]),rep(0,N[2])),nrow=sum(N))))

 }
 
 data_list = list(m=n_admin1,
                  m_data = n_admin1,
                  data_areas = 1:n_admin1,
                  len_q = length(q),
                  y=admin1.dir$mean,
                  v_hat_scaled = v_hat_scaled,
                  # for chi square approx
                  q = q,  nu_main = nu_main, nu_urban = nu_urban, q_id = q_id,
                  # number of qs in each area
                  q_per_area = as.vector(table(q_id)),
                  # start points in q for each area
                  start_q =  start_q,
                  # constant for likelihood
                  Cons = Cons,
                  # covariates in mean model
                  p_mean = 2,
                  X = cbind(rep(1,n_admin1),admin1.dir$urb_frac),
                  # covariates in log-variance model (none for now)
                  p_var = 1,
                  Z = matrix(1,n_admin1,1)) 
 
 # Fit the model
 fit <- mod$sample(
   data = data_list,
   chains = 4,
   parallel_chains = 4,
   iter_warmup = 1000,
   iter_sampling = 1000,
   refresh=200
 )
 
theta_draws <- as.matrix(unclass(fit$draws(variables = c("theta"), format = "matrix")))
plot(fh1$summary.fitted.values$mean,colMeans(theta_draws))
abline(0,1)

sig2.draws.satt <- exp(fit$draws(c("log_sig2")))
plot((admin1.dir$variance),posterior::summarise_draws(sig2.draws.satt)$mean*data_list$Cons)
abline(0,1)
 

## ADMIN2 -----------------------

data_areas <- admin2.dir.stable[!is.na(admin2.dir.stable$mean),]$admin2

Cons <- v_hat_scaled <-rep(NA,n_admin2)
q <- q_id <-  nu_main <- nu_urban <- start_q <- NULL
for(area in data_areas){
 # tmp <- dir.dat[dir.dat$admin2==area,] %>% group_by(urban, cluster) %>% reframe(n=n(),wt=unique(v005/1e7)) %>% arrange(desc(urban))
  
  # number of clusters in the larger, planned domain
  tmp <- dir.dat[dir.dat$admin1==admin.key[admin.key$admin2==area,]$admin1,] %>% group_by(admin2, urban, cluster) %>% reframe(n=n(),wt=unique(v005/1e7)) %>% arrange(desc(urban))
  tmp[tmp$admin2 != area,]$wt <- 0
  
  N <- c(sum(tmp$urban), sum(1-tmp$urban))
  
  omega <- tmp$wt*tmp$n
  
  M <- diag(omega) %*% (diag(1,sum(N)) - (matrix(1,sum(N),ncol=1) %*% t(omega))/sum(omega))
  
  # if both strata have samples
  if(sum(N==0)==0){
    B_comps <- lapply(1:2,function(i){
      if(N[i]>1){
        (N[i]/(N[i]-1))*(diag(1,N[i]) - matrix(1,N[i],1)%*%matrix(1,1,N[i])/N[i])
      }else if(N[i]==1){
        matrix(1,1,1)
      }else{NULL} # shouldnt get here
    })
    
    B <- bdiag(B_comps)
  }else{
    which.strat <- which(N!=0) # which strata has clusters?
    B <- (N[which.strat]/(N[which.strat]-1))*(diag(1,N[which.strat]) - matrix(1,N[which.strat],1)%*%matrix(1,1,N[which.strat])/N[which.strat])
  }
  
  C <- t(M)%*%B%*%M 
  
  S = diag(1/tmp$n) # structure of the covariance matrix
  out <- eigen(sqrtm(S)%*%C%*%sqrtm(S))
  
  q_tmp <- out$values
  U_tmp <- out$vectors
  keep_id <- which(abs(q_tmp)>1e-10) # only keep non-zero eigenvalues and corresponding eigenvectors
  
  start_q <- c(start_q, length(q)+1)
  q <- c(q,q_tmp[keep_id])
  q_id <- c(q_id,rep(area,length(keep_id)))
  U <- U_tmp[,keep_id]
  
  # outputs
  Cons[area] <- sum(tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2 #multiplier for variance term in weighted mean distribution
  
  # recale to be in terms of chi square approximation
  v_hat_scaled[area] <- sum(tmp$wt*tmp$n)^2*admin2.dir$variance[area]
  
  nu_main <- c(nu_main, as.vector(t(U)%*%sqrtm(S)%*%C%*%matrix(1,nrow=sum(N))))
  nu_urban <- c(nu_urban, as.vector(t(U)%*%sqrtm(S)%*%C%*%matrix(c(rep(1,N[1]),rep(0,N[2])),nrow=sum(N))))
  
}

data_list = list(m=n_admin2,
                 m_data = length(data_areas),
                 data_areas = data_areas,
                 len_q = length(q),
                 y=admin2.dir$mean[data_areas],
                 v_hat_scaled = v_hat_scaled[data_areas],
                 # for chi square approx
                 q = q,  nu_main = nu_main, nu_urban = nu_urban, q_id = q_id,
                 # number of qs in each area
                 q_per_area = as.vector(table(q_id)),
                 # start points in q for each area
                 start_q =  start_q,
                 # constant for likelihood
                 Cons = Cons[data_areas],
                 # covariates in mean model
                 p_mean = 2,
                 X = cbind(rep(1,n_admin2),admin2.dir$urb_frac),
                 # covariates in log-variance model (none for now)
                 p_var = 1,
                 Z = matrix(1,n_admin2,1)) 

# Fit the model
fit <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh=200
)

theta_draws <- as.matrix(unclass(fit$draws(variables = c("theta"), format = "matrix")))
plot(fh2$summary.fitted.values$mean,colMeans(theta_draws))
abline(0,1)

sig2.draws.satt <- exp(fit$draws(c("log_sig2_data")))
plot((admin2.dir.stable[data_areas,]$variance),
     posterior::summarise_draws(sig2.draws.satt)$mean*data_list$Cons)
abline(0,1)




# BYM2, stratified, satterwhaite -------------
setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
mod <- cmdstan_model("BYM2_FH_Satt_Strat.stan")

####
## ADMIN1 -----

dir.dat$urban <- ifelse(dir.dat$v025==2,0,dir.dat$v025)

v_hat_scaled <-rep(NA,n_admin1)
Cons <- matrix(NA,n_admin1,2)
q <- q_id <-  nu_main <- nu_urban <- start_q <- NULL
kappa_seq <- c(seq(-0.9,4,0.1))
for(area in 1:n_admin1){
  tmp <- dir.dat[dir.dat$admin1==area,] %>% group_by(urban, cluster) %>% reframe(n=n(),wt=unique(v005/1e7)) %>% arrange(desc(urban))
  
  Cons[area,] <- c(sum((1-tmp$urban)*tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2,sum(tmp$urban*tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2) #multiplier for variance term in weighted mean distribution
  
  # recale to be in terms of chi square approximation
  v_hat_scaled[area] <- sum(tmp$wt*tmp$n)^2*admin1.dir$variance[area]
  
  N <- c(sum(tmp$urban), sum(1-tmp$urban))
  
  omega <- tmp$wt*tmp$n
  
  M <- diag(omega) %*% (diag(1,sum(N)) - (matrix(1,sum(N),ncol=1) %*% t(omega))/sum(omega))
  
  # if both strata have samples
  if(sum(N==0)==0){
    B_comps <- lapply(1:2,function(i){
      if(N[i]>1){
        (N[i]/(N[i]-1))*(diag(1,N[i]) - matrix(1,N[i],1)%*%matrix(1,1,N[i])/N[i])
      }else if(N[i]==1){
        matrix(1,1,1)
      }else{NULL} # shouldnt get here
    })
    
    B <- bdiag(B_comps)
  }else{
    which.strat <- which(N!=0) # which strata has clusters?
    B <- (N[which.strat]/(N[which.strat]-1))*(diag(1,N[which.strat]) - matrix(1,N[which.strat],1)%*%matrix(1,1,N[which.strat])/N[which.strat])
  }
  
  C <- t(M)%*%B%*%M 
  
  S = diag(1/tmp$n) # structure of the covariance matrix
  
  rank = nrow(tmp)-sum(N>0)
  q_id <- c(q_id,rep(area,rank))
  start_q <- c(start_q,max(nrow(q) +1,1))
  q_tmp <- nu_main_tmp <- nu_urban_tmp <- matrix(NA,rank,length(kappa_seq))
  for(kappa in kappa_seq){
    kappa_id <- which(kappa==kappa_seq)
    
    perm <- (diag(1,nrow(tmp)) + kappa*diag(tmp$urban))
    
    out <- eigen((perm%*%S)^(1/2)%*%C%*%(perm%*%S)^(1/2))
    
    q_tmp[,kappa_id] <- out$values[1:rank]
    
    U <- out$vectors[,1:rank]
    nu_main_tmp[,kappa_id] <- as.vector(t(U)%*%sqrtm(perm%*%S)%*%C%*%matrix(1,nrow=sum(N))) 
    nu_urban_tmp[,kappa_id] <- as.vector(t(U)%*%sqrtm(perm%*%S)%*%C%*%matrix(c(rep(1,N[1]),rep(0,N[2])),nrow=sum(N)))
    
  }
  
  q <- rbind(q,q_tmp)
  nu_main <- rbind(nu_main, nu_main_tmp)
  nu_urban <- rbind(nu_urban, nu_urban_tmp)
  
}

data_list = list(m=n_admin1,
                 m_data = n_admin1,
                 data_areas = 1:n_admin1,
                 len_q = nrow(q),
                 y=admin1.dir$mean,
                 v_hat_scaled = v_hat_scaled,
                 len_kappa_seq = length(kappa_seq),
                 kappa_seq=kappa_seq,
                 # for chi square approx
                 q = q,  nu_main = nu_main, nu_urban = nu_urban, q_id = q_id,
                 # number of qs in each area
                 q_per_area = as.vector(table(q_id)),
                 # start points in q for each area
                 start_q =  start_q,
                 # constant for likelihood
                 Cons = Cons,
                 # covariates in mean model
                 p_mean = 2,
                 X = cbind(rep(1,n_admin1),admin1.dir$urb_frac),
                 # covariates in log-variance model (none for now)
                 p_var = 1,
                 Z = matrix(1,n_admin1,1),
                 N_edges = nrow(nodes1),
                 node1 = nodes1$node1,
                 node2 = nodes1$node2,
                 car_scale = Q1_scaled[1,1]/Q.admin1[1,1]) 

# Fit the model
fit <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh=200
)

theta_draws <- as.matrix(unclass(fit$draws(variables = c("theta"), format = "matrix")))
plot(fh1$summary.fitted.values$mean,colMeans(theta_draws))
abline(0,1)

sig2.draws.satt <- exp(fit$draws(c("log_sig2")))
kappa.draws <- as.matrix(unclass(fit$draws(variables = c("kappa"), format = "matrix")))
plot((admin1.dir$variance),posterior::summarise_draws(sig2.draws.satt)$mean*data_list$Cons[,1] + (1+mean(kappa_draws))*data_list$Cons[,2])
abline(0,1)

## ADMIN2 -----------------------

data_areas <- admin2.dir.stable[!is.na(admin2.dir.stable$mean),]$admin2

Cons <- matrix(NA,n_admin2,2)
v_hat_scaled <-rep(NA,n_admin2)
q <- q_id <-  nu_main <- nu_urban <- start_q <- NULL
kappa_seq <- c(seq(-0.9,4,0.1))
for(area in data_areas){
  # clusters in the larger, planned domain
  tmp <- dir.dat[dir.dat$admin1==admin.key[admin.key$admin2==area,]$admin1,] %>% group_by(admin2, urban, cluster) %>% reframe(n=n(),wt=unique(v005/1e7)) %>% arrange(desc(urban))
  tmp[tmp$admin2 != area,]$wt <- 0
  
  Cons[area,] <- c(sum((1-tmp$urban)*tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2,sum(tmp$urban*tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2) #multiplier for variance term in weighted mean distribution
  
  # recale to be in terms of chi square approximation
  v_hat_scaled[area] <- sum(tmp$wt*tmp$n)^2*admin2.dir$variance[area]
  
  
  N <- c(sum(tmp$urban), sum(1-tmp$urban))
  
  omega <- tmp$wt*tmp$n
  
  M <- diag(omega) %*% (diag(1,sum(N)) - (matrix(1,sum(N),ncol=1) %*% t(omega))/sum(omega))
  
  # if both strata have samples
  if(sum(N==0)==0){
    B_comps <- lapply(1:2,function(i){
      if(N[i]>1){
        (N[i]/(N[i]-1))*(diag(1,N[i]) - matrix(1,N[i],1)%*%matrix(1,1,N[i])/N[i])
      }else if(N[i]==1){
        matrix(1,1,1)
      }else{NULL} # shouldnt get here
    })
    
    B <- bdiag(B_comps)
  }else{
    which.strat <- which(N!=0) # which strata has clusters?
    B <- (N[which.strat]/(N[which.strat]-1))*(diag(1,N[which.strat]) - matrix(1,N[which.strat],1)%*%matrix(1,1,N[which.strat])/N[which.strat])
  }
  
  C <- t(M)%*%B%*%M 
  S = diag(1/tmp$n) # structure of the covariance matrix
    
  rank = nrow(tmp[tmp$wt>0,]) - 2 + I(sum(tmp$wt==0)>0)
  q_id <- c(q_id,rep(area,rank))
  start_q <- c(start_q,max(nrow(q) +1,1))
  q_tmp <- nu_main_tmp <- nu_urban_tmp <- matrix(NA,rank,length(kappa_seq))
  for(kappa in kappa_seq){
     kappa_id <- which(kappa==kappa_seq)
     
     perm <- (diag(1,nrow(tmp)) + kappa*diag(tmp$urban))
     
     out <- eigen((perm%*%S)^(1/2)%*%C%*%(perm%*%S)^(1/2))
     
     q_tmp[,kappa_id] <- out$values[1:rank]
     
     U <- out$vectors[,1:rank]
     nu_main_tmp[,kappa_id] <- as.vector(t(U)%*%sqrtm(perm%*%S)%*%C%*%matrix(1,nrow=sum(N))) 
     nu_urban_tmp[,kappa_id] <- as.vector(t(U)%*%sqrtm(perm%*%S)%*%C%*%matrix(c(rep(1,N[1]),rep(0,N[2])),nrow=sum(N)))
     
   }
   
   q <- rbind(q,q_tmp)
   nu_main <- rbind(nu_main, nu_main_tmp)
   nu_urban <- rbind(nu_urban, nu_urban_tmp)
  
}

data_list = list(m=n_admin2,
                 m_data = length(data_areas),
                 data_areas = data_areas,
                 len_q = nrow(q),
                 y=admin2.dir$mean[data_areas],
                 v_hat_scaled = v_hat_scaled[data_areas],
                 len_kappa_seq = length(kappa_seq),
                 kappa_seq=kappa_seq,
                 # for chi square approx
                 q = q,  nu_main = nu_main, nu_urban = nu_urban, q_id = q_id,
                 # number of qs in each area
                 q_per_area = as.vector(table(q_id)),
                 # start points in q for each area
                 start_q =  start_q,
                 # constant for likelihood
                 Cons = Cons[data_areas,],
                 # covariates in mean model
                 p_mean = 2,
                 X = cbind(rep(1,n_admin2),admin2.dir$urb_frac),
                 # covariates in log-variance model (none for now)
                 p_var = 1,
                 Z = matrix(1,n_admin2,1),
                 N_edges = nrow(nodes2),
                 node1 = nodes2$node1,
                 node2 = nodes2$node2,
                 car_scale = Q2_scaled[1,1]/Q.admin2[1,1]) 

# Fit the model
fit <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh=200
)

theta_draws <- as.matrix(unclass(fit$draws(variables = c("theta"), format = "matrix")))
plot(fh2$summary.fitted.values$mean,colMeans(theta_draws))
abline(0,1)

sig2.draws.satt <- exp(fit$draws(c("log_sig2_data")))

kappa.draws <- as.matrix(unclass(fit$draws(variables = c("kappa"), format = "matrix")))
plot((admin2.dir.stable[data_areas,]$variance),posterior::summarise_draws(sig2.draws.satt)$mean*data_list$Cons[,1] + (1+mean(kappa_draws))*data_list$Cons[,2])
abline(0,1)



## properties of the data ------

# are variances different between areas -- yes!
dir.dat %>% group_by(admin1) %>% summarise(v = var(value)/n()) %>% ggplot() +geom_histogram(aes(v))

# is variance different across strata -- yes, within area
dir.dat %>% group_by(v025) %>% summarise(sd = var(value)) 

dir.dat %>% group_by(admin1,v025) %>% summarise(v = var(haz)) %>% pivot_wider(values_from = v, names_from = v025) %>% summarise(val = `2`/`1`)

### testing -------------





