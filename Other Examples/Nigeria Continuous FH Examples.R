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

# read in survey data -------

load(file='/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/NGA_2024_cont_outcomes.rda')

# geo <- getDHSgeo(country='Nigeria',year=2024)
# cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, by.adm1 = "NAME_1",by.adm2 = "NAME_2")
# cluster.info$data <- cluster.info$data %>%
#   rename(NAME_1 = admin1.name, NAME_2 = admin2.name) %>%
#   dplyr::select(-c(LONGNUM,LATNUM,geometry,admin2.name.full))
# cluster.info$data <- merge(cluster.info$data,admin.key,by=c('NAME_1', 'NAME_2'))
# 
# set_rdhs_config(email = "amcgov@uw.edu",
#                 project = "Spatial Modeling for Subnational Administrative Level 2 Small-Area Estimation - Under 5 Mortality Rate")
# 
# survey_codes <- dhs_datasets(countryIds = "NG") %>%
#   dplyr::filter(SurveyId == 'NG2018DHS' & FileFormat=='Stata dataset (.dta)')
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
# save(dat,file='/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/NGA_2024_cont_outcomes.rda')


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
  poly.adm1 <- st_transform(poly.adm1, 4326)
}
centroids <- st_centroid(poly.adm1)
coords <- st_coordinates(centroids)

cmat_admin1$long <- coords[,('X')]
cmat_admin1$lat <- coords[,('Y')]

# add log transforms of variance covariates
cmat_admin2 <- cmat_admin2 %>% mutate(nt_lights_var_log = log(nt_lights_var+0.001),
                                      tthc_var_log = log(tthc_var),
                                      precip_var_log = log(precip_var),
                                      temp_var_log = log(temp_var),
                                      elev_var_log = log(elev_var))

cmat_admin1 <- cmat_admin1 %>% mutate(nt_lights_var_log = log(nt_lights_var+0.001),
                                      tthc_var_log = log(tthc_var),
                                      precip_var_log = log(precip_var),
                                      temp_var_log = log(temp_var),
                                      elev_var_log = log(elev_var))


# standardize covariates
cmat_admin2 <- cmat_admin2 %>% mutate(ntlights_log = (nt_lights_log - mean(nt_lights_log))/sd(nt_lights_log),
                                      elev = (elev - mean(elev))/sd(elev),
                                      tthc_log = (tthc_log - mean(tthc_log))/sd(tthc_log),
                                      precip = (precip - mean(precip))/sd(precip),
                                      temp = (temp - mean(temp))/sd(temp),
                                      density_log = (density_log - mean(density_log))/sd(density_log),
                                      long = (long - mean(long))/sd(long),
                                      lat = (lat - mean(lat))/sd(lat),
                                      area_log = (area_log-mean(area_log))/sd(area_log),
                                      pop_var_log = (pop_var_log-mean(pop_var_log))/sd(pop_var_log),
                                      nt_lights_var_log = (nt_lights_var_log-mean(nt_lights_var_log))/sd(nt_lights_var_log),
                                      tthc_var_log = (tthc_var_log-mean(tthc_var_log))/sd(tthc_var_log),
                                      precip_var_log = (precip_var_log-mean(precip_var_log))/sd(precip_var_log),
                                      temp_var_log = (temp_var_log-mean(temp_var_log))/sd(temp_var_log),
                                      elev_var_log = (elev_var_log-mean(elev_var_log))/sd(elev_var_log))


cmat_admin1 <- cmat_admin1 %>% mutate(ntlights_log = (nt_lights_log - mean(nt_lights_log))/sd(nt_lights_log),
                                      elev = (elev - mean(elev))/sd(elev),
                                      tthc_log = (tthc_log - mean(tthc_log))/sd(tthc_log),
                                      precip = (precip - mean(precip))/sd(precip),
                                      temp = (temp - mean(temp))/sd(temp),
                                      density_log = (density_log - mean(density_log))/sd(density_log),
                                      long = (long - mean(long))/sd(long),
                                      lat = (lat - mean(lat))/sd(lat),
                                      area_log = (area_log-mean(area_log))/sd(area_log),
                                      pop_var_log = (pop_var_log-mean(pop_var_log))/sd(pop_var_log),
                                      nt_lights_var_log = (nt_lights_var_log-mean(nt_lights_var_log))/sd(nt_lights_var_log),
                                      tthc_var_log = (tthc_var_log-mean(tthc_var_log))/sd(tthc_var_log),
                                      precip_var_log = (precip_var_log-mean(precip_var_log))/sd(precip_var_log),
                                      temp_var_log = (temp_var_log-mean(temp_var_log))/sd(temp_var_log),
                                      elev_var_log = (elev_var_log-mean(elev_var_log))/sd(elev_var_log))
# start loop -----
setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
var_seq <- c('haz','waz','whz')
var_t <- 'haz'
#for(var_t in var_seq){
 # pdf(paste('Nigeria',var_t,'.pdf'))
# direct estimate --------

dir.dat <- dat[!is.na(dat[,c(var_t)]),]
dir.dat$value <- dir.dat[,c(var_t)]
#hist(dir.dat$value)
#qqnorm(dir.dat$value)
#qqline(dir.dat$value)

options(survey.lonely.psu = "adjust")
my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"), 
                                  strata = ~v023, nest = T, weights = ~v005, data = dir.dat)

## admin1
admin1.HT.withNA <- function(which.area) {
  admin1 <- NULL
  tmp <- subset(my.svydesign, (admin1 == as.character(which.area)))
  
  if (dim(tmp)[1] == 0) {
    return(rep(NA, 2))
  } else {
    lm.ob <- survey::svymean(value ~ 1, design = tmp)
    return(c(which.area, lm.ob[1], vcov(lm.ob)))
  }
}

x <- mapply(admin1.HT.withNA, which.area = 1:n_admin1)
admin1.dir <- data.frame(t(x))
colnames(admin1.dir) <- c('admin1','mean','variance')

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

admin2.dir$bad.var <- ifelse(admin2.dir$variance<1e-10,T,F)

admin2.dir <- merge(admin2.dir,admin.key)
admin2.dir <- merge(admin2.dir,cmat_admin2)
admin2.dir <- admin2.dir[order(admin2.dir$admin2),]

admin2.dir.stable <- admin2.dir
admin2.dir.stable$variance <- ifelse(admin2.dir.stable$variance < 1e-10,NA,admin2.dir.stable$variance)
admin2.dir.stable$mean <- ifelse(is.na(admin2.dir.stable$variance),NA,admin2.dir.stable$mean)

# setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
# save(admin1.dir,file = 'Nigeria_haz_admin1_weighted_estimates.rda')
# save(admin2.dir,file = 'Nigeria_haz_admin2_weighted_estimates.rda')


# get maps for direct estimates -----------
m1 <- clean_map_theme + geom_sf(data = merge(poly.adm1,admin1.dir),aes(fill=mean)) + ggtitle('Admin1') +
  scale_fill_viridis_c()
v1 <- clean_map_theme + geom_sf(data = merge(poly.adm1,admin1.dir),aes(fill=sqrt(variance)))+
  scale_fill_viridis_c(name = 'SE',direction = -1)


m2 <- clean_map_theme + geom_sf(data = merge(poly.adm2,admin2.dir),aes(fill=mean))  + ggtitle('Admin2') +
  scale_fill_viridis_c()

v2 <- clean_map_theme + geom_sf(data = merge(poly.adm2,admin2.dir),aes(fill = sqrt(variance))) +
  geom_sf(data = subset(merge(poly.adm2,admin2.dir), bad.var), fill = "white") +
  scale_fill_viridis_c(name = 'SE',direction = -1) 

p <- ggarrange(plotlist = list(m1,v1,m2,v2))

missing_adm2 <- round(100*sum(is.na(admin2.dir$mean))/n_admin2,2)
lonely_adm2 <- round(100*sum(admin2.dir$bad.var,na.rm=T)/n_admin2,2)
print(annotate_figure(p, 
                      bottom = text_grob(paste("Admin 2 areas with: no data = ", missing_adm2,"%; 1 cluster = ", lonely_adm2,'%'), 
                                         face = "bold", size = 14)))

# admin2 scatter plots with covariates -------
# add long and lat of center point
admin2.dir <- merge(admin2.dir,admin.key)
admin2.dir <- merge(admin2.dir,cmat_admin2)
admin2.dir <- admin2.dir[order(admin2.dir$admin2),]

# compare mean to covariate
plot_list <- lapply(mean_covariates,function(cov){
  admin2.dir$cov <- admin2.dir[,c(cov)]
  p <- admin2.dir %>% ggplot(aes(mean, cov)) +
    geom_point(alpha = 0.5,size=0.5) + theme_minimal() +
    geom_smooth(method = "loess", span = 0.5, se = FALSE, color = "blue") +
    labs(title = cov) 
})

print(ggarrange(plotlist = plot_list,nrow=3,ncol=2))

# fit with lat/long and compare residuals of means to covariates
admin2.dir.stable <- admin2.dir
admin2.dir.stable$variance <- ifelse(admin2.dir.stable$variance < 1e-10, NA, admin2.dir.stable$variance)
admin2.dir.stable$mean <- ifelse(is.na(admin2.dir.stable$variance),NA,admin2.dir.stable$mean)

fit2 <- inla(mean ~ long + lat,
             data = admin2.dir.stable,
             family='gaussian',
             scale = 1/admin2.dir.stable$variance,
             control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))), 
             control.compute = list(config = TRUE))
admin2.dir$mean.residual <- admin2.dir$mean - fit2$summary.fitted.values$mean

plot_list <- lapply(covariates,function(cov){
  admin2.dir$cov <- admin2.dir[,c(cov)]
  p <- admin2.dir %>% ggplot(aes(mean.residual, cov)) +
    geom_point(alpha = 0.5,size=0.5) + theme_minimal() +
    geom_smooth(method = "loess", span = 0.5, se = FALSE, color = "blue") +
    labs(title = cov) 
})

print(ggarrange(plotlist = plot_list,nrow=3,ncol=2))


# compare variance to variance of covariate
covariates <- c('elev_var','precip_var','temp_var','nt_lights_var','tthc_var','area_log')

plot_list <- lapply(covariates,function(cov){
  admin2.dir$cov <- admin2.dir[,c(cov)]
  p <- admin2.dir %>% filter(!bad.var) %>% ggplot(aes(log(variance), log(cov))) +
    geom_point(alpha = 0.5,size=0.5) + theme_minimal() + ylab('log(var(cov))') +
    geom_smooth(method = "loess", span = 0.5, se = FALSE, color = "blue") +
    labs(title = cov)
})

print(ggarrange(plotlist = plot_list,nrow=3,ncol=2))

# admin1 scatter plots with covariates -------
admin1.dir <- merge(admin1.dir,admin.key)
admin1.dir <- merge(admin1.dir,cmat_admin1)

# compare mean to covariate
covariates <- c('elev','precip','temp','nt_lights_log','tthc_log')

plot_list <- lapply(covariates,function(cov){
  admin1.dir$cov <- admin1.dir[,c(cov)]
  p <- admin1.dir %>% ggplot(aes(mean, cov)) +
    geom_point(alpha = 0.5,size=0.5) + theme_minimal() +
    geom_smooth(method = "loess", span = 1, se = FALSE, color = "blue") +
    labs(title = cov) 
})

print(ggarrange(plotlist = plot_list,nrow=3,ncol=2))

# adjust for lat/long and compare mean residual to covariate
fit1 <- inla(mean ~ long + lat,
             data = admin1.dir,
             family='gaussian',
             scale = 1/admin1.dir$variance,
             control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))), 
             control.compute = list(config = TRUE))
admin1.dir$mean.residual <- admin1.dir$mean - fit1$summary.fitted.values$mean

plot_list <- lapply(covariates,function(cov){
  admin1.dir$cov <- admin1.dir[,c(cov)]
  p <- admin1.dir %>% ggplot(aes(mean.residual, cov)) +
    geom_point(alpha = 0.5,size=0.5) + theme_minimal() +
    geom_smooth(method = "loess", span =1, se = FALSE, color = "blue") +
    labs(title = cov) 
})

print(ggarrange(plotlist = plot_list,nrow=3,ncol=2))

# compare variance to variance of covariate
covariates <- c('elev_var','precip_var','temp_var','nt_lights_var','tthc_var','area_log')

plot_list <- lapply(covariates,function(cov){
  admin1.dir$cov <- admin1.dir[,c(cov)]
  p <- admin1.dir %>% ggplot(aes(log(variance), log(cov))) +
    geom_point(alpha = 0.5,size=0.5) + theme_minimal() + ylab('log(var(cov))') +
    geom_smooth(method = "loess", span = 1, se = FALSE, color = "blue") +
    labs(title = cov)
})

print(ggarrange(plotlist = plot_list,nrow=3,ncol=2))


dev.off()


#}
# get variograms for direct estimates -------
# get_variogram <- function(Amat,z,thresh = 50){
#   # Create graph
#   g <- graph_from_adjacency_matrix(I(Amat >0), mode = "undirected")
#   
#   # Compute shortest path distances between all nodes
#   dist_mat <- distances(g)
#   
#   max_dist <- max(dist_mat[is.finite(dist_mat)])  # ignore Inf for disconnected nodes
#   variogram_values <- n <- rep(NA,max_dist)
#   for(d in 1:max_dist){
#     pairs <- which(dist_mat == d, arr.ind = TRUE)
#     pairs <- pairs[pairs[,'row'] < pairs[,'col'], ]
#     if(!is.null(nrow(pairs))){
#       vals <- 0.5 * (z[pairs[,'row']] - z[pairs[,'col']])^2
#       variogram_values[d] <- mean(vals,na.rm=T)
#       n[d] <- sum(!is.na(vals))
#     }
#   }
#   
#   return(variogram_values[n > thresh])
# }
# 
# par(mfrow=c(2,2))
# 
# # ADMIN 1
# out <- get_variogram(admin1.mat,admin1.dir$mean)
# plot(out, ylim = c(0,max(out)),
#      xlab = '', ylab = '',
#      main = paste('Admin1 weighted estimates for', var_t))
# 
# out <- get_variogram(admin1.mat,log(admin1.dir$variance))
# plot(out, ylim = c(0,max(out)),
#      xlab = '', ylab = 'log scale',
#      main = paste('Admin1 variance estimates for', var_t))
# 
# # ADMIN2
# out <- get_variogram(admin2.mat,admin2.dir$mean)
# plot(out, ylim = c(0,max(out,na.rm=T)),
#      xlab = '', ylab = '',
#      main = paste('Admin2 weighted estimates for', var_t))
# 
# var.censored <- admin2.dir$variance
# var.censored[var.censored < 1e-10] <- NA
# out <- get_variogram(admin2.mat,log(var.censored))
# plot(out, ylim = c(0,max(out,na.rm=T)),
#      xlab = '', ylab = 'log scale',
#      main = paste('Admin2 variance estimates for', var_t))


#}

# unit level (for simulation settings) ------------
hyperpc.iid = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

dir.dat <- merge(dir.dat,cmat_admin2)

fit1 <- inla(value ~ factor(admin2), #+ f(cluster,model='iid',hyper = hyperpc.iid),
             data = dir.dat,
             family = 'gaussian')

fit1 <- inla(value ~ f(admin2,model='iid',hyper = hyperpc.iid), #+ f(cluster,model='iid',hyper = hyperpc.iid),
             data = dir.dat,
             family = 'gaussian')
summary(fit1)
# total variance = 1.33
sum(1/fit1$summary.hyperpar$mean)
# cluster + individual effect
sum(1/fit1$summary.hyperpar$mean[c(1,3)])
# rho
1/fit1$summary.hyperpar$mean[1]/sum(1/fit1$summary.hyperpar$mean[c(1,3)])


fit2 <- inla(value ~  ntlights_log + tthc_log + elev + density_log + precip + temp + 
               lat + long + 
               f(admin2,model='iid',hyper = hyperpc.iid) + f(cluster,model='iid',hyper = hyperpc.iid),
             data = dir.dat,
             family = 'gaussian')
summary(fit2)
# total variance
sum(1/fit2$summary.hyperpar$mean)
# cluster + individual effect
sum(1/fit2$summary.hyperpar$mean[c(1,3)])
# rho
1/fit2$summary.hyperpar$mean[1]/sum(1/fit2$summary.hyperpar$mean[c(1,3)])


fit3 <- inla(value ~  precip + temp + 
               factor(admin2) + f(cluster,model='iid',hyper = hyperpc.iid),
             data = dir.dat,
             family = 'gaussian')
# total variance = 1.29
sum(1/fit3$summary.hyperpar$mean)
# cluster + individual effect
sum(1/fit3$summary.hyperpar$mean)
# rho
1/fit3$summary.hyperpar$mean[1]/sum(1/fit3$summary.hyperpar$mean)


## which covariates should we include? -------
hyperpc.iid = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

fh1 <- inla(mean ~ nt_lights_log+tthc_log+precip+temp+elev+density_log+long+lat+ 
              f(admin1,model='iid',hyper=hyperpc.iid),
            family = "gaussian",
            data = admin1.dir,
            scale = 1/admin1.dir$variance,
            quantiles = c(0.05,0.95),
            control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
            control.compute = list(config = TRUE))

admin1.dir$log_variance <- log(admin1.dir$variance)
varfit1 <- inla(log_variance ~ area_log + pop_var_log + nt_lights_var_log + tthc_var_log + precip_var_log + temp_var_log + elev_var_log,
            family = "gaussian",
            data = admin1.dir,
            quantiles = c(0.05,0.95),
            control.compute = list(config = TRUE))

summary(varfit1)

lm1 <- lm(log(variance) ~ area_log + pop_var_log + nt_lights_var_log + tthc_var_log + precip_var_log + temp_var_log + elev_var_log,
   data=admin1.dir)
summary(lm1)


mean_covariates <- c('nt_lights_log','tthc_log','precip','temp','elev','density_log','long','lat')#[c(1,4,5,8)]

var_covariates <- c('area_log','pop_var_log','nt_lights_var_log','tthc_var_log','precip_var_log','temp_var_log','elev_var_log')[c(1,3,5:7)]

## standard FH ------------------

hyperpc.iid = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

fh1 <- inla(mean ~ #nt_lights_log+tthc_log+precip+temp+elev+density_log+long+lat+ 
              f(admin1,model='iid',hyper=hyperpc.iid),
            family = "gaussian",
            data = admin1.dir,
            scale = 1/admin1.dir$variance,
            control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
            control.compute = list(config = TRUE))
samples <- get_inla_samples_local(fh1, 1000)
fh1.est <- data.frame(admin1 = 1:n_admin1,
                      mean = colMeans(samples[,1:n_admin1] + samples[,n_admin1+1]),
                      lower90 = apply(samples[,1:n_admin1] + samples[,n_admin1+1],2,quantile,0.05),
                      upper90 = apply(samples[,1:n_admin1] + samples[,n_admin1+1],2,quantile,0.95))


fh2 <- inla(mean ~ f(admin2,model='iid',hyper=hyperpc.iid),
            family = "gaussian",
            data = admin2.dir.stable,
            scale = 1/admin2.dir.stable$variance,
            control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
            control.compute = list(config = TRUE))
samples <- get_inla_samples_local(fh2, 1000)
fh2.est <- data.frame(admin2 = 1:n_admin2,
                      mean = colMeans(samples[,1:n_admin2] + samples[,n_admin2+1]),
                      lower90 = apply(samples[,1:n_admin2] + samples[,n_admin2+1],2,quantile,0.05),
                      upper90 = apply(samples[,1:n_admin2] + samples[,n_admin2+1],2,quantile,0.95))


## FH with smoothed variance (N-1 approximation) ---------------

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
mod0 <- cmdstan_model("IID_FH_VarSmooth.stan")

## ADMIN1 (one strata) ------------
data_list = list(m=n_admin1,
                 m_data = n_admin1,
                 data_areas = 1:n_admin1,
                 y=admin1.dir$mean,
                 v=admin1.dir$variance*10,
                #  N=(dir.dat %>% group_by(admin1) %>% summarise(val = length(unique(cluster))))$val, # sampled in strata
                #  N_D=(dir.dat %>% group_by(admin1) %>% summarise(val = length(unique(cluster))))$val,
                N=rep(3,n_admin1),
                 N_D=rep(3,n_admin1),
                 # covariates in mean model
                 p_mean = 1,
                 X = matrix(1,n_admin1,1),
                # p_mean = length(mean_covariates) +1,
                # X = cbind(1,admin1.dir[,mean_covariates]),
                 # covariates in log-variance model (none for now)
                  p_var = 1,
                  Z = matrix(1,n_admin1,1)
                ) #sampled in intersection of strata and domain


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

plot(fh1.est$mean,fh1.est.smooth$mean)
abline(0,1)

sig2.draws <- exp(fit$draws(c("log_sig2")))
plot((admin1.dir$variance),posterior::summarise_draws(sig2.draws)$mean/data_list$N_D)
abline(0,1)

plot(fh1.est$upper90 - fh1.est$lower90,fh1.est.smooth$q95 - fh1.est.smooth$q5)

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
                 # p_mean = length(mean_covariates) +1,
                 # X = cbind(1,admin1.dir[,mean_covariates]),
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

plot(fh2.est$mean,fh2.est.smooth$mean)
abline(0,1)

sig2.draws.smooth <- exp(fit$draws(c("log_sig2")))
plot((admin2.dir$variance),posterior::summarise_draws(sig2.draws.satt)$mean*data_list$Cons)
abline(0,1)

plot(fh2.est$upper90 - fh2.est$lower90,fh2.est.smooth.satt$q95 - fh2.est.smooth.satt$q5)
abline(0,1)



## FH with smoothed variance (Satterwhaite) ------------

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
mod <- cmdstan_model("IID_FH_VarSmooth_Satt.stan")

## ADMIN1 -----

Cons <- v_hat_scaled <- df <- rep(NA,n_admin1)
for(area in 1:n_admin1){
  tmp <- dir.dat[dir.dat$admin1==area,] %>% group_by(cluster) %>% reframe(n=n(),wt=unique(v005/1e7))
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
                 m_data = n_admin1,
                 data_areas = 1:n_admin1,
                 y=admin1.dir$mean,
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

plot(fh1.est$mean,fh1.est.smooth.satt$mean)
abline(0,1)

sig2.draws.satt <- exp(fit$draws(c("log_sig2")))
plot((admin1.dir$variance),posterior::summarise_draws(sig2.draws.satt)$mean*data_list$Cons)
abline(0,1)

plot(fh1.est$upper90 - fh1.est$lower90,fh1.est.smooth.satt$q95 - fh1.est.smooth.satt$q5)
abline(0,1)

## ADMIN 2 ----

data_areas <- admin2.dir.stable[!is.na(admin2.dir.stable$mean),]$admin2

Cons <- v_hat_scaled <- df <- rep(NA,length(data_areas))
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
  q <- eigen(S^(1/2)%*%t(M)%*%D%*%M%*%S^(1/2))$values

  # outputs
  Cons[i] <- c(sum(tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2) #multiplier for variance term in weighted mean distribution

  # rescale to be in terms of chi square approximation
  v_hat_scaled[i] <- sum(q)/(sum(q^2)*N_p/(sum(tmp$wt*tmp$n)^2*(N_p-1)))*admin2.dir.stable[admin2.dir.stable$admin2==data_areas[i],]$variance

  df[i] <- sum(q)^2/sum(q^2)
}

data_list = list(m=n_admin2,
                 m_data = length(data_areas),
                 data_areas = data_areas,
                 y=admin2.dir.stable[!is.na(admin2.dir.stable$mean),]$mean,
                 v_hat_scaled = v_hat_scaled,
                 # df for chi square approx
                 df = df,
                 # constant for likelihood
                 Cons = Cons,
                 # covariates in mean model
                 p_mean = 1,
                 X = matrix(1,n_admin2,1),
                 # p_mean = length(mean_covariates) +1,
                 # X = cbind(1,admin1.dir[,mean_covariates]),
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

sig2.draws.satt <- exp(fit$draws(c("log_sig2")))
plot((admin2.dir$variance),posterior::summarise_draws(sig2.draws.satt)$mean*data_list$Cons)
abline(0,1)

plot(fh2.est$upper90 - fh2.est$lower90,fh2.est.smooth.satt$q95 - fh2.est.smooth.satt$q5)
abline(0,1)


## FH with smoothed variance (Pearson approximation) -- convergence issues ---------------
# problem for areas where 4/delta_3^2 is less than 1 
# v_raw goes negative for large sigma^2
# try with satterwhaite

setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
mod <- cmdstan_model("IID_FH_VarSmooth_Pearson.stan")

Cons <- matrix(NA,n_admin1,2)
delta <- matrix(NA,n_admin1,3)
A <- B <- df <- rep(NA,n_admin1)
for(area in 1:n_admin1){
  tmp <- dir.dat[dir.dat$admin1==area,] %>% group_by(cluster) %>% reframe(n=n(),wt=unique(v005/1e7))
  N = nrow(tmp)
  
  omega <- tmp$n*tmp$wt ## Nx1
  D = diag(omega^2)
  M = diag(1,N) - (rep(1,N) %*% t(omega))/as.numeric(t(rep(1,N)) %*% t(t(omega)))
  S = diag(1/tmp$n) # structure of the covariance matrix
  q <- eigen(S^(1/2)%*%t(M)%*%D%*%M%*%S^(1/2))$values
  
  # outputs
  Cons[area,] <- c(sum(tmp$n*tmp$wt^2)/sum(tmp$n*tmp$wt)^2, #multiplier for variance term in weighted mean distribution
                   N/(sum(tmp$wt*tmp$n)^2*(N-1))) #multiplier for variance in sampling variance distribution
  
  delta[area,] <- c(sum(q),
             sqrt(2*sum(q^2)),
             sqrt(8)*sum(q^3)/(sum(q^2)^(3/2)))
  
  A[area] <- 4/(Cons[area,2]*delta[area,2]*delta[area,3])
  B[area] <- -4*delta[area,1]/(delta[area,2]*delta[area,3]) + 8/delta[area,3]^2
  df[area] <- 8/delta[area,3]^2
}

## ADMIN1 (one strata) 
keep_id <- (1:n_admin1)[!(1:n_admin1 %in% c(17,25))]

data_list = list(m=n_admin1,
                 y=admin1.dir$mean,
                 v_hat=admin1.dir$variance,
                 # constants for approximation
                 Cons = Cons[,1],
                 A=A,
                 B=B,
                 df=df,
                 # covariates in mean model
                 p_mean = 1,
                 X = matrix(1,n_admin1,1),
                 # p_mean = length(mean_covariates) +1,
                 # X = cbind(1,admin1.dir[,mean_covariates]),
                 # covariates in log-variance model (none for now)
                 p_var = 1,
                 Z = matrix(1,n_admin1,1)) 



# Fit the model
fit <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)


draws <- fit$draws(c("theta"))
fh1.est.smooth <- posterior::summarise_draws(draws)

plot(fh1.est$mean,fh1.est.smooth$mean)
abline(0,1)

sig2.draws <- exp(fit$draws(c("log_sig2")))
plot((admin1.dir$variance),posterior::summarise_draws(sig2.draws)$mean*data_list$Cons[,1])
abline(0,1)

plot(fh1.est$upper90 - fh1.est$lower90,fh1.est.smooth$q95 - fh1.est.smooth$q5)
abline(0,1)

