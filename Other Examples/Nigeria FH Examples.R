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

geo <- getDHSgeo(country='Nigeria',year=2018)
cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, by.adm1 = "NAME_1",by.adm2 = "NAME_2")
cluster.info$data <- cluster.info$data %>%
  rename(NAME_1 = admin1.name, NAME_2 = admin2.name) %>%
  dplyr::select(-c(LONGNUM,LATNUM,geometry,admin2.name.full))
cluster.info$data <- merge(cluster.info$data,admin.key,by=c('NAME_1', 'NAME_2'))

# read in survey data -------
# set_rdhs_config(email = "amcgov@uw.edu",
#                 project = "Spatial Modeling for Subnational Administrative Level 2 Small-Area Estimation - Under 5 Mortality Rate")

survey_codes <- dhs_datasets(countryIds = "NG") %>%
  dplyr::filter(SurveyId == 'NG2018DHS' & FileFormat=='Stata dataset (.dta)')

path.children <- get_datasets(survey_codes[survey_codes$FileType=="Children's Recode",]$FileName, clear_cache = T)
#path.birth <- get_datasets(survey_codes[survey_codes$FileType=="Births Recode",]$FileName, clear_cache = T)
#paths.household <- get_datasets(survey_codes[survey_codes$FileType=="Household Recode",]$FileName, clear_cache = T)

raw.dat.tmp <- readRDS(paste0(path.children))
raw.dat.tmp <- raw.dat.tmp %>% 
  dplyr::select(v001,v005,v023,
                h2,h3,h5,h7,h4,h6,h8,h9) %>% 
  rename(cluster = v001, strata = v023, wt = v005)
#get_variable_labels(raw.dat.tmp)

dat <- merge(raw.dat.tmp,cluster.info$data,by='cluster') %>% 
  mutate(bcg = if_else(h2 %in% 1:3,1,ifelse(h2 == 8,NA,h2)),
         dpt1 = if_else(h3 %in% 1:3,1,ifelse(h3 == 8,NA,h3)),
         dpt2 = if_else(h5 %in% 1:3,1,ifelse(h5 == 8,NA,h5)),
         dpt3 = if_else(h7 %in% 1:3,1,ifelse(h7 == 8,NA,h7)),
         polio1 = if_else(h4 %in% 1:3,1,ifelse(h4 == 8,NA,h4)),
         polio2 = if_else(h6 %in% 1:3,1,ifelse(h6 == 8,NA,h6)),
         polio3 = if_else(h8 %in% 1:3,1,ifelse(h8 == 8,NA,h8)),
         measles = if_else(h9 %in% 1:3,1,ifelse(h9 == 8,NA,h9)),
         dpt_all = ifelse(dpt1 + dpt2 + dpt3 == 3,1,ifelse(dpt1 + dpt2 + dpt3 < 3,0,NA)),
         polio_all = ifelse(polio1 + polio2 + polio3 == 3,1,ifelse(polio1 + polio2 + polio3 < 3,0,NA)),
         full_vax = ifelse(dpt_all + polio_all + bcg + measles == 4,1,
                           ifelse(dpt_all + polio_all + bcg + measles < 4,0,NA))) 
dat <- dat[!is.na(dat$full_vax),]

cluster_dat <- dat %>% group_by(admin1,admin1.char,NAME_1,NAME_2,admin2,admin2.char,strata,cluster,wt) %>% reframe(bcg = sum(bcg), measles = sum(measles), full_vax = sum(full_vax),
                                                                                      dpt1 = sum(dpt1),dpt2 = sum(dpt2),dpt3 = sum(dpt3),dpt_all = sum(dpt_all),
                                                                                      polio1 = sum(polio1),polio2 = sum(polio2),polio3 = sum(polio3),polio_all = sum(polio_all),
                                                                                      total = n())
# load covariates ----------------
load(file = "/Users/alanamcgovern/Desktop/Research/NGA_Covariates/Nigera_admin1_covariates.rda")
load(file = "/Users/alanamcgovern/Desktop/Research/NGA_Covariates/Nigera_admin2_covariates.rda")

# standardize covariates
cmat_admin2 <- cmat_admin2 %>% mutate(ntlights_log = (nt_lights_log - mean(nt_lights_log))/sd(nt_lights_log),
                        elev = (elev - mean(elev))/sd(elev),
                        tthc_log = (tthc_log - mean(tthc_log))/sd(tthc_log),
                        precip = (precip - mean(precip))/sd(precip),
                        temp = (temp - mean(temp))/sd(temp))

cmat_admin1 <- cmat_admin1 %>% mutate(ntlights_log = (nt_lights_log - mean(nt_lights_log))/sd(nt_lights_log),
                                      elev = (elev - mean(elev))/sd(elev),
                                      tthc_log = (tthc_log - mean(tthc_log))/sd(tthc_log),
                                      precip = (precip - mean(precip))/sd(precip),
                                      temp = (temp - mean(temp))/sd(temp))


# get direct estimates -------------

indicator_seq <- c('bcg','measles','full_vax','dpt_all','polio_all')
indicator <- indicator_seq[3]

 setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
# pdf("Nigeria Variograms.pdf")
#for(indicator in indicator_seq){

cluster_dat[,'value'] <- cluster_dat[,indicator]

options(survey.lonely.psu = "adjust")
my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"), 
                                  strata = ~strata, nest = T, weights = ~wt, data = cluster_dat)

admin1.HT.withNA <- function(which.area) {
  admin1 <- NULL
  tmp <- subset(my.svydesign, (admin1 == as.character(which.area)))
  
  if (dim(tmp)[1] == 0 | sum(tmp$variables$value) == 0 | sum(tmp$variables$value) == sum(tmp$variables$total)) {
    return(c(which.area,rep(NA, 2)))
  } else {
    glm.ob <- survey::svyglm(value/total ~ 1, design = tmp, 
                             family = stats::quasibinomial, maxit = 50,
                             weights = tmp$variables$total)
    var.est <- stats::vcov(glm.ob)
    logit.mean <- glm.ob$coefficients
    return(c(which.area, logit.mean, var.est))
  }
}

x <- mapply(admin1.HT.withNA, which.area = 1:n_admin1)
admin1.dir <- data.frame(t(x))
colnames(admin1.dir) <- c('admin1','logit.mean','var.logit.est')


admin2.HT.withNA <- function(which.area) {
  admin2 <- NULL
  tmp <- subset(my.svydesign, (admin2 == as.character(which.area)))
  
  if (dim(tmp)[1] == 0 | sum(tmp$variables$value) == 0 | sum(tmp$variables$value) == sum(tmp$variables$total)) {
    return(c(which.area,rep(NA, 2)))
  } else {
    glm.ob <- survey::svyglm(value/total ~ 1, design = tmp, 
                             family = stats::quasibinomial, maxit = 50,
                             weights = tmp$variables$total)
    var.est <- stats::vcov(glm.ob)
    logit.mean <- glm.ob$coefficients
    return(c(which.area, logit.mean, var.est))
  }
}

x <- mapply(admin2.HT.withNA, which.area = 1:n_admin2)
admin2.dir <- data.frame(t(x))
colnames(admin2.dir) <- c('admin2','logit.mean','var.logit.est')

admin2.dir.stable <- admin2.dir[!is.na(admin2.dir$logit.mean),]
admin2.dir.stable <- admin2.dir.stable[admin2.dir.stable$var.logit.est > 1e-10,]

# get maps for direct estimates -----------
m1 <- clean_map_theme + geom_sf(data = merge(poly.adm1,admin1.dir),aes(fill=logit.mean)) + ggtitle(paste(indicator))
v1 <- clean_map_theme + geom_sf(data = merge(poly.adm1,admin1.dir),aes(fill=log(var.logit.est)))

m2 <- clean_map_theme + geom_sf(data = merge(poly.adm2,admin2.dir),aes(fill=logit.mean))
v2 <- clean_map_theme + geom_sf(data = merge(poly.adm2,admin2.dir),aes(fill=log(var.logit.est)))

print(ggarrange(plotlist = list(m1,v1,m2,v2)))

# scatter plots with covariates -------

admin2.dir <- merge(admin2.dir,admin.key)
admin2.dir$full_name <- paste0(admin2.dir$NAME_1,'_',admin2.dir$NAME_2)
admin2.dir <- merge(admin2.dir,cmat[,c('elev_log','precip','temp','ntlights_log','tthc_log','malaria','LGA_correct')],
                    by.x = 'full_name',by.y = 'LGA_correct')

covariates <- c('elev_log','precip','temp','ntlights_log','tthc_log','malaria')

plot_list <- lapply(covariates,function(cov){
  admin2.dir$cov <- admin2.dir[,c(cov)]
  p <- admin2.dir %>% ggplot(aes((logit.mean), cov)) +
    geom_point(alpha = 0.5) + theme_minimal() +
    geom_smooth(method = "loess", span = 0.25, se = FALSE, color = "blue") +
    labs(title = cov)
})

annotate_figure(ggarrange(plotlist = plot_list),
                top = text_grob('Admin2 weighted estimates v covariates',face='bold',size=14))

plot_list <- lapply(covariates,function(cov){
  admin2.dir$cov <- admin2.dir[,c(cov)]
  p <- admin2.dir %>% filter(var.logit.est > 1e-10) %>% ggplot(aes(log(var.logit.est), cov)) +
    geom_point(alpha = 0.5) + theme_minimal() +
    geom_smooth(method = "loess", span = 0.25, se = FALSE, color = "blue") +
    labs(title = cov)
})

annotate_figure(ggarrange(plotlist = plot_list),
                top = text_grob('Admin2 design-based variance estimates v covariates',face='bold',size=14))


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
# out <- get_variogram(admin1.mat,admin1.dir$logit.mean)
# plot(out, ylim = c(0,max(out)),
#      xlab = '', ylab = '',
#      main = paste('Admin1 weighted estimates for', indicator))
# 
# out <- get_variogram(admin1.mat,log(admin1.dir$var.logit.est))
# plot(out, ylim = c(0,max(out)),
#      xlab = '', ylab = 'log scale',
#      main = paste('Admin1 variance estimates for', indicator))
# 
# # ADMIN2
# out <- get_variogram(admin2.mat,admin2.dir$logit.mean)
# plot(out, ylim = c(0,max(out,na.rm=T)),
#      xlab = '', ylab = '',
#      main = paste('Admin2 weighted estimates for', indicator))
# 
# var.censored <- admin2.dir$var.logit.est
# var.censored[var.censored < 1e-10] <- NA
# out <- get_variogram(admin2.mat,log(var.censored))
# plot(out, ylim = c(0,max(out,na.rm=T)),
#      xlab = '', ylab = 'log scale',
#      main = paste('Admin2 variance estimates for', indicator))

#}
#dev.off()


# get variograms for IID FH ------
hyperpc.iid = list(prec = list(prior = "pc.prec", param = c(5, 0.01)))

fh1.iid.fit <- inla(logit.mean ~ f(admin1,model='iid',hyper = hyperpc.iid),
                    family = 'gaussian',
                    data = admin1.dir,
                    scale = 1/admin1.dir$var.logit.est,
                    control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))), 
                    control.compute = list(config = TRUE))

out <- get_variogram(admin1.mat,fh1.iid.fit$summary.linear.predictor$`0.5quant`)
plot(out, ylim = c(0,max(out)),
     xlab = '', ylab = '',
     main = 'Variogram for admin1 smoothed estimates')

fh2.iid.fit <- inla(logit.mean ~ f(admin2,model='iid',hyper = hyperpc.iid),
                    family = 'gaussian',
                    data = admin2.dir.stable,
                    scale = 1/admin2.dir.stable$var.logit.est,
                    control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))), 
                    control.compute = list(config = TRUE))

# variogram goes flat because of shrinkage?
out <- get_variogram(admin2.mat,fh2.iid.fit$summary.linear.predictor$`0.5quant`)
plot(out, ylim = c(0,max(out,na.rm=T)),
     xlab = '', ylab = '',
     main = 'Variogram for admin1 smoothed estimates')

# covariate models ------------

admin2.dir <- merge(admin2.dir,admin.key)
admin2.dir$full_name <- paste0(admin2.dir$NAME_1,'_',admin2.dir$NAME_2)
admin2.dir <- merge(admin2.dir,cmat[,c('elev_log','precip','temp','ntlights_log','tthc_log','malaria','LGA_correct')],
      by.x = 'full_name',by.y = 'LGA_correct')

admin2.dir.stable <- admin2.dir[!is.na(admin2.dir$logit.mean),]
admin2.dir.stable <- admin2.dir.stable[admin2.dir.stable$var.logit.est > 1e-10,]

covariates <- c('elev_log','precip','temp','ntlights_log','tthc_log','malaria')
cov2.fit <- inla(as.formula(paste('logit.mean ~',paste(covariates,collapse = '+'))),
                 family = 'gaussian',
                 data = admin2.dir.stable,
                 scale = 1/admin2.dir.stable$var.logit.est,
                 control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))), 
                 control.compute = list(config = TRUE))
summary(cov2.fit)
samples <- get_inla_samples_local(cov2.fit,1000)

# errors for some reason
cbind(1,admin2.dir[,covariates]) %*% t(samples)


# once I have this, make variogram and compare to raw direct estimates

