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

load(file=paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Kenya_Example/KEN_2022_cont_outcomes.rda'))

geo <- getDHSgeo(country='Kenya',year=2022)
cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, by.adm1 = "NAME_1",by.adm2 = "NAME_2")
# cluster.info$data <- cluster.info$data %>%
#   rename(NAME_1 = admin1.name, NAME_2 = admin2.name) %>%
#   dplyr::select(-c(LONGNUM,LATNUM,geometry,admin2.name.full))
# cluster.info$data <- merge(cluster.info$data,admin.key,by=c('NAME_1', 'NAME_2'))

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

cluster.info$data <- merge(cluster.info$data, dat %>% dplyr::select(cluster,v025) %>% unique())

clean_map_theme + geom_sf(data=poly.adm1,lwd=0.5,fill='transparent') +  geom_sf(data=poly.adm2,color='grey60',fill='transparent') +
  xlab('') + ylab('') +
  geom_point(data=cluster.info$data,aes(LONGNUM,LATNUM,col=factor(v025)),size=0.25,pch=3) +
  scale_color_manual(name='',labels=c('1'='Urban','2'='Rural'),values=c('1'='red','2'='blue')) +
  theme(legend.position = 'bottom',legend.text = element_text(size=12),panel.border = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)))

# load covariates (country specific) ----------------
setwd('/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Kenya_Example')

ur_weights_adm1 <- readRDS('Kenya_2022_admin1_u5_ur_weights.rds')
ur_weights_adm2 <- readRDS('Kenya_2022_admin2_u5_ur_weights.rds')

load(file='KEN_Covariates/Kenya_admin2_covariates.rda')
load(file='KEN_Covariates/Kenya_admin1_covariates.rda')

cmat_admin1 <- merge(cmat_admin1,ur_weights_adm1)
cmat_admin2 <- merge(cmat_admin2,ur_weights_adm2)

mean_covariates <- c('urb_frac','density_log', 'nt_lights_log', 'tthc_log', 'precip', 'temp', 'elev')
var_covariates <- c('urb_frac', 'pop_var_log', 'nt_lights_var_log', 'tthc_var_log', 'precip_var_log', 'temp_var_log', 'elev_var_log', 'area_log')

mean_model_matrix <- cbind(rep(1,n_admin2),cmat_admin2[,mean_covariates])
var_model_matrix <- cbind(rep(1,n_admin2),cmat_admin2[,var_covariates])


# specify outcome and get direct estimates ----------
var_t <- 'haz'

dir.dat <- dat[!is.na(dat[,c(var_t)]),]
dir.dat$value <- dir.dat[,c(var_t)]
dir.dat$wt <- dir.dat$v005/1e6
options(survey.lonely.psu = "adjust")
my.svydesign <- survey::svydesign(ids = stats::formula("~cluster"), 
                                  strata = ~v023, nest = T, weights = ~wt, data = dir.dat)

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

# setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance")
# save(admin1.dir,file = paste0(country,'_',var_t,'_admin1_weighted_estimates.rda'))
# save(admin2.dir,file = paste0(country,'_',var_t,'_admin2_weighted_estimates.rda'))


# standard FH ----------
 
hyperpc.iid = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
hyperpc.bym2 = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                    phi = list(prior = "pc", param = c(0.5, 0.5)))
 
 fh2.bym2.cov <- inla(mean ~ f(admin2,model='bym2',hyper=hyperpc.bym2,graph = admin2.mat,scale.model = T)+ urb_frac + 
                    density_log+ nt_lights_log+ tthc_log+ 
                      precip+ temp+ elev,
             family = "gaussian",
             quantiles = c(0.05,0.95),
             data = admin2.dir.stable,
             scale = 1/admin2.dir.stable$variance,
             control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
             control.compute = list(config = TRUE))
 
 std.fit <- fh2.bym2.cov$summary.fitted.values
 std.fit$admin2 <- 1:n_admin2
 
 
# Naive approximation -------
 mod <- cmdstan_model(stan_file = "/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing/Stan/NonStrat.stan")
 
 data_areas <- which(!is.na(admin2.dir.stable$mean))
 df_naive <- v_hat_scaled_naive <- rep(NA,n_admin2)
 
 for(area in data_areas){
   tmp_D <- dir.dat %>% filter(admin2 ==area) %>% group_by(v025,cluster) %>% summarise(n())
   N_D <- c(sum(tmp_D$v025==2),sum(tmp_D$v025==1))
   
   df_naive[area] <- pmax(1,sum(N_D) - sum(N_D>0))
   v_hat_scaled_naive[area] <- admin2.dir.stable$variance[area]*pmax(1,sum(N_D) - sum(N_D>0)) # naive scaling
 }
 
 data_list_naive = list(m=n_admin2,
                        m_data=length(data_areas),
                        data_areas = data_areas,
                        y=admin2.dir.stable$mean[data_areas],
                        v_hat_scaled = v_hat_scaled_naive[data_areas],
                        df=df_naive[data_areas],
                        p_mean = ncol(mean_model_matrix),
                        X = mean_model_matrix,
                        bym2_mean = 1,
                        p_var = ncol(var_model_matrix),
                        Z = var_model_matrix,
                        bym2_var = 1,
                        # all bym2 stuff
                        N_edges = nrow(nodes2),
                        node_1 = nodes2$node1,
                        node_2 = nodes2$node2,
                        car_scale = Q2_scaled[1,1]/Q.admin2[1,1]) 
 
 fit1 <- mod$sample(
   data = data_list_naive,
   chains = 4,
   parallel_chains = 4,
   iter_warmup = 1000,
   iter_sampling = 1000,
  # show_messages = F,
   #show_exceptions = F,
   refresh=200)
 
 theta_fit <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
 
 naive.fit <- data.frame(admin2 = 1:n_admin2,
            mean = apply(theta_fit,2,mean),
            upper = apply(theta_fit,2,quantile,probs=0.95),
            lower = apply(theta_fit,2,quantile,probs=0.05))
 
 # compare models ------
 par(mfrow=c(1,1))
 
 plot((std.fit$mean), (naive.fit$mean),
      xlab = 'No variance smoothing',ylab='Naive variance smoothing')
 abline(0,1)
 
 plot((std.fit$`0.95quant` - std.fit$`0.05quant`), (naive.fit$upper - naive.fit$lower),
      xlab = 'No variance smoothing',ylab='Naive variance smoothing')
 abline(0,1)
 
 merge(poly.adm2,naive.fit) %>% ggplot() + geom_sf(aes(fill=mean))
 merge(poly.adm2,std.fit) %>% ggplot() + geom_sf(aes(fill=mean))
 
 range_mean <- range(naive.fit$mean,std.fit$mean)
 m1 <- clean_map_theme + geom_sf(data = merge(poly.adm2,std.fit) ,aes(fill=mean)) + 
   scale_fill_continuous(lim = range_mean) +
   ggtitle('No variance smoothing')
 
 m2 <- clean_map_theme + geom_sf(data = merge(poly.adm2,naive.fit) ,aes(fill=mean)) + 
   scale_fill_continuous(lim = range_mean) +
   ggtitle('Variance smoothing (naive approximation)')
 
 ggarrange(plotlist = list(m1,m2),common.legend=T)
 
 
 range_width <- range(naive.fit$upper - naive.fit$lower, std.fit$`0.95quant` - std.fit$`0.05quant`)
 m1 <- clean_map_theme + geom_sf(data = merge(poly.adm2,std.fit) ,aes(fill=`0.95quant` - `0.05quant`)) + 
   scale_fill_continuous(name = 'Credible interval width',lim = range_width) +
   ggtitle('No variance smoothing')
 
 m2 <- clean_map_theme + geom_sf(data = merge(poly.adm2,naive.fit) ,aes(fill=upper - lower)) + 
   scale_fill_continuous(name = 'Credible interval width',lim = range_width) +
   ggtitle('Variance smoothing (naive approximation)')
 
 ggarrange(plotlist = list(m1,m2),common.legend=T)
 
## properties of the data ------

# are variances different between areas -- yes!
dir.dat %>% group_by(admin1) %>% summarise(v = var(value)/n()) %>% ggplot() +geom_histogram(aes(v))

# is variance different across strata -- yes, within area
dir.dat %>% group_by(v025) %>% summarise(sd = var(value)) 

t <- dir.dat %>% group_by(admin2,v025) %>% summarise(v = var(haz)) %>% pivot_wider(values_from = v, names_from = v025) %>% summarise(val = `1`/`2`)
t2 <- dir.dat %>% group_by(admin2,v025) %>% summarise(m = mean(haz)) %>% pivot_wider(values_from = m, names_from = v025) %>% summarise(val = `1`-`2`)

plot(t2$val,t$val)


dir.dat %>% group_by(admin2) %>% summarise(v = log(var(haz))) %>% summarise(var(v))
dir.dat %>% group_by(admin2) %>% summarise(v = log(var(haz))) %>% summarise(mean(v))

dir.dat %>% group_by(admin2) %>% summarise(v = mean(haz)) %>% summarise(var(v))

# normality of cluster means within strata

tmp <- dir.dat %>% group_by(cluster,admin1,v023,v025) %>% summarise(val = mean(haz),n=n()) %>% 
  mutate(t=sqrt(n)*val) %>% ggplot(aes(sample=t)) + geom_qq() + geom_qq_line() 

tmp <- dir.dat %>% ggplot(aes(sample=haz)) + geom_qq() + geom_qq_line() 


#pdf('Normality of cluster means.pdf')
tmp + ggforce::facet_wrap_paginate(~v023,nrow=5,ncol=5,page=1)
tmp + ggforce::facet_wrap_paginate(~v023,nrow=5,ncol=5,page=2)
tmp + ggforce::facet_wrap_paginate(~v023,nrow=5,ncol=5,page=3)
tmp + ggforce::facet_wrap_paginate(~v023,nrow=5,ncol=5,page=4)


tmp + ggforce::facet_wrap_paginate(~admin1,nrow=5,ncol=5,page=1)
tmp + ggforce::facet_wrap_paginate(~admin1,nrow=5,ncol=5,page=2)

out <- rep(NA,n_admin2)
for(area in 1:n_admin2){
  t <- (dir.dat[dir.dat$admin2==area,]$haz - mean(dir.dat[dir.dat$admin2==area,]$haz))/sd(dir.dat[dir.dat$admin2==area,]$haz)
  if(length(t)>0){
    out[area] <- ks.test(t, "pnorm", mean = 0, sd = 1)$p.value
  }
}

### testing -------------





