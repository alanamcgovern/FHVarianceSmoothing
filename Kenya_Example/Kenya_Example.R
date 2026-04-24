library(rdhs)
library(haven)
library(surveyPrev)
library(sf)
library(spdep)
library(tidyverse)
library(INLA)
#library(igraph)
library(ggpubr)
library(ggpattern)
library(cmdstanr)
library(expm)
library(patchwork)
library(scico)
source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")

country <- 'Kenya'
country.abbrev <- 'KEN'
dhs.abbrev <- 'KE'
survey_year <- 2022

# POLYGONS -----
setwd('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm')
poly.adm2 <- st_read(dsn = paste0('gadm41_',country.abbrev,'_shp'), layer = paste0("gadm41_",country.abbrev,"_2"), options = "ENCODING=UTF-8")

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


# DATA AND COVARIATES -------
load(file=paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Kenya_Example/KEN_2022_cont_outcomes.rda'))

setwd('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Kenya_Example')

ur_weights_adm1 <- readRDS('Kenya_2022_admin1_u5_ur_weights.rds')
ur_weights_adm2 <- readRDS('Kenya_2022_admin2_u5_ur_weights.rds')

load(file='KEN_Covariates/Kenya_admin2_covariates.rda')
load(file='KEN_Covariates/Kenya_admin1_covariates.rda')

cmat_admin1 <- merge(cmat_admin1,ur_weights_adm1)
cmat_admin2 <- merge(cmat_admin2,ur_weights_adm2)

mean_covariates <- c('urb_frac','density_log', 'nt_lights_log', 'tthc_log', 'precip', 'temp', 'elev')
var_covariates <- c('urb_frac', #'pop_var_log', 
                    'nt_lights_var_log', 'tthc_var_log', 'precip_var_log', 'temp_var_log', 'elev_var_log', 'area_log')

mean_model_matrix <- cbind(rep(1,n_admin2),cmat_admin2[,mean_covariates])
var_model_matrix <- cbind(rep(1,n_admin2),cmat_admin2[,var_covariates])

# SPECIFY outcome and get DIRECT estimates ---------------

var_t <- 'haz'

dir.dat <- dat[!is.na(dat[,c(var_t)]),]
dir.dat$value <- dir.dat[,c(var_t)]
dir.dat$urban <- ifelse(dir.dat$v025==2,0,dir.dat$v025)
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
admin2.dir.stable$variance <- ifelse(admin2.dir.stable$variance < 1e-5,NA,admin2.dir.stable$variance)
admin2.dir.stable$mean <- ifelse(is.na(admin2.dir.stable$variance),NA,admin2.dir.stable$mean)

# COVARIATE PLOTS ---------
for(cov in mean_covariates){
  plot(admin2.dir.stable$mean,admin2.dir.stable[,cov],main=cov)
  print(cor(admin2.dir.stable$mean,admin2.dir.stable[,cov],use='pairwise.complete.obs'))
}
plot(admin2.dir.stable$mean,admin2.dir.stable[,mean_covariates[2]])

c4 <- clean_map_theme + geom_sf(data=merge(poly.adm2,admin2.dir),aes(fill=urb_frac),color='grey20',lwd=0.25) + 
  ggtitle('Urban population proportion') + 
  scale_fill_viridis_c(name = '') + theme(legend.position = 'bottom',title = element_text(size=8))

# standardized, log-transformed population density of each area
c5 <- clean_map_theme + geom_sf(data=merge(poly.adm2,admin2.dir),aes(fill=density_log),color='grey20',lwd=0.25) + 
  ggtitle('Population density') + 
  scale_fill_viridis_c(name = '') + theme(legend.position = 'bottom',title = element_text(size=8))


c1 <- clean_map_theme + geom_sf(data=merge(poly.adm2,admin2.dir),aes(fill=temp),color='grey20',lwd=0.25) + 
  ggtitle('Temperature') + 
  scale_fill_viridis_c(name = '') + theme(legend.position = 'bottom',title = element_text(size=8))

# average across area, standardized
c2 <- clean_map_theme + geom_sf(data=merge(poly.adm2,admin2.dir),aes(fill=precip),color='grey20',lwd=0.25) + 
  ggtitle('Precipitation') + 
  scale_fill_viridis_c(name = '') + theme(legend.position = 'bottom',title = element_text(size=8))

# average across area, standardized
c3 <- clean_map_theme + geom_sf(data=merge(poly.adm2,admin2.dir),aes(fill=elev),color='grey20',lwd=0.25) + 
  ggtitle('Elevation') + 
  scale_fill_viridis_c(name = '') + theme(legend.position = 'bottom',title = element_text(size=8))

# average across area, log-transformed and then standardized
c6 <- clean_map_theme + geom_sf(data=merge(poly.adm2,admin2.dir),aes(fill=nt_lights_log),color='grey20',lwd=0.25) + 
  ggtitle('Nighttime lights') + 
  scale_fill_viridis_c(name = '') + theme(legend.position = 'bottom',title = element_text(size=8))

# average across area, log-transformed and then standardized
c7 <- clean_map_theme + geom_sf(data=merge(poly.adm2,admin2.dir),aes(fill=tthc_log),color='grey20',lwd=0.25) + 
  ggtitle('Motorized travel time \n to healthcare') + 
  scale_fill_viridis_c(name = '') + theme(legend.position = 'bottom',title = element_text(size=8))

# 680 x 890
ggarrange(plotlist = list(c4,c5,c1,c2,c3,c6,c7),nrow=3,ncol=3)



# STANDARD FH for comparison ------------
# mod_std <- cmdstan_model(stan_file = "/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Stan/Standard.stan")
# 
# data_areas <- which(!is.na(admin2.dir.stable$mean))
# data_list = list(m=n_admin2,
#                  m_data=length(data_areas),
#                  data_areas = data_areas,
#                  y=admin2.dir.stable$mean[data_areas],
#                  v_hat = admin2.dir.stable$variance[data_areas],
#                  p_mean = ncol(mean_model_matrix),
#                  X = mean_model_matrix,
#                  bym2_mean = 1,
#                  # all bym2 stuff
#                  N_edges = nrow(nodes2),
#                  node_1 = nodes2$node1,
#                  node_2 = nodes2$node2,
#                  car_scale = Q2_scaled[1,1]/Q.admin2[1,1]) 
# 
# fit0 <- mod_std$sample(
#   data = data_list,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1000,
#   iter_sampling = 1000,
#   adapt_delta = 0.99,
#   # show_messages = F,
#   #  show_exceptions = F,
#   refresh=200)


hyperpc.bym2 = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                    phi = list(prior = "pc", param = c(0.5, 0.5)))

formula_str <- paste0("mean ~ ", paste(mean_covariates, collapse = " + "), 
                      " + f(admin2, model='bym2', graph=admin2.mat, hyper=hyperpc.bym2, scale.model=T, adjust.for.con.comp=T)")

start.time <- Sys.time()
std.fit <- inla(as.formula(formula_str),
  family = "gaussian",
  data = admin2.dir.stable,
  scale = 1 / variance,
  control.fixed = list(quantiles=c(0.05,0.5,0.95)),
  control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
  control.compute = list(config = TRUE)
)
Sys.time() - start.time

std.samples <- get_inla_samples_local(std.fit, 1000)
X <- model.matrix(as.formula(paste0("~ ", paste(mean_covariates, collapse = " + "))), data = admin2.dir.stable)
n_fixed <- ncol(X)
intercept_col <- 2*n_admin2 + 1
fixed_cols <- (intercept_col):(intercept_col + (n_fixed - 1))
eta.samples <- std.samples[, fixed_cols]%*% t(X[, , drop = FALSE]) + std.samples[, 1:n_admin2]

std.res <- data.frame(mean  = apply(eta.samples,2,mean),
           sd = apply(eta.samples,2,sd),
           lower = apply(eta.samples,2,quantile,prob=0.05),
           upper = apply(eta.samples,2,quantile,prob=0.95),
           area = 1:n_admin2,
           model='std')

std.res$stunt <- apply(eta.samples,2,function(x){mean(x < -2)})
std.res$lessneg1 <- apply(eta.samples,2,function(x){mean(x < -1)})
std.res$less0 <- apply(eta.samples,2,function(x){mean(x < 0)})

ranks <- t(apply(eta.samples, 1, rank, ties.method = "average"))
prop_greater <- (ranks - 1) / n_admin2
std.res$lower.1 <- colMeans(prop_greater < 0.1)
std.res$lower.25 <- colMeans(prop_greater < 0.25)


# CHOOSE covariates for variance model ---------
formula_str <- paste0("log(variance) ~ ", paste(var_covariates, collapse = " + "), 
                      " + f(admin2, model='bym2', graph=admin2.mat, hyper=hyperpc.bym2, scale.model=T, adjust.for.con.comp=T)")

test.fit <- inla(as.formula(formula_str),
                family = "gaussian",
                data = admin2.dir.stable)

# NAIVE distribution models --------------
mod_naive <- cmdstan_model(stan_file = "/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Stan/NonStrat.stan")

data_areas <- which(!is.na(admin2.dir.stable$mean))
Cons_naive <- scale_naive <- df_naive <- v_hat_scaled_naive <- rep(NA,n_admin2)

for(area in data_areas){
  tmp_D <- dir.dat %>% filter(admin2==area) %>% group_by(urban,cluster) %>% summarise(n=n(),wt = unique(v005)/1e6)
  N_D <- c(sum(1-tmp_D$urban),sum(tmp_D$urban))
  
  df_naive[area] <- pmax(1,sum(N_D) - sum(N_D>0))
  v_hat_scaled_naive[area] <- admin2.dir.stable$variance[area]*pmax(1,sum(N_D) - sum(N_D>0))
  Cons_naive[area] <- 1/(sum(N_D)*mean(tmp_D$n))
  scale_naive[area] <- 1/(sum(N_D)*mean(tmp_D$n)*pmax(1,sum(N_D) - sum(N_D>0)))
}

data_list = list(m=n_admin2,
                 m_data=length(data_areas),
                 data_areas = data_areas,
                 y=admin2.dir.stable$mean[data_areas],
                 v_hat_scaled = v_hat_scaled_naive[data_areas],
                 df=df_naive[data_areas],
                 Cons = Cons_naive[data_areas],
                 p_mean = ncol(mean_model_matrix),
                 X = mean_model_matrix,
                 bym2_mean = 1,
                 # all bym2 stuff
                 N_edges = nrow(nodes2),
                 node_1 = nodes2$node1,
                 node_2 = nodes2$node2,
                 car_scale = Q2_scaled[1,1]/Q.admin2[1,1]) 

### with BYM2 and covariates ------------
data_list$bym2_var <- 1
data_list$p_var = ncol(var_model_matrix)
data_list$Z = var_model_matrix

# 38.6 seconds
fit1 <- mod_naive$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
 # show_messages = F,
#  show_exceptions = F,
  refresh=200)

theta_fit <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
mod1_res <- data.frame(area = 1:n_admin2,
                        model = 'Mod1',
                        mean = apply(theta_fit,2,mean),
                        upper = apply(theta_fit,2,quantile,probs=0.95),
                        lower = apply(theta_fit,2,quantile,probs=0.05))

ranks <- t(apply(theta_fit, 1, rank, ties.method = "average"))
prop_greater <- (ranks - 1) / n_admin2
mod1_res$lower.1 <- colMeans(prop_greater < 0.1)
mod1_res$lower.25 <- colMeans(prop_greater < 0.25)
mod1_res$stunt <- apply(theta_fit,2,function(x){mean(x < -2)})

## diagnostics
summary_df <- fit1$summary()
diag_df <- fit1$diagnostic_summary()

mod1_diag <- data.frame(model = 'Mod1',
                         divergences = sum(diag_df$num_divergent),
                         max_treedepth = sum(diag_df$num_max_treedepth),
                         min_efbmi = min(diag_df$ebfmi),
                         max_rhat = max(summary_df$rhat,na.rm=T),
                         bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                         tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail)
rm(fit1)
rm(theta_fit)

### with BYM2 and no covariates ------------
data_list$bym2_var <- 1
data_list$p_var = 1
data_list$Z = cbind(rep(1,n_admin2))

# 29 seconds
fit1b <- mod_naive$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
  # show_messages = F,
  #  show_exceptions = F,
  refresh=200)

theta_fit <- as.matrix(unclass(fit1b$draws(variables = c("theta"), format = "matrix")))
mod1b_res <- data.frame(area = 1:n_admin2,
                       model = 'Mod1b',
                       mean = apply(theta_fit,2,mean),
                       sd = apply(theta_fit,2,sd),
                       upper = apply(theta_fit,2,quantile,probs=0.95),
                       lower = apply(theta_fit,2,quantile,probs=0.05))

ranks <- t(apply(theta_fit, 1, rank, ties.method = "average"))
prop_greater <- (ranks - 1) / n_admin2
mod1b_res$lower.1 <- colMeans(prop_greater < 0.1)
mod1b_res$lower.25 <- colMeans(prop_greater < 0.25)

mod1b_res$stunt <- apply(theta_fit,2,function(x){mean(x < -2)})
mod1b_res$less0 <- apply(theta_fit,2,function(x){mean(x < 0)})
mod1b_res$lessneg1 <- apply(theta_fit,2,function(x){mean(x < -1)})

## diagnostics
summary_df <- fit1b$summary()
diag_df <- fit1b$diagnostic_summary()

mod1b_diag <- data.frame(model = 'Mod1b',
                        divergences = sum(diag_df$num_divergent),
                        max_treedepth = sum(diag_df$num_max_treedepth),
                        min_efbmi = min(diag_df$ebfmi),
                        max_rhat = max(summary_df$rhat,na.rm=T),
                        bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                        tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail)

rm(fit1b)
rm(theta_fit)

### with IID and no covariates (except urban prop) -------

data_list$bym2_var <- 0
data_list$p_var = 2
data_list$Z = cbind(rep(1,n_admin2),var_model_matrix[,c('urb_frac')])

# 22.4 s
fit2 <- mod_naive$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
 # show_messages = F,
#  show_exceptions = F,
  refresh=200)


theta_fit <- as.matrix(unclass(fit2$draws(variables = c("theta"), format = "matrix")))
mod1a_res <- data.frame(area = 1:n_admin2,
                     model = 'Mod1a',
                     mean = apply(theta_fit,2,mean),
                     upper = apply(theta_fit,2,quantile,probs=0.95),
                     lower = apply(theta_fit,2,quantile,probs=0.05))

## diagnostics
summary_df <- fit2$summary()
diag_df <- fit2$diagnostic_summary()

mod1a_diag <- data.frame(model = 'Mod1a',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail)

rm(fit2)
rm(theta_fit)

# SASW distribution models ---------------------------------
mod_sasw <- cmdstan_model(stan_file = "/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Stan/Strat_noKappa.stan")

data_areas <- which(!is.na(admin2.dir.stable$mean))
cons_exact <- v_hat_scaled_exact <- q_start <- rep(NA,n_admin2)
q <- q_id <- nu <- NULL
tol <- 1e-10
for(area in data_areas){
  tmp <- dir.dat %>% filter(admin1 == admin.key[admin.key$admin2 == area,]$admin1) %>% group_by(urban,cluster,admin2) %>% summarise(n=n(),wt = unique(v005)/1e6)
  tmp$wt <- pmin(tmp$wt,quantile(tmp$wt,0.95))
  
  tmp[tmp$admin2 != area,]$wt <- 0
  N <- c(sum(1-tmp$urban),sum(tmp$urban))
  
  omega <- tmp$n*tmp$wt
  
  B_comps <- lapply(1:2,function(h){
    if(N[h]>0){
      return(N[h]/(N[h]-1)*(diag(1,N[h]) - 1/N[h]*matrix(1,N[h],N[h])))
    }
  })
  
  which.not.null <- c(1:2)[!sapply(B_comps, is.null)]
  B <- bdiag(B_comps[which.not.null])
  W <- diag(omega)%*%(diag(1,sum(N)) - (rep(1,sum(N))%*%t(omega))/sum(omega))
  C <- t(W)%*%B%*%W
  S <- diag(1/tmp$n)
  
  L <- chol(S)              # S = L' L
  Li <- backsolve(L, diag(nrow(S)))  # L^{-1}
  
  # Form scale-free matrix A = S^{1/2} C S^{1/2}
  A <- t(L) %*% C %*% L
  
  out <- eigen(A, symmetric = TRUE)
  V <- Li %*% out$vectors
  
  keep_id <- which(out$values > tol)
  q_tmp <- out$values[keep_id]
  q_start[area] <- length(q) + 1
  q_id <- c(q_id,rep(area,length(keep_id)))
  q <- c(q,q_tmp)
  
  V <- V[,keep_id]
  nu <- c(nu,t(V)%*%tmp$urban)
  
  v_hat_scaled_exact[area] <- admin2.dir.stable$variance[area]*(t(omega)%*%S%*%omega)
  cons_exact[area] <- (t(omega)%*%S%*%omega)/sum(omega)^2

}

data_list = list(m=n_admin2,
                 m_data=length(data_areas),
                 data_areas = data_areas,
                 y=admin2.dir.stable$mean[data_areas],
                 v_hat_scaled = v_hat_scaled_exact[data_areas],
                 Cons = cons_exact[data_areas],
                 lenq = length(q),
                 q_id = q_id,
                 q_start = q_start[data_areas],
                 q_per_area = as.vector(table(q_id)),
                 q = q, nu = nu,
                 p_mean = ncol(mean_model_matrix),
                 X = mean_model_matrix,
                 bym2_mean = 1,
                 # all bym2 stuff
                 N_edges = nrow(nodes2),
                 node_1 = nodes2$node1,
                 node_2 = nodes2$node2,
                 car_scale = Q2_scaled[1,1]/Q.admin2[1,1],
                 bias_adj = 0) 

### ASSESS SATT APPROX -----------

mu_diff <- c(1,2)[1]
sig2 <- c(1,4)[1]

delta <- (nu*mu_diff)^2/sig2

# 830 x 400
set.seed(260216)
par(mfrow=c(2,4),
    mar = c(1, 1, 1, 1))
for(a in sort(sample(data_areas,8))){
  q_tmp <- q[q_id==a]
  delta_tmp <- delta[q_id==a]
  r <- length(q_tmp)
  
  dist_exact <- replicate(1e4,sum(q_tmp*sapply(1:r,function(j){rchisq(1,1,delta_tmp[j])})))
  q_exact <- quantile(dist_exact,seq(0.01,0.99,0.01))
  
  Q1 <- sum(q_tmp*(1+delta_tmp))
  Q2 <- sum(q_tmp^2*(1+2*delta_tmp))
  dist_satt <- Q2/Q1*rchisq(1e4,Q1^2/Q2)
  q_satt <- quantile(dist_satt,seq(0.01,0.99,0.01))
  
  plot(q_exact,q_satt,main = paste0('Area ',a),axes=F,cex=0.5)
  box()
  abline(0,1)
}

### fit with BYM2 and covariates -----------

data_list$p_var = ncol(var_model_matrix)
data_list$Z = var_model_matrix
data_list$bym2_var = 1

# 109.2 s
fit3 <- mod_sasw$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
 # show_messages = F,
#  show_exceptions = F,
  refresh=100)

theta_fit <- as.matrix(unclass(fit3$draws(variables = c("theta"), format = "matrix")))
mod2_res <- data.frame(area = 1:n_admin2,
                     model = 'Mod2',
                     mean = apply(theta_fit,2,mean),
                     upper = apply(theta_fit,2,quantile,probs=0.95),
                     lower = apply(theta_fit,2,quantile,probs=0.05))

ranks <- t(apply(theta_fit, 1, rank, ties.method = "average"))
prop_greater <- (ranks - 1) / n_admin2
mod2_res$lower.1 <- colMeans(prop_greater < 0.1)
mod2_res$lower.25 <- colMeans(prop_greater < 0.25)
mod2_res$stunt <- apply(theta_fit,2,function(x){mean(x < -2)})

summary_df <- fit3$summary()
diag_df    <- fit3$diagnostic_summary()

mod2_diag <- data.frame(model = 'Mod2',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail)

rm(theta_fit)
rm(fit3)

### fit with BYM2 and no covariates -----------

data_list$p_var = 1
data_list$Z = cbind(rep(1,n_admin2))
data_list$bym2_var = 1

# 69 s
fit3 <- mod_sasw$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
  # show_messages = F,
  #  show_exceptions = F,
  refresh=100)

theta_fit <- as.matrix(unclass(fit3$draws(variables = c("theta"), format = "matrix")))
mod2b_res <- data.frame(area = 1:n_admin2,
                       model = 'Mod2b',
                       mean = apply(theta_fit,2,mean),
                       sd = apply(theta_fit,2,sd),
                       upper = apply(theta_fit,2,quantile,probs=0.95),
                       lower = apply(theta_fit,2,quantile,probs=0.05))

ranks <- t(apply(theta_fit, 1, rank, ties.method = "average"))
prop_greater <- (ranks - 1) / n_admin2
mod2b_res$lower.1 <- colMeans(prop_greater < 0.1)
mod2b_res$lower.25 <- colMeans(prop_greater < 0.25)

mod2b_res$stunt <- apply(theta_fit,2,function(x){mean(x < -2)})
mod2b_res$less0 <- apply(theta_fit,2,function(x){mean(x < 0)})
mod2b_res$lessneg1 <- apply(theta_fit,2,function(x){mean(x < -1)})

summary_df <- fit3$summary()
diag_df    <- fit3$diagnostic_summary()

mod2b_diag <- data.frame(model = 'Mod2b',
                        divergences = sum(diag_df$num_divergent),
                        max_treedepth = sum(diag_df$num_max_treedepth),
                        min_efbmi = min(diag_df$ebfmi),
                        max_rhat = max(summary_df$rhat,na.rm=T),
                        bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                        tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail)
rm(fit3)
rm(theta_fit)

### fit with IID and no covariates (except urban prop) -----------

data_list$p_var = 2
data_list$Z = cbind(rep(1,n_admin2),var_model_matrix[,c('urb_frac')])
data_list$bym2_var = 0

# 63.5 s
fit4 <- mod_sasw$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
 # show_messages = F,
#  show_exceptions = F,
  refresh=100)

theta_fit <- as.matrix(unclass(fit4$draws(variables = c("theta"), format = "matrix")))
mod2a_res <- data.frame(area = 1:n_admin2,
                     model = 'Mod2a',
                     mean = apply(theta_fit,2,mean),
                     upper = apply(theta_fit,2,quantile,probs=0.95),
                     lower = apply(theta_fit,2,quantile,probs=0.05))

summary_df <- fit4$summary()
diag_df <- fit4$diagnostic_summary()

mod2a_diag <- data.frame(model = 'Mod2a',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail)

rm(fit4)

# COMBINE AND SAVE MODELS ----------

mods_res <- rbind(std.res[,c('area','model','mean','sd','upper','lower','lower.1','lower.25','stunt','less0','lessneg1')],
                  mod1b_res,
                  #mod1a_res,
                  mod2b_res)
                 # mod2a_res)
mods_diag <- rbind(mod1b_diag,
                   #mod1a_diag,
                   mod2b_diag)
                  # mod2a_diag)

setwd("/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Kenya_Example")
save(mods_res,file='model_results.rda')
save(mods_diag,file='model_diagnostics.rda')

load(file='model_results.rda')
load(file='model_diagnostics.rda')
mods_res <- mods_res %>% mutate(ci.width = upper - lower) %>% rename(admin2 = area)

# PLOTS ---------

model_labs <- c(
  Mod1b = "Simple",
  Mod2b = "SASW",
  std = "Standard")

### assessing normality assumption -------
hist(dir.dat$value)
dir.dat %>% ggplot() + geom_histogram(aes(value)) + facet_grid(~urban)
dir.dat %>% ggplot() + geom_histogram(aes(value)) + facet_wrap(~v023)
dir.dat %>% ggplot() + geom_histogram(aes(value)) + facet_wrap(~admin1)

### distribution of cluster sizes ---------
cluster_sizes <- dir.dat %>% group_by(admin1,admin2,v025,cluster) %>% summarise(n= n())

par(mfrow=c(1,2),mar = c(5, 1, 4, 1),   # bottom, left, top, right
    oma = c(0, 0, 0, 0))
hist(cluster_sizes[cluster_sizes$v025==1,]$n,prob=T,
     main = 'Urban clusters',ylab='',yaxt='n',xlab='sampled individuals')
lines(1:50-0.5,dnbinom(1:50,size=8,mu=9),col='red',lwd=2)

hist(cluster_sizes[cluster_sizes$v025==2,]$n,prob=T,
     main = 'Rural clusters',ylab='',yaxt='n',xlab='sampled individuals')
lines(1:50-0.5,dnbinom(1:50,size=4,mu=11),col='red',lwd=2)


### sampled clusters map ---------
cluster_sizes <- dir.dat %>% group_by(admin1,admin2,v025,cluster) %>% summarise(n= n())

admin1_sizes <- cluster_sizes %>% group_by(admin1) %>% summarise(m=n())
admin2_sizes <- cluster_sizes %>% group_by(admin2) %>% summarise(m=n())

lims = c(0,max(admin1_sizes$m))

d1 <- clean_map_theme + geom_sf(data = merge(poly.adm1,admin1_sizes),aes(fill=m),color='grey20',lwd=0.25) + 
  scale_fill_viridis_c(name='Number of sampled clusters',direction = -1,limits = c(0,max(admin1_sizes$m))) + 
  theme(legend.position = 'bottom') +
  ggtitle('First administrative areas')


adm2_merge <- left_join(poly.adm2,admin2_sizes)
adm2_merge[is.na(adm2_merge$m),]$m <- 0
adm2_merge$low_clusters <- I(adm2_merge$m<5)

d2 <- clean_map_theme + geom_sf(data = adm2_merge,aes(fill=m),color='grey20',lwd=0.25) + 
  geom_sf(data = filter(adm2_merge,m==0),fill='grey70') + 
  geom_sf_pattern(data = filter(adm2_merge,low_clusters==1),
                  pattern_density=0.01,
                  fill=NA,
                  pattern = "stripe",    # stripe, crosshatch, or other patterns
                  pattern_angle = 45,
                  color = "black",       # color of hatch lines
                  pattern_spacing = 0.01) +
  scale_fill_viridis_c(name='Number of sampled clusters',direction = -1, limits = c(0,max(admin1_sizes$m))) + 
  theme(legend.position = 'bottom')+
  ggtitle('Second administrative areas')

# 780 x 580
ggarrange(plotlist = list(d1,d2),nrow=1,common.legend = T)

### design-based estimates and variances maps ----------
poly.adm2$admin2_mean <- admin2.dir.stable$mean[poly.adm2$admin2]
poly.adm1$admin1_mean <- admin1.dir$mean[poly.adm1$admin1]

lims <- range(admin1.dir$mean,admin2.dir.stable$mean,na.rm=T)

m1 <- clean_map_theme + geom_sf(data=poly.adm1,aes(fill=admin1_mean),color='grey20',lwd=0.25) + 
  ggtitle('First administrative areas') + 
  scale_fill_viridis_c(name = 'Design-based mean estimate',limits=lims,direction = -1)


m2 <- clean_map_theme + geom_sf(data=poly.adm2,aes(fill=admin2_mean),color='grey80',lwd=0.01) + 
  geom_sf(fill = "transparent", size=1, color = "grey20", lwd=0.25, data = poly.adm2 %>% group_by(NAME_1) %>% summarise()) +
  ggtitle('Second administrative areas') +
  scale_fill_viridis_c(name = 'Design-based mean estimate',limits=lims,direction = -1)

poly.adm2$v <- admin2.dir.stable$variance[poly.adm2$admin2]
poly.adm1$v <- admin1.dir$variance[poly.adm1$admin1]

lims <- range(sqrt(admin1.dir$variance),sqrt(admin2.dir.stable$variance),na.rm=T)

v1 <- clean_map_theme + geom_sf(data=poly.adm1,aes(fill=sqrt(v)),color='grey20',lwd=0.25) + 
  ggtitle('First administrative areas') + 
  scale_fill_scico(name = 'Design-based standard deviation estimate',limits=lims,direction = -1,
                   palette = 'lajolla',
                   trans='log',breaks=c(0.025,0.1,0.4))


v2 <- clean_map_theme + geom_sf(data=poly.adm2,aes(fill=sqrt(v)),color='grey80',lwd=0.01) + 
  geom_sf(fill = "transparent", size=1, color = "grey20", lwd=0.25, data = poly.adm2 %>% group_by(NAME_1) %>% summarise()) +
  ggtitle('Second administrative areas') +
  scale_fill_scico(name = 'Design-based standard deviation estimate',limits=lims,direction = -1,
                   palette = 'lajolla',
                   trans='log',breaks=c(0.025,0.1,0.4))

top_row <- (m1 | m2) + plot_layout(guides = "collect")
bottom_row <- (v1 | v2) + plot_layout(guides = "collect")

# 755 x 850
top_row / bottom_row &
  theme(legend.position = "bottom")




### model results maps -------

lims <- range(mods_res$mean)

m1 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='std',]),aes(fill=mean)) +
  ggtitle('Standard') +
  scale_fill_viridis_c(name = 'Mean',direction = -1,limits = lims)
m2 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod1b',]),aes(fill=mean)) +
  ggtitle('Simple') +
  scale_fill_viridis_c(name = 'Mean',direction = -1,limits = lims)
m3 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod2b',]),aes(fill=mean)) +
  ggtitle('SASW') +
  scale_fill_viridis_c(name = 'Mean',direction = -1,limits = lims)

row1 <- ggarrange(plotlist = list(m1,m2,m3),common.legend = T,nrow=1)

lims <- range(mods_res$ci.width)

v1 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='std',]),aes(fill=ci.width)) +
  ggtitle('Standard') +
  scale_fill_scico(name = 'Width of 90% credible intervals',palette = 'lajolla',direction = -1,limits = lims)
v2 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod1b',]),aes(fill=ci.width)) +
  ggtitle('Simple') +
  scale_fill_scico(name = 'Width of 90% credible intervals',palette = 'lajolla',direction = -1,limits = lims)
v3 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod2b',]),aes(fill=ci.width)) +
  ggtitle('SASW') +
  scale_fill_scico(name = 'Width of 90% credible intervals',palette = 'lajolla',direction = -1,limits = lims)

row2 <- ggarrange(plotlist = list(v1,v2,v3),common.legend = T,nrow=1)

lims <- range(mods_res$lower.1)

d1 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='std',]),aes(fill=lower.1)) +
  ggtitle('Standard') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims) +
  theme(legend.position = 'none')

d2 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod1b',]),aes(fill=lower.1)) +
  ggtitle('Simple') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims)+
  theme(legend.position = 'none')

d3 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod2b',]),aes(fill=lower.1)) +
  ggtitle('SASW') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims)+
  theme(legend.position = 'none')

row3 <- ggarrange(plotlist = list(d1,d2,d3),nrow=1)

lims <- range(mods_res$lower.25)

q1 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='std',]),aes(fill=lower.25)) +
  ggtitle('Standard') + scale_fill_scico(name = 'Probability ',palette = 'devon',direction = -1,limits = lims)

q2 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod1b',]),aes(fill=lower.25)) +
  ggtitle('Simple') + scale_fill_scico(name = 'Probability ',palette = 'devon',direction = -1,limits = lims)

q3 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod2b',]),aes(fill=lower.25)) +
  ggtitle('SASW') + scale_fill_scico(name = 'Probability ',palette = 'devon',direction = -1,limits = lims)

row4 <- ggarrange(plotlist = list(q1,q2,q3),common.legend = T,nrow=1,legend = 'bottom')

# 750 x 700
ggarrange(row1,row2,ncol = 1)

row3 <- annotate_figure(row3,
                        top = text_grob("Posterior probability of membership in lowest decile", 
                                        size = 16, hjust = 0, x=0))
row4 <- annotate_figure(row4,
                        top = text_grob("Posterior probability of membership in lowest quartile", 
                                        size = 16, hjust = 0, x = 0))
ggarrange(row3,row4,ncol = 1,common.legend = T,heights = c(1,1.15))

pdf('test.pdf')
## STUNTING (near 0)
lims <- range(mods_res$stunt)

s1 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='std',]),aes(fill=stunt)) +
  ggtitle('Standard') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims) +
  theme(legend.position = 'none')

s2 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod1b',]),aes(fill=stunt)) +
  ggtitle('Simple') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims)+
  theme(legend.position = 'none')

s3 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod2b',]),aes(fill=stunt)) +
  ggtitle('SASW') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims)+
  theme(legend.position = 'none')

print(ggarrange(plotlist = list(s1,s2,s3),common.legend = T,nrow=1,legend = 'bottom'))

## HAZ < 0 (near 1)
lims <- range(mods_res$less0)

s1 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='std',]),aes(fill=less0)) +
  ggtitle('Standard') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims) +
  theme(legend.position = 'none')

s2 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod1b',]),aes(fill=less0)) +
  ggtitle('Simple') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims)+
  theme(legend.position = 'none')

s3 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod2b',]),aes(fill=less0)) +
  ggtitle('SASW') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims)+
  theme(legend.position = 'none')

print(ggarrange(plotlist = list(s1,s2,s3),common.legend = T,nrow=1,legend = 'bottom'))

## HAZ < -1 
lims <- range(mods_res$lessneg1)

s1 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='std',]),aes(fill=lessneg1)) +
  ggtitle('Standard') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims) +
  theme(legend.position = 'none')

s2 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod1b',]),aes(fill=lessneg1)) +
  ggtitle('Simple') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims)+
  theme(legend.position = 'none')

s3 <- clean_map_theme + geom_sf(data = merge(poly.adm2, mods_res[mods_res$model=='Mod2b',]),aes(fill=lessneg1)) +
  ggtitle('SASW') + scale_fill_scico(name = '',palette = 'devon',direction = -1,limits = lims)+
  theme(legend.position = 'none')


print(ggarrange(plotlist = list(s1,s2,s3),common.legend = T,nrow=1,legend = 'bottom'))

dev.off()


### model results scatter plots -------

dir.est.tmp <- admin2.dir.stable %>% 
  rename(direct.est = mean) %>% mutate(direct.sd = sqrt(variance)) %>%
  dplyr::select(admin2,direct.est,direct.sd)

lims <- range(mods_res$mean,dir.est.tmp$direct.est,na.rm=T)

mean_scatter <- merge(mods_res,dir.est.tmp,by='admin2') %>% ggplot() + geom_point(aes(direct.est,mean),size=0.5) + 
  xlim(lims) + ylim(lims) +ggtitle('A. Mean') +
  xlab('Direct') + ylab('Model-based') +
  facet_wrap(~model,labeller=as_labeller(model_labs)) + geom_abline(intercept = 0,slope=1)

lims <- range(mods_res$sd,dir.est.tmp$direct.sd,na.rm=T)

var_scatter <- merge(mods_res,dir.est.tmp,by='admin2') %>% ggplot() + geom_point(aes(direct.sd,sd),size=0.5) + 
  xlim(lims) + ylim(lims) + ggtitle('B. Standard deviation') +
  xlab('Direct') + ylab('Model-based') +
  facet_wrap(~model,labeller=as_labeller(model_labs)) + geom_abline(intercept = 0,slope=1)

# 760 x 555
ggarrange(plotlist = list(mean_scatter,var_scatter),nrow=2)

### line plots ------

tmp <- mods_res %>% filter(model=='Mod1b')
tmp <- tmp[order(tmp$mean),]
tmp$area_order <- 1:n_admin2

order_key <- tmp[,c('admin2','area_order')]

mods_res <- merge(mods_res,order_key,by='admin2')




mods_res <- mods_res %>% 
  mutate(special_areas = #ifelse(admin2 %in% c(90,200,217,218),1,
                        ifelse(admin2 %in% mods_res[mods_res$model=='std' & mods_res$ci.width<0.15,]$admin2,1,
                                ifelse(admin2 %in% mods_res[mods_res$model=='std' & mods_res$ci.width>1,]$admin2,2,0)))
# 700 x 740
mods_res %>% ggplot() + 
  geom_segment(aes(x = area_order, xend = area_order, y=lower, yend = upper,colour = factor(special_areas)),lwd=0.5) + theme_bw() +
  geom_point(aes(area_order,mean),color='grey40',size=0.25)+ 
  xlab('Second adminstrative area (ordered by mean estimate under Simple model)') + ylab('Height-for-age z-score') +
  scale_color_manual(values = c('0'='grey50','1'='purple3','2'='orange3')) +
  theme(axis.text.x = element_blank(),legend.position = 'none') +
  facet_wrap(~model,labeller=as_labeller(model_labs),nrow=3)

# CALIBRATE parameters for simulations ----------

hyperpc.bym2 = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                    phi = list(prior = "pc", param = c(0.5, 0.5)))

fh2.bym2 <- inla(mean ~ f(admin2,model='bym2',hyper=hyperpc.bym2,graph = admin2.mat,scale.model = T)+ urb_frac,
                 family = "gaussian",
                 quantiles = c(0.05,0.95),
                 data = admin2.dir.stable,
                 scale = 1/admin2.dir.stable$variance,
                 control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
                 control.compute = list(config = TRUE))

std.fit <- fh2.bym2$summary.fitted.values
std.fit$admin2 <- 1:n_admin2

# urban effect on mean
dir.dat$urban <- ifelse(dir.dat$v025==2,0,dir.dat$v025)
tmp <- merge(dir.dat[,c('admin1','admin1.char','admin2','admin2.char','cluster','urban','v005','value')],std.fit[,c('mean','admin2')])

tmp %>% group_by(urban) %>% summarise(obs_mean = mean(value))
tmp %>% group_by(urban,admin2) %>% summarise(obs_mean = mean(value)) %>% 
  pivot_wider(id_cols=admin2,values_from = obs_mean,names_from=urban) %>% mutate(diff=`1`-`0`) %>% 
  summarise(mean(diff,na.rm=T),var(diff,na.rm=T))

# overall variance of residuals
tmp %>% summarise(log_sig2 = log(var(mean-value)))
tmp %>% group_by(urban) %>% summarise(log_sig2 = log(var(mean-value)))

# urban effect on variance
tmp %>% group_by(urban,admin2) %>% summarise(log_sig2 = log(var(mean-value))) %>% 
  pivot_wider(id_cols=admin2,values_from = log_sig2,names_from=urban) %>% mutate(diff=`1`-`0`)%>% 
  summarise(mean(diff,na.rm=T),var(diff,na.rm=T))

unit.model <- inla(value ~ f(admin2,model='bym2',hyper=hyperpc.bym2,graph = admin2.mat,scale.model = T)+
                     f(cluster,model='iid'),
                   data = dir.dat)
summary(unit.model)

#estimate cluster correlation:
(1/16.21)/(1/0.67 + 1/16.21)

