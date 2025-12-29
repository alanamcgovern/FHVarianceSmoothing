
suppressMessages({
  library(dplyr)
  library(cmdstanr)
  library(Matrix)
  library(readr)
})

args <- commandArgs(trailingOnly=TRUE)

setting <- as.numeric(args[1])
k <- as.numeric(args[2])

mod <- cmdstan_model(stan_file = "NonStrat.stan")
##mod <- cmdstan_model(stan_file = "Stan/NonStrat.stan")

setwd(paste0('Sim',setting))

load(file='objects.rda')
load(file='cmat_admin2.rda')
load(file='sampled_clusters.rda')
load(file='direct.rda')

folders <- c('Mod1.5','Mod1.5a')
for(f in folders){
  if(!dir.exists(paste0(f,'/Hyper'))){
    dir.create(paste0(f,'/Hyper'), recursive = TRUE)
  }
}

admin.key <- objects$admin.key
n_admin2 <- nrow(admin.key)
n_admin1 <- length(unique(admin.key$admin1))

mean_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac','X1','X2','X3')])
var_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac','X4','X5','X6')])

dir.dat <- direct[[k]]
dir.dat <- merge(dir.dat,cmat_admin2)

change_id <- which(dir.dat$variance < 1e-10)
if(length(change_id)>0){
  dir.dat[change_id,]$variance <- NA
  dir.dat[change_id,]$mean <- NA
}

# get input 
data_areas <- which(!is.na(dir.dat$mean))
Cons <- df <- v_hat_scaled <- rep(NA,n_admin2)
tol <- 1e-10

for(area in data_areas){
  tmp <- sampled_clusters[[k]] %>% filter(admin2 == area) %>% arrange(urban)# %>% mutate(wt=wt/1e3)
  N <- c(sum(1-tmp$urban),sum(tmp$urban))
  
  tmp$omega <- tmp$n*tmp$wt
  if(N[1]>0)
    tmp[tmp$urban==0,]$omega <- median(tmp[tmp$urban==0,]$omega)
  if(N[2]>0)
    tmp[tmp$urban==1,]$omega <- median(tmp[tmp$urban==1,]$omega)
 
  Q_list <- lapply(1:2,function(h){
    if(N[h]>1){
      M <- diag(1,N[h]) - 1/N[h]*matrix(1,N[h],N[h])
      S <- diag(1/tmp[tmp$urban==h-1,]$n)
      out <- eigen(sqrt(S)%*%M%*%sqrt(S), symmetric = TRUE)
      q <- out$values[out$values>tol]
      Q1 <- N[h]/(N[h]-1)*unique(tmp[tmp$urban==h-1,]$omega)^2*sum(q)
      Q2 <- N[h]^2/(N[h]-1)^2*unique(tmp[tmp$urban==h-1,]$omega)^4*sum(q^2)
      return(c(Q1,Q2))
    }else{
      return(c(0,0))
    }
  })
  Q_mat <- do.call(rbind,Q_list)

  Cons_comp <- sapply(1:2,function(h){
    if(N[h]>0){
      return(unique(tmp[tmp$urban==h-1,]$omega)^2*sum(1/tmp[tmp$urban==h-1,]$n))
    }else{
      return(0)
    }
  })
  
  Cons[area] <- sum(Cons_comp)/sum(tmp$omega)^2
  if(sum(Q_mat > 0)==0){
    df[area] <- 1
    v_hat_scaled[area] <- dir.dat$variance[area]*sum(Cons_comp)
  }else{
    df[area] <- sum(Q_mat[,1])^2/sum(Q_mat[,2])
    v_hat_scaled[area] <- dir.dat$variance[area]*sum(Cons_comp)*sum(Q_mat[,1])/sum(Q_mat[,2])
  }

}

data_list = list(m=n_admin2,
                 m_data=length(data_areas),
                 data_areas = data_areas,
                 y=dir.dat$mean[data_areas],
                 v_hat_scaled = v_hat_scaled[data_areas],
                 df=df[data_areas],
                 Cons = Cons[data_areas],
                 p_mean = ncol(mean_model_matrix_admin2),
                 X = mean_model_matrix_admin2,
                 bym2_mean = 1,
                 # all bym2 stuff
                 N_edges = nrow(objects$nodes2),
                 node_1 = objects$nodes2$node1,
                 node_2 = objects$nodes2$node2,
                 car_scale = objects$car_scale2) 

# use BYM2 and covariates -------

data_list$bym2_var <- 1
data_list$p_var = ncol(var_model_matrix_admin2)
data_list$Z = var_model_matrix_admin2

fit1 <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  show_messages = F,
  show_exceptions = F,
  refresh=0)

theta_fit1 <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod1.5',
                     mean = apply(theta_fit1,2,mean),
                     upper = apply(theta_fit1,2,quantile,probs=0.95),
                     lower = apply(theta_fit1,2,quantile,probs=0.05)),file=paste0('Mod1.5/Result',k,'.csv'))

v_fit1 <- (as.matrix(unclass(fit1$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod1.5',
                     mean = apply(v_fit1,2,mean)),file=paste0('Mod1.5/Hyper/Var_Result',k,'.csv'))

hyper_fit1 <- as.matrix(unclass(fit1$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod1.5',
                     param = colnames(hyper_fit1),
                     mean = colMeans(hyper_fit1)),
          file=paste0('Mod1.5/Hyper/Hyper_Result',k,'.csv'))

rm(fit1,theta_fit1,v_fit1,hyper_fit1)


# use IID and no covariates (except urban/rural) -------

data_list$bym2_var <- 0
data_list$p_var = 2
data_list$Z = cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac')])

fit1 <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  show_messages = F,
  show_exceptions = F,
  refresh=0)

theta_fit1 <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod1.5a',
                     mean = apply(theta_fit1,2,mean),
                     upper = apply(theta_fit1,2,quantile,probs=0.95),
                     lower = apply(theta_fit1,2,quantile,probs=0.05)),file=paste0('Mod1.5a/Result',k,'.csv'))

v_fit1 <- (as.matrix(unclass(fit1$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod1.5a',
                     mean = apply(v_fit1,2,mean)),file=paste0('Mod1.5a/Hyper/Var_Result',k,'.csv'))

hyper_fit1 <- as.matrix(unclass(fit1$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod1.5a',
                     param = colnames(hyper_fit1),
                     mean = colMeans(hyper_fit1)),
          file=paste0('Mod1.5a/Hyper/Hyper_Result',k,'.csv'))

rm(fit1,theta_fit1,v_fit1,hyper_fit1)



