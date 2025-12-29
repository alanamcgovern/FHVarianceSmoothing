suppressMessages({
  library(dplyr)
  library(cmdstanr)
  library(Matrix)
  library(readr)
})

args <- commandArgs(trailingOnly=TRUE)

setting <- as.numeric(args[1])
k <- as.numeric(args[2])

mod1 <- cmdstan_model(stan_file = "Strat_noKappa.stan")
#mod1 <- cmdstan_model(stan_file = "Stan/Strat_noKappa.stan")

setwd(paste0('Sim',setting))

load(file='objects.rda')
load(file='cmat_admin2.rda')
load(file='sampled_clusters.rda')
load(file='direct.rda')

folders <- c('Mod2','Mod2a')
for(f in folders){
  if(!dir.exists(paste0(f,'/Hyper'))){
    dir.create(paste0(f,'/Hyper'), recursive = TRUE)
  }
}

admin.key <- objects$admin.key
n_admin2 <- nrow(admin.key)

mean_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac','X1','X2','X3')])
var_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac','X4','X5','X6')])

dir.dat <- direct[[k]]
dir.dat <- merge(dir.dat,cmat_admin2)

change_id <- which(dir.dat$variance < 1e-10)
if(length(change_id)>0){
  dir.dat[change_id,]$variance <- NA
  dir.dat[change_id,]$mean <- NA
}

# get input for models --------
data_areas <- which(!is.na(dir.dat$mean))
Cons_exact <- v_hat_scaled_exact <- q_start <- rep(NA,n_admin2)
q <- q_id <- nu <- NULL
tol <- 1e-10
for(area in data_areas){
  tmp <- sampled_clusters[[k]] %>% filter(admin1 == admin.key[admin.key$admin2 == area,]$admin1) %>% arrange(urban)# %>% mutate(wt=wt/1e3)
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
  
  v_hat_scaled_exact[area] <- dir.dat$variance[area]*as.numeric(t(omega)%*%S%*%omega)
  
  Cons_exact[area] <- as.numeric(t(omega)%*%S%*%omega)/sum(omega)^2

}

data_list_exact = list(m=n_admin2,
                       m_data=length(data_areas),
                       data_areas = data_areas,
                       y=dir.dat$mean[data_areas],
                       v_hat_scaled = v_hat_scaled_exact[data_areas],
                       lenq = length(q),
                       Cons = Cons_exact[data_areas],
                       q_id = q_id,
                       q_start = q_start[data_areas],
                       q_per_area = as.vector(table(q_id)),
                       q = q, nu = nu,
                       p_mean = ncol(mean_model_matrix_admin2),
                       X = mean_model_matrix_admin2,
                       bym2_mean = 1,
                       # all bym2 stuff
                       N_edges = nrow(objects$nodes2),
                       node_1 = objects$nodes2$node1,
                       node_2 = objects$nodes2$node2,
                       car_scale = objects$car_scale2) 

# fit with BYM2 and covariates -----------

data_list_exact$p_var = ncol(var_model_matrix_admin2)
data_list_exact$Z = var_model_matrix_admin2
data_list_exact$bym2_var = 1

fit4 <- mod1$sample(
  data = data_list_exact,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  show_messages = F,
  show_exceptions = F,
  refresh=00)

theta_fit4 <- as.matrix(unclass(fit4$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod2',
                     mean = apply(theta_fit4,2,mean),
                     upper = apply(theta_fit4,2,quantile,probs=0.95),
                     lower = apply(theta_fit4,2,quantile,probs=0.05)),file=paste0('Mod2/Result',k,'.csv'))

v_fit4 <- (as.matrix(unclass(fit4$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod2',
                     mean = apply(v_fit4,2,mean)),file=paste0('Mod2/Hyper/Var_Result',k,'.csv'))

hyper_fit4 <- as.matrix(unclass(fit4$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2','df','Q'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod2',
                     param = colnames(hyper_fit4),
                     mean = colMeans(hyper_fit4)),
          file=paste0('Mod2/Hyper/Hyper_Result',k,'.csv'))

rm(fit4,theta_fit4,v_fit4,hyper_fit4)



# fit with IID and no covariates (except urban prop) -----------

data_list_exact$p_var = 2
data_list_exact$Z = cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac')])
data_list_exact$bym2_var = 0

fit4 <- mod1$sample(
  data = data_list_exact,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  show_messages = F,
  show_exceptions = F,
  refresh=00)

theta_fit4 <- as.matrix(unclass(fit4$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod2a',
                     mean = apply(theta_fit4,2,mean),
                     upper = apply(theta_fit4,2,quantile,probs=0.95),
                     lower = apply(theta_fit4,2,quantile,probs=0.05)),file=paste0('Mod2a/Result',k,'.csv'))

v_fit4 <- (as.matrix(unclass(fit4$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod2a',
                     mean = apply(v_fit4,2,mean)),file=paste0('Mod2a/Hyper/Var_Result',k,'.csv'))

hyper_fit4 <- as.matrix(unclass(fit4$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2','df','Q'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod2a',
                     param = colnames(hyper_fit4),
                     mean = colMeans(hyper_fit4)),
          file=paste0('Mod2a/Hyper/Hyper_Result',k,'.csv'))

rm(fit4,theta_fit4,v_fit4,hyper_fit4)


