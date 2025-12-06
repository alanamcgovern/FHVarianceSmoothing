suppressMessages({
  library(dplyr)
  library(cmdstanr)
  library(Matrix)
  library(readr)
})

args <- commandArgs(trailingOnly=TRUE)

setting <- as.numeric(args[1])
k <- as.numeric(args[2])

mod1 <- cmdstan_model(stan_file = "Strat_oneKappa.stan")

setwd(paste0('Sim',setting))

load(file='objects.rda')
load(file='cmat_admin2.rda')
load(file='sampled_clusters.rda')
load(file='direct.rda')

folders <- c('Mod3_1','Mod3_2')
for(f in folders){
  if(!dir.exists(paste0(f,'/Hyper'))){
    dir.create(paste0(f,'/Hyper'), recursive = TRUE)
  }
}

admin.key <- objects$admin.key
n_admin2 <- nrow(admin.key)
n_admin1 <- length(unique(admin.key$admin1))

kappa_seq <- seq(-0.9,4,0.1)

mean_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac','X1','X2','X3')])
mean_model_matrix_admin1 <- cbind(rep(1,n_admin1),objects$A_2to1 %*% as.matrix(cmat_admin2[,c('urb_frac','X1','X2','X3')]))

var_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('X4','X5','X6')])
var_model_matrix_admin1 <- cbind(rep(1,n_admin1),objects$A_2to1 %*% as.matrix(cmat_admin2[,c('X4','X5','X6')]))

dir.dat <- direct[[k]]
dir.dat <- merge(dir.dat,cmat_admin2)

change_id <- which(dir.dat$variance < 1e-10)
if(length(change_id)>0){
  dir.dat[change_id,]$variance <- NA
  dir.dat[change_id,]$mean <- NA
}

# get input for models --------
data_areas <- which(!is.na(dir.dat$mean))
v_hat_scaled_exact <- q_start <- rep(NA,n_admin2)
Cons_exact <- matrix(NA,n_admin2,2)
q <- q_id <- nu_main <- nu_urban <- NULL
for(area in data_areas){
  tmp <- sampled_clusters[[k]] %>% filter(admin1 == admin.key[admin.key$admin2 == area,]$admin1) %>% arrange(urban) #%>% mutate(wt=wt/1e3)
 # tmp$wt <- pmin(tmp$wt,quantile(tmp$wt,0.95))
  tmp[tmp$admin2 != area,]$wt <- 0
  N_D <- c(sum((1-tmp$urban)*tmp$wt>0),sum(tmp$urban*tmp$wt>0))
  N <- c(sum(1-tmp$urban),sum(tmp$urban))
  
  omega <- tmp$n*tmp$wt
  Cons_exact[area,] <- c(sum(tmp$wt^2*tmp$n), sum(tmp$urban*tmp$wt^2*tmp$n))
  v_hat_scaled_exact[area] <- dir.dat$variance[area]
  
  B_comps <- lapply(1:2,function(h){
    if(N[h]>0){
      return(N[h]/(N[h]-1)*(diag(1,N[h]) - 1/N[h]*matrix(1,N[h],N[h])))
    }
  })
  
  which.not.null <- c(1:2)[!sapply(B_comps, is.null)]
  B <- bdiag(B_comps[which.not.null])
  W <- diag(omega)*(diag(1,sum(N)) - (rep(1,sum(N))%*%t(omega))/sum(omega))
  C <- t(W)%*%B%*%W
  S <- diag(1/tmp$n)
  
  if(is.null(q)){
    q_start[area] <- 1
  }else{
    q_start[area] <- nrow(q) + 1
  }
  r <- rankMatrix(sqrt(S)%*%C%*%sqrt(S))
  q_id <- c(q_id,rep(area,r))
  q_tmp <- nu_main_tmp <- nu_urban_tmp <- matrix(NA,r,length(kappa_seq))
  for(kappa in kappa_seq){
    kappa_id <- which(kappa==kappa_seq)
    perm <- diag(1,sum(N)) + kappa*diag(tmp$urban)
    out <- eigen(sqrt(perm%*%S)%*%C%*%sqrt(perm%*%S))
   
    q_tmp[,kappa_id] <- out$values[1:r]
    U <- out$vectors[,1:r]
    
    nu_main_tmp[,kappa_id] <- t(U)%*%sqrt(perm%*%S)%*%rep(1,sum(N))
    nu_urban_tmp[,kappa_id] <- t(U)%*%sqrt(perm%*%S)%*%tmp$urban
  }
  
  q <- rbind(q,q_tmp)
  nu_main <- rbind(nu_main,nu_main_tmp)
  nu_urban <- rbind(nu_urban,nu_urban_tmp)
  
}

# "Exact" smoothing ----------

data_list_exact = list(m=n_admin2,
                       m_data=length(data_areas),
                       data_areas = data_areas,
                       y=dir.dat$mean[data_areas],
                       v_hat_scaled = v_hat_scaled_exact[data_areas],
                       Cons=Cons_exact[data_areas,],
                       lenq = nrow(q),
                       q_id = q_id,
                       q_start = q_start[data_areas],
                       q_per_area = as.vector(table(q_id)),
                       len_kappa_seq = length(kappa_seq),
                       kappa_seq = kappa_seq,
                       q_mat = q, nu_main = nu_main, nu_urban = nu_urban,
                       p_mean = ncol(mean_model_matrix_admin2),
                       X = mean_model_matrix_admin2,
                       bym2_mean = 1,
                       # all bym2 stuff
                       N_edges1 = nrow(objects$nodes1),
                       node1_1 = objects$nodes1$node1,
                       node1_2 = objects$nodes1$node2,
                       car_scale1 = objects$car_scale1,
                       N_edges2 = nrow(objects$nodes2),
                       node2_1 = objects$nodes2$node1,
                       node2_2 = objects$nodes2$node2,
                       car_scale2 = objects$car_scale2) 
### Type 3_1 ----------

data_list_exact$sig2_level = 1
data_list_exact$len_sig2 = n_admin1
data_list_exact$sig2_id = admin.key[order(admin.key$admin2),]$admin1[data_areas]
data_list_exact$p_var = ncol(var_model_matrix_admin1)
data_list_exact$Z = var_model_matrix_admin1
data_list_exact$bym2_var = 0

fit1 <- mod1$sample(
  data = data_list_exact,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  show_messages = F,
  show_exceptions = F,
  refresh=00)

theta_fit1 <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod3_1',
                     mean = apply(theta_fit1,2,mean),
                     upper = apply(theta_fit1,2,quantile,probs=0.95),
                     lower = apply(theta_fit1,2,quantile,probs=0.05)),file=paste0('Mod3_1/Result',k,'.csv'))

sig2_fit1 <- exp(as.matrix(unclass(fit1$draws(variables = c("log_sig2"), format = "matrix"))))
kappa_fit1 <- as.matrix(unclass(fit1$draws(variables = c("kappa"), format = "matrix")))

write.csv(data.frame(sim=k,
                     area = 1:n_admin1,
                     model = 'Mod3_1',
                     mean = apply(sig2_fit1,2,mean)),
          file=paste0('Mod3_1/Hyper/Sig2_Result',k,'.csv'))

hyper_fit1 <- as.matrix(unclass(fit1$draws(variables = c("beta",'gamma','sig_u','phi1','kappa','df'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod3_1',
                     param = colnames(hyper_fit1),
                     mean = colMeans(hyper_fit1)),
          file=paste0('Mod3_1/Hyper/Hyper_Result',k,'.csv'))

rm(fit1,theta_fit1,sig2_fit1,hyper_fit1)

### Type 3_2 ----------

data_list_exact$sig2_level = 2
data_list_exact$len_sig2 = n_admin2
data_list_exact$sig2_id = data_areas
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
                     model = 'Mod3_2',
                     mean = apply(theta_fit4,2,mean),
                     upper = apply(theta_fit4,2,quantile,probs=0.95),
                     lower = apply(theta_fit4,2,quantile,probs=0.05)),file=paste0('Mod3_2/Result',k,'.csv'))

sig2_fit4 <- exp(as.matrix(unclass(fit4$draws(variables = c("log_sig2"), format = "matrix"))))
kappa_fit4 <- as.matrix(unclass(fit4$draws(variables = c("kappa"), format = "matrix")))

write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod3_2',
                     mean = apply(sig2_fit4,2,mean)),
          file=paste0('Mod3_2/Hyper/Sig2_Result',k,'.csv'))

hyper_fit4 <- as.matrix(unclass(fit4$draws(variables = c("beta",'gamma','sig_u','phi1','phi2','kappa','df'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod3_2',
                     param = colnames(hyper_fit4),
                     mean = colMeans(hyper_fit4)),
          file=paste0('Mod3_2/Hyper/Hyper_Result',k,'.csv'))

rm(fit4,theta_fit4,sig2_fit4,hyper_fit4)




