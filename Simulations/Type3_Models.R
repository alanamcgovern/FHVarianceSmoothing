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
#mod1 <- cmdstan_model(stan_file = "Stan/Strat_oneKappa.stan")

setwd(paste0('Sim',setting))

load(file='objects.rda')
load(file='cmat_admin2.rda')
load(file='sampled_clusters.rda')
load(file='direct.rda')

folders <- c('Mod3')
for(f in folders){
  if(!dir.exists(paste0(f,'/Hyper'))){
    dir.create(paste0(f,'/Hyper'), recursive = TRUE)
  }
}

admin.key <- objects$admin.key
n_admin2 <- nrow(admin.key)

kappa_seq <- seq(-0.9,4,0.1)

mean_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac','X1','X2','X3')])
var_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('X4','X5','X6')])

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
q <- q_id <- nu <- NULL
tol <- 1e-10
for(area in data_areas){
  tmp <- sampled_clusters[[k]] %>% filter(admin1 == admin.key[admin.key$admin2 == area,]$admin1) %>% arrange(urban) #%>% mutate(wt=wt/1e3)
  tmp$wt <- pmin(tmp$wt,quantile(tmp$wt,0.95))
  tmp[tmp$admin2 != area,]$wt <- 0
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
  W <- diag(omega)%*%(diag(1,sum(N)) - (rep(1,sum(N))%*%t(omega))/sum(omega))
  C <- t(W)%*%B%*%W
  S <- diag(1/tmp$n)
  
  if(is.null(q)){
    q_start[area] <- 1
  }else{
    q_start[area] <- nrow(q) + 1
  }
  
 t <-  lapply(kappa_seq,function(kappa){
    perm <- diag(1,sum(N)) + kappa*diag(tmp$urban)
    
    L <- chol(perm %*% S)              # S = L' L
    Li <- backsolve(L, diag(nrow(S)))  # L^{-1}
    
    # Form scale-free matrix A = S^{1/2} C S^{1/2}
    A <- t(L) %*% C %*% L
    
    out <- eigen(A, symmetric = TRUE)
    V <- Li %*% out$vectors
    keep_id <- which(out$values > tol)
   
    q_tmp <- out$values[keep_id]
    V <- V[,keep_id]
    
    nu_tmp <- sqrt(sum(tmp$wt^2*tmp$n))/sum(omega)*t(V)%*%tmp$urban
    
    return(list(q = out$values[keep_id],nu = nu_tmp))
  })
  
  q_id <- c(q_id,rep(area,length(t[[1]]$q)))
  
  q <- rbind(q,do.call(cbind,lapply(t,function(x){x$q})))
  nu <- rbind(nu,do.call(cbind,lapply(t,function(x){x$nu})))
  
}

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
                       q_mat = q, nu = nu,
                      #len_kappa_seq = 1,
                      #kappa_seq = 0,
                       #q_mat = q[,10], nu = nu[,10],
                       p_mean = ncol(mean_model_matrix_admin2),
                       X = mean_model_matrix_admin2,
                       bym2_mean = 1,
                       p_var = ncol(var_model_matrix_admin2),
                       Z = var_model_matrix_admin2,
                       bym2_var = 1,
                       # all bym2 stuff
                       N_edges = nrow(objects$nodes2),
                       node_1 = objects$nodes2$node1,
                       node_2 = objects$nodes2$node2,
                       car_scale = objects$car_scale2) 

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
                     model = 'Mod3',
                     mean = apply(theta_fit4,2,mean),
                     upper = apply(theta_fit4,2,quantile,probs=0.95),
                     lower = apply(theta_fit4,2,quantile,probs=0.05)),file=paste0('Mod3/Result',k,'.csv'))

v_fit4 <- exp(as.matrix(unclass(fit4$draws(variables = c("log_v"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod3',
                     mean = apply(v_fit4,2,mean)),
          file=paste0('Mod3/Hyper/Var_Result',k,'.csv'))

hyper_fit4 <- as.matrix(unclass(fit4$draws(variables = c("beta",'gamma','sig_u','phi1','phi2','kappa','df','delta'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod3',
                     param = colnames(hyper_fit4),
                     mean = colMeans(hyper_fit4)),
          file=paste0('Mod3/Hyper/Hyper_Result',k,'.csv'))

rm(fit4,theta_fit4,v_fit4,hyper_fit4)




