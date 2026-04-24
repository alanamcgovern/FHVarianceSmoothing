suppressMessages({
  library(dplyr)
  library(cmdstanr)
  library(Matrix)
  library(readr)
  library(data.table)
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
#load(file='cluster_frame.rda')
load(file='direct.rda')

#domain_probs <- cluster_frame[, .N, by = .(strata, admin2, urban)]
#domain_probs[, prop := N / sum(N), by = strata]

folders <- c('Mod4','Mod4a')
for(f in folders){
  if(!dir.exists(paste0(f,'/Hyper'))){
    dir.create(paste0(f,'/Hyper'), recursive = TRUE)
  }
  if(!dir.exists(paste0(f,'/Diagnostic'))){
    dir.create(paste0(f,'/Diagnostic'), recursive = TRUE)
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
cons_exact <- v_hat_scaled_exact <- q_start <- rep(NA,n_admin2)
q <- q_id <- nu <- NULL
tol <- 1e-10
for(area in data_areas){
  tmp_D <- sampled_clusters[[k]] %>% filter(admin2 == area) %>% group_by(urban) %>% summarise(N_D = n())
  #  if(sum(tmp_D$N_D>0)==2){
  
  tmp <- sampled_clusters[[k]] %>% filter(admin1 == admin.key[admin.key$admin2 == area,]$admin1) %>% arrange(urban)# %>% mutate(wt=wt/1e3)
  tmp$wt <- pmin(tmp$wt,quantile(tmp$wt,0.95))
  
  # omega_full <- tmp$n*tmp$wt
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
  
  v_hat_scaled_exact[area] <- dir.dat$variance[area]
  
  #  p_h <- c(rep(domain_probs[admin2==area & urban== 0,]$prop,N[1]),rep(domain_probs[admin2==area & urban== 1,]$prop,N[2]))
  
  # cons_exact[area] <- (t(omega_full)%*%diag(p_h)%*%S%*%omega_full)/sum(p_h*omega_full)^2
  cons_exact[area] <- (t(omega)%*%S%*%omega)/sum(omega)^2
  # }
  
}

data_areas <- which(!is.na(cons_exact))

data_list = list(m=n_admin2,
                       m_data=length(data_areas),
                       data_areas = data_areas,
                       y=dir.dat$mean[data_areas],
                       v_hat_scaled = v_hat_scaled_exact[data_areas],
                       Cons = cons_exact[data_areas],
                       lenq = length(q),
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
                       car_scale = objects$car_scale2,
                       bias_adj = 1) 

# fit with BYM2 and covariates -----------

data_list$p_var = ncol(var_model_matrix_admin2)
data_list$Z = var_model_matrix_admin2
data_list$bym2_var = 1

fit1 <- mod1$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
  show_messages = F,
  show_exceptions = F,
  refresh=00)

theta_fit1 <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod4',
                     mean = apply(theta_fit1,2,mean),
                     upper = apply(theta_fit1,2,quantile,probs=0.95),
                     lower = apply(theta_fit1,2,quantile,probs=0.05)),file=paste0('Mod4/Result',k,'.csv'))

v_fit1 <- (as.matrix(unclass(fit1$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod4',
                     mean = apply(v_fit1,2,mean)),file=paste0('Mod4/Hyper/Var_Result',k,'.csv'))

hyper_fit1 <- as.matrix(unclass(fit1$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2','Q'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod4',
                     param = colnames(hyper_fit1),
                     mean = colMeans(hyper_fit1)),
          file=paste0('Mod4/Hyper/Hyper_Result',k,'.csv'))

df_fit1 <- as.matrix(unclass(fit1$draws(variables = c('df'), format = "matrix")))
v_raw <- as.matrix(unclass(fit1$draws(variables = c('v_raw'), format = "matrix")))

scale_fit1 <- sapply(1:ncol(v_raw),function(j){dir.dat$variance[data_areas[j]]/v_raw[,j]})

write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod4',
                     df = colMeans(df_fit1),
                     scale = colMeans(scale_fit1)),
          file=paste0('Mod4/Hyper/Chisq_params',k,'.csv'))

# diagnostics
summary_df <- fit1$summary()
diag_df    <- fit1$diagnostic_summary()

write.csv(data.frame(sim=k,
                     model = 'Mod4',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail),
          file=paste0('Mod4/Diagnostic/Diag_Summary_',k,'.csv'))

prob_params <- summary_df %>% filter((!is.na(rhat) & rhat > 1.01) | ess_bulk < 400 | ess_tail < 100)
if(nrow(prob_params)>0){
  write.csv(prob_params,file=paste0('Mod4/Diagnostic/Diag_detect_',k,'.csv'))
}




# fit with IID and no covariates (except urban prop) -----------

data_list$p_var = 2
data_list$Z = cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac')])
data_list$bym2_var = 0

fit2 <- mod1$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
  show_messages = F,
  show_exceptions = F,
  refresh=00)

theta_fit2 <- as.matrix(unclass(fit2$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod4a',
                     mean = apply(theta_fit2,2,mean),
                     upper = apply(theta_fit2,2,quantile,probs=0.95),
                     lower = apply(theta_fit2,2,quantile,probs=0.05)),file=paste0('Mod4a/Result',k,'.csv'))

v_fit2 <- (as.matrix(unclass(fit2$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod4a',
                     mean = apply(v_fit2,2,mean)),file=paste0('Mod4a/Hyper/Var_Result',k,'.csv'))

hyper_fit2 <- as.matrix(unclass(fit2$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2','Q'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod4a',
                     param = colnames(hyper_fit2),
                     mean = colMeans(hyper_fit2)),
          file=paste0('Mod4a/Hyper/Hyper_Result',k,'.csv'))

df_fit2 <- as.matrix(unclass(fit2$draws(variables = c('df'), format = "matrix")))
v_raw <- as.matrix(unclass(fit2$draws(variables = c('v_raw'), format = "matrix")))

scale_fit2 <- sapply(1:ncol(v_raw),function(j){dir.dat$variance[data_areas[j]]/v_raw[,j]})

write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod4a',
                     df = colMeans(df_fit2),
                     scale = colMeans(scale_fit2)),
          file=paste0('Mod4a/Hyper/Chisq_params',k,'.csv'))

# diagnostics
summary_df <- fit2$summary()
diag_df    <- fit2$diagnostic_summary()

write.csv(data.frame(sim=k,
                     model = 'Mod4a',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail),
          file=paste0('Mod4a/Diagnostic/Diag_Summary_',k,'.csv'))

prob_params <- summary_df %>% filter((!is.na(rhat) & rhat > 1.01) | ess_bulk < 400 | ess_tail < 100)
if(nrow(prob_params)>0){
  write.csv(prob_params,file=paste0('Mod4a/Diagnostic/Diag_detect_',k,'.csv'))
}

