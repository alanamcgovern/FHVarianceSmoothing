setwd('/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing')

setting <- 3
mod <- cmdstan_model(stan_file = "Stan/NonStrat_fixedareas.stan")
mod2 <- cmdstan_model(stan_file = "Stan/Strat_fixedareas.stan")

setwd(paste0('Simulations/Sim',setting))

load(file='objects.rda')
load(file='cmat_admin2.rda')
load(file='sampled_clusters.rda')
load(file='direct.rda')
load(file='params.rda')
admin.key <- objects$admin.key

all_areas <- sort(sample.int(300,10))

true_v = sapply(all_areas,function(a){var(unlist(lapply(direct,function(x){x$mean[a]})))})

theta_est <- beta_est <- sig2_est <- v_est <- NULL
for(k in 1:50){
  cat(k,'\n')
  dir.dat <- direct[[k]]
  
  change_id <- which(dir.dat$variance < 1e-10)
  if(length(change_id)>0){
    dir.dat[change_id,]$mean <- NA
  }
  
  data_areas <- all_areas[which(!is.na(dir.dat[all_areas,]$mean))]
  
  mean_model_matrix_admin2 <- cbind(rep(1,length(data_areas)),cmat_admin2[data_areas,c('urb_frac','X1')])
  var_model_matrix_admin2 <- matrix(1,nrow=length(data_areas),ncol=1)
  
  # get input ----
  Cons_naive <- df_naive <- v_hat_scaled_naive <- cons_satt <- v_hat_scaled_satt <- q_start <-  rep(NA,length(all_areas))
  q <- q_id <- nu <- NULL
  tol <- 1e-10
  for(area in data_areas){
    tmp_D <- sampled_clusters[[k]] %>% filter(admin2 == area)
    N_D <- c(sum(1-tmp_D$urban),sum(tmp_D$urban))
    
    df_naive[area] <- pmax(1,sum(N_D) - sum(N_D>0))
    v_hat_scaled_naive[area] <- dir.dat$variance[area]*pmax(1,sum(N_D) - sum(N_D>0))
    Cons_naive[area] <- 1/(sum(N_D)*mean(tmp_D$n))
    
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
    q_id <- c(q_id,rep(which(area==data_areas),length(keep_id)))
    q <- c(q,q_tmp)
    
    V <- V[,keep_id]
    nu <- c(nu,t(V)%*%tmp$urban)
    
    v_hat_scaled_satt[area] <- dir.dat$variance[area]*as.numeric(t(omega)%*%S%*%omega)
    
    cons_satt[area] <- as.numeric(t(omega)%*%S%*%omega)/sum(omega)^2
  }
  
  # naive -----
  data_list = list(m = length(data_areas),
                   y=dir.dat$mean[data_areas],
                   v_hat_scaled = v_hat_scaled_naive[data_areas],
                   df=df_naive[data_areas],
                   Cons = Cons_naive[data_areas],
                   p_mean = ncol(mean_model_matrix_admin2), 
                   p_var = ncol(var_model_matrix_admin2),
                   X = mean_model_matrix_admin2,
                   Z = var_model_matrix_admin2)
  
  fit1 <- mod$sample(
    data = data_list,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    show_messages = F,
    show_exceptions = F,
    refresh=00)
  
  naive_theta_fit <- as.matrix(unclass(fit1$draws(variables = c('theta'),format = "matrix")))
  naive_beta_fit <- as.matrix(unclass(fit1$draws(variables = c('beta'),format = "matrix")))
  naive_sig2_fit <- exp(as.matrix(unclass(fit1$draws(variables = c('log_sig2'),format = "matrix"))))
  naive_v_fit <- as.matrix(unclass(fit1$draws(variables = c('v'),format = "matrix")))
  
  # satterwhaite ------
  data_list = list(m = length(data_areas),
                   y=dir.dat$mean[data_areas],
                   v_hat_scaled = v_hat_scaled_satt[data_areas],
                   lenq = length(q),
                   Cons =cons_satt[data_areas],
                   q_id = q_id,
                   q_start = q_start[data_areas],
                   q_per_area = as.vector(table(q_id)),
                   q = q, nu = nu,
                   p_mean = ncol(mean_model_matrix_admin2), 
                   p_var = ncol(var_model_matrix_admin2),
                   X = mean_model_matrix_admin2,
                   Z = var_model_matrix_admin2)
  
  fit1 <- mod2$sample(
    data = data_list,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    show_messages = F,
    show_exceptions = F,
    refresh=00)
  
  satt_theta_fit <- as.matrix(unclass(fit1$draws(variables = c('theta'),format = "matrix")))
  satt_beta_fit <- as.matrix(unclass(fit1$draws(variables = c('beta'),format = "matrix")))
  satt_sig2_fit <- exp(as.matrix(unclass(fit1$draws(variables = c('log_sig2'),format = "matrix"))))
  satt_v_fit <- as.matrix(unclass(fit1$draws(variables = c('v'),format = "matrix")))
  
  # record ----
  theta_est <- rbind(theta_est,
                     data.frame(sim = k, method = 'naive', area = data_areas,
                                mean = apply(naive_theta_fit,2,mean),
                                lower = apply(naive_theta_fit,2,quantile,0.05),
                                upper = apply(naive_theta_fit,2,quantile,0.95)),
                     data.frame(sim = k, method = 'satt', area = data_areas,
                                mean = apply(satt_theta_fit,2,mean),
                                lower = apply(satt_theta_fit,2,quantile,0.05),
                                upper = apply(satt_theta_fit,2,quantile,0.95)))
  
  beta_est <- rbind(beta_est,
                     data.frame(sim = k, method = 'naive', param = 1:3,
                                mean = apply(naive_beta_fit,2,mean),
                                lower = apply(naive_beta_fit,2,quantile,0.05),
                                upper = apply(naive_beta_fit,2,quantile,0.95)),
                     data.frame(sim = k, method = 'satt', param = 1:3,
                                mean = apply(satt_beta_fit,2,mean),
                                lower = apply(satt_beta_fit,2,quantile,0.05),
                                upper = apply(satt_beta_fit,2,quantile,0.95)))
  
  sig2_est <- rbind(sig2_est,
                 data.frame(sim = k, method = 'naive', area = data_areas,
                            mean = apply(naive_sig2_fit,2,mean),
                            lower = apply(naive_sig2_fit,2,quantile,0.05),
                            upper = apply(naive_sig2_fit,2,quantile,0.95)),
                 data.frame(sim = k, method = 'satt', area = data_areas,
                            mean = apply(satt_sig2_fit,2,mean),
                            lower = apply(satt_sig2_fit,2,quantile,0.05),
                            upper = apply(satt_sig2_fit,2,quantile,0.95)))
  
  v_est <- rbind(v_est,
                     data.frame(sim = k, method = 'naive', area = data_areas,
                                mean = apply(naive_v_fit,2,mean),
                                lower = apply(naive_v_fit,2,quantile,0.05),
                                upper = apply(naive_v_fit,2,quantile,0.95)),
                     data.frame(sim = k, method = 'satt', area = data_areas,
                                mean = apply(satt_v_fit,2,mean),
                                lower = apply(satt_v_fit,2,quantile,0.05),
                                upper = apply(satt_v_fit,2,quantile,0.95)))
  
}

theta_est$truth <- params$admin2_mean[all_areas]
beta_est$truth <- 0
sig2_est$truth <- params$sig2[all_areas]
#v_est$truth <- true_v[v_est$area]

theta_est %>% group_by(method,area,truth) %>% summarise(val = mean(mean)) %>% ggplot() + geom_point(aes(truth,val)) + 
  facet_grid(~method) +geom_abline(intercept = 0,slope = 1)
sig2_est %>% group_by(method,area,truth) %>% summarise(val = mean(mean)) %>% ggplot() + geom_point(aes(truth,val)) + 
  facet_grid(~method) +geom_abline(intercept = 0,slope = 1)
# v_est %>% group_by(method,area,truth) %>% summarise(val = mean(mean)) %>% ggplot() + geom_point(aes(truth,val)) + 
#   facet_grid(~method) +geom_abline(intercept = 0,slope = 1)

theta_est %>% mutate(cov= lower<truth & upper > truth,
                     width = upper - lower) %>% group_by(method) %>% summarise(mean(cov),mean(width))

beta_est %>% group_by(method,param) %>% summarise(mean(mean))
