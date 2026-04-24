suppressMessages({
  library(dplyr)
  library(cmdstanr)
  library(Matrix)
  library(readr)
  library(stringr)
})

setwd("/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing")

mod1 <- cmdstan_model(stan_file = "Stan/NonStrat.stan")
mod2 <- cmdstan_model(stan_file = "Stan/Strat_noKappa.stan")

time_log <- NULL
for(setting in 1:5){
  for(k in 1:5){
    cat(setting,'--',k,'\n')
    start.time0 <- Sys.time()
    
    setwd(paste0('/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/Sim',setting))
    
    load(file='objects.rda')
    load(file='cmat_admin2.rda')
    load(file='sampled_clusters.rda')
    load(file='direct.rda')
    
    
    admin.key <- objects$admin.key
    n_admin2 <- nrow(admin.key)
    n_admin1 <- length(unique(admin.key$admin1))
    
    mean_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac','X1','X2','X3')])
    var_model_matrix_admin2 <- cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac','X4','X5','X6')])
    
    dir.dat <- direct[[k]]
    dir.dat <- merge(dir.dat,cmat_admin2)
    
    change_id <- which(dir.dat$variance < 1e-5)
    if(length(change_id)>0){
      dir.dat[change_id,]$variance <- NA
      dir.dat[change_id,]$mean <- NA
    }
    
    ## NAIVE --------------------
    ## get input -------------------
    data_areas <- which(!is.na(dir.dat$mean))
    Cons_naive <- scale_naive <- df_naive <- v_hat_scaled_naive <- rep(NA,n_admin2)
    
    for(area in data_areas){
      
      tmp_D <- sampled_clusters[[k]] %>% filter(admin2 == area)
      N_D <- c(sum(1-tmp_D$urban),sum(tmp_D$urban))
      
      df_naive[area] <- pmax(1,sum(N_D) - sum(N_D>0))
      v_hat_scaled_naive[area] <- dir.dat$variance[area]*pmax(1,sum(N_D) - sum(N_D>0))
      Cons_naive[area] <- 1/(sum(N_D)*mean(tmp_D$n))
      scale_naive[area] <- 1/(sum(N_D)*mean(tmp_D$n)*pmax(1,sum(N_D) - sum(N_D>0)))
      
    }
    
    data_areas <- which(!is.na(df_naive))
    
    data_list = list(m=n_admin2,
                     m_data=length(data_areas),
                     data_areas = data_areas,
                     y=dir.dat$mean[data_areas],
                     v_hat_scaled = v_hat_scaled_naive[data_areas],
                     df=df_naive[data_areas],
                     Cons = Cons_naive[data_areas],
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
    
    start.time <- Sys.time()
    fit <- mod1$sample(
      data = data_list,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.99,
      show_messages = F,
      show_exceptions = F,
      refresh=0)
    
    time_log <- rbind(time_log,
                      data.frame(setting = setting, model = 'Mod1', time = Sys.time() - start.time))
    
    # use IID and no covariates (except urban/rural) -------
    
    data_list$bym2_var <- 0
    data_list$p_var = 2
    data_list$Z = cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac')])
    
    start.time <- Sys.time()
    fit <- mod1$sample(
      data = data_list,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.99,
      show_messages = F,
      show_exceptions = F,
      refresh=0)
    
    time_log <- rbind(time_log,
                      data.frame(setting = setting, model = 'Mod1a', time = Sys.time() - start.time))
    
    
    ## SASW -------------
    ## get input -----------
    data_areas <- which(!is.na(dir.dat$mean))
    cons_exact <- v_hat_scaled_exact <- q_start <-  omegaSomega <- rep(NA,n_admin2)
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
      
      omegaSomega[area] <- t(omega)%*%S%*%omega
      
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
      
      v_hat_scaled_exact[area] <- dir.dat$variance[area]*(t(omega)%*%S%*%omega)
      
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
                     bias_adj = 0) 
    
    # fit with BYM2 and covariates -----------
    
    data_list$p_var = ncol(var_model_matrix_admin2)
    data_list$Z = var_model_matrix_admin2
    data_list$bym2_var = 1
    
    start.time <- Sys.time()
    fit1 <- mod2$sample(
      data = data_list,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.99,
      show_messages = F,
      show_exceptions = F,
      refresh=0)
    Sys.time() - start.time
    time_log <- rbind(time_log,
                      data.frame(setting = setting, model = 'Mod2', time = Sys.time() - start.time))
    
    # fit with IID and no covariates (except urban prop) -----------
    
    data_list$p_var = 2
    data_list$Z = cbind(rep(1,n_admin2),cmat_admin2[,c('urb_frac')])
    data_list$bym2_var = 0
    
    start.time <- Sys.time()
    fit2 <- mod2$sample(
      data = data_list,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.99,
      show_messages = F,
      show_exceptions = F,
      refresh=00)
    
    time_log <- rbind(time_log,
                      data.frame(setting = setting, model = 'Mod2a', time = Sys.time() - start.time))
    
    
    print(Sys.time() - start.time0)
  }
}

save(time_log, file = '/Users/alanamcgovern/Desktop/Research/FHVariance_Smoothing/Simulations/time_log.rda')

time_log %>% filter(model == 'Mod1') %>% group_by(setting == 4) %>% reframe(mean(time))
time_log %>% filter(model == 'Mod1a') %>% group_by(setting == 4) %>% reframe(mean(time))
time_log %>% filter(model == 'Mod2') %>% group_by(setting == 4) %>% reframe(mean(time))
time_log %>% filter(model == 'Mod2a') %>% group_by(setting == 4) %>% reframe(mean(time))
