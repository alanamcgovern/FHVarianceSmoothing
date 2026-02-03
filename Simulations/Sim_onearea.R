setwd('/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing')

setting <- 3
mod <- cmdstan_model(stan_file = "Stan/NonStrat_onearea.stan")

setwd(paste0('Simulations/Sim',setting))

load(file='objects.rda')
load(file='sampled_clusters.rda')
load(file='direct.rda')
load(file='params.rda')

area <- 4
# check sample sizes
hist(unlist(lapply(sampled_clusters,function(x){nrow(x[x$admin2==area,])})))

true_theta = params$admin2_mean[area]
true_logsig2 = log(params$sig2[area])
true_v = var(unlist(lapply(direct,function(x){x$mean[area]})))

admin.key <- objects$admin.key

theta_est <- logsig2_est <- v_est <- NULL
for(k in 1:100){
  cat(k,'\n')
  dir.dat <- direct[[k]]
  
  change_id <- which(dir.dat$variance < 1e-10)
  if(length(change_id)>0){
    dir.dat[change_id,]$mean <- NA
  }
  
  if(!is.na(dir.dat$mean[area])){
    # get input ----
    tmp_D <- sampled_clusters[[k]] %>% filter(admin2 == area)
    N_D <- c(sum(1-tmp_D$urban),sum(tmp_D$urban))
    
    df_naive <- pmax(1,sum(N_D) - sum(N_D>0))
    v_hat_scaled_naive <- dir.dat$variance[area]*pmax(1,sum(N_D) - sum(N_D>0))
    Cons_naive <- 1/(sum(N_D)*mean(tmp_D$n))
   scale_naive <- 1/(sum(N_D)*mean(tmp_D$n)*pmax(1,sum(N_D) - sum(N_D>0)))
    
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
    
    keep_id <- which(out$values > 1e-10)
    q <- out$values[keep_id]
    df_exact <- (sum(q))^2/sum(q^2)
    
    v_hat_scaled_exact <- dir.dat$variance[area]*as.numeric(t(omega)%*%S%*%omega)*sum(q)/sum(q^2)
    cons_exact <- as.numeric((t(omega)%*%S%*%omega)/sum(omega)^2)
    
    # naive -----
    data_list = list(y=dir.dat$mean[area],
                     v_hat_scaled = v_hat_scaled_naive,
                     df=df_naive,
                     Cons = Cons_naive)

    fit1 <- mod$sample(
      data = data_list,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      show_messages = F,
      show_exceptions = F,
      refresh=00)

    naive_param_fit <- as.matrix(unclass(fit1$draws(format = "matrix")))

    # satterwhaite ----
    data_list = list(y=dir.dat$mean[area],
                     v_hat_scaled = v_hat_scaled_exact,
                     Cons = cons_exact,
                     df = df_exact)

    fit2 <-  mod$sample(
      data = data_list,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      show_messages = F,
      show_exceptions = F,
      refresh=0)

    satt_param_fit <- as.matrix(unclass(fit2$draws(format = "matrix")))

    # record ----
    theta_est <- rbind(theta_est,
                       data.frame(sim = k,
                                  method = c('naive','satt'),
                                  mean = c(mean(naive_param_fit[,2]),mean(satt_param_fit[,2])),
                                  lower = c(quantile(naive_param_fit[,2],0.05),quantile(satt_param_fit[,2],0.05)),
                                  upper = c(quantile(naive_param_fit[,2],0.95),quantile(satt_param_fit[,2],0.95))))

    logsig2_est <- rbind(logsig2_est,
                         data.frame(sim = k,
                                    method = c('naive','satt'),
                                    mean = c(mean(naive_param_fit[,3]),mean(satt_param_fit[,3])),
                                    lower = c(quantile(naive_param_fit[,3],0.05),quantile(satt_param_fit[,3],0.05)),
                                    upper = c(quantile(naive_param_fit[,3],0.95),quantile(satt_param_fit[,3],0.95))))

    v_est <- rbind(v_est,
                         data.frame(sim = k,
                                    method = c('naive','satt'),
                                    mean = c(mean(naive_param_fit[,4]),mean(satt_param_fit[,4])),
                                    lower = c(quantile(naive_param_fit[,4],0.05),quantile(satt_param_fit[,4],0.05)),
                                    upper = c(quantile(naive_param_fit[,4],0.95),quantile(satt_param_fit[,4],0.95)),
                                    emp = cons_exact*exp(true_logsig2)))
    
  }
}

theta_est %>% ggplot() + geom_histogram(aes(mean)) + facet_grid(~method) +geom_vline(xintercept = true_theta,col='blue')
logsig2_est %>% ggplot() + geom_histogram(aes(mean)) + facet_grid(~method) +geom_vline(xintercept = true_logsig2,col='blue')
v_est %>% ggplot() + geom_point(aes(emp,mean))+ facet_grid(~method) +geom_abline(intercept = 0,slope=1,col='red')

plot(v_est[v_est$method=='naive',]$mean,v_est[v_est$method=='satt',]$mean)
abline(0,1)

theta_est %>% mutate(cov= lower<true_theta & upper > true_theta,
                     width = upper - lower) %>% group_by(method) %>% summarise(mean(cov),mean(width))

