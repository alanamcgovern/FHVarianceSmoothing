suppressMessages({
  library(dplyr)
  library(cmdstanr)
  library(Matrix)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly=TRUE)

setting <- as.numeric(args[1])
k <- as.numeric(args[2])

mod <- cmdstan_model(stan_file = "CSS.stan")
##mod <- cmdstan_model(stan_file = "Stan/CSS.stan")

setwd(paste0('Sim',setting))

load(file='objects.rda')
load(file='cmat_admin2.rda')
load(file='sampled_clusters.rda')
load(file='direct.rda')

folders <- c('Mod3','Mod3a')
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

# get input 
data_areas <- which(!is.na(dir.dat$mean))
Cons <- scale <- df <- v_hat_scaled <- rep(NA,n_admin2)

for(area in data_areas){
  
  n_D <- as.numeric(sampled_clusters[[k]] %>% filter(admin2 == area) %>% summarise(sum(n)))
  
  df[area] <- n_D - 1
  v_hat_scaled[area] <- dir.dat$variance[area]*(n_D - 1)
  Cons[area] <- 1/n_D
  scale[area] <- 1/(n_D*(n_D - 1))
}

data_areas <- which(!is.na(df))

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
  adapt_delta = 0.99,
  max_treedepth = 15,
  show_messages = F,
  show_exceptions = F,
  refresh=00)

theta_fit1 <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod3',
                     mean = apply(theta_fit1,2,mean),
                     upper = apply(theta_fit1,2,quantile,probs=0.95),
                     lower = apply(theta_fit1,2,quantile,probs=0.05)),file=paste0('Mod3/Result',k,'.csv'))

v_fit1 <- (as.matrix(unclass(fit1$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod3',
                     mean = apply(v_fit1,2,mean),
                     median = apply(v_fit1,2,median)),file=paste0('Mod3/Hyper/Var_Result',k,'.csv'))

hyper_fit1 <- as.matrix(unclass(fit1$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod3',
                     param = colnames(hyper_fit1),
                     mean = colMeans(hyper_fit1),
                     median = apply(hyper_fit1,2,median)),
          file=paste0('Mod3/Hyper/Hyper_Result',k,'.csv'))

params_fit1 <- data.frame(sim=k,
                          model='Mod3',
                          area = 1:n_admin2,
                          scale = exp(colMeans(hyper_fit1)[str_detect(colnames(hyper_fit1),'log_sig2')])*scale,
                          df = df)
write.csv(params_fit1,
          file=paste0('Mod3/Hyper/Chisq_params',k,'.csv'))

# diagnostics
summary_df <- fit1$summary()
diag_df    <- fit1$diagnostic_summary()

write.csv(data.frame(sim=k,
                     model = 'Mod3',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail),
          file=paste0('Mod3/Diagnostic/Diag_Summary_',k,'.csv'))

prob_params <- summary_df %>% filter((!is.na(rhat) & rhat > 1.01) | ess_bulk < 400 | ess_tail < 100)
if(nrow(prob_params)>0){
  write.csv(prob_params,file=paste0('Mod3/Diagnostic/Diag_detect_',k,'.csv'))
}

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
  max_treedepth = 15,
  adapt_delta = 0.99,
  show_messages = F,
  show_exceptions = F,
  refresh=00)

theta_fit1 <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod3a',
                     mean = apply(theta_fit1,2,mean),
                     upper = apply(theta_fit1,2,quantile,probs=0.95),
                     lower = apply(theta_fit1,2,quantile,probs=0.05)),file=paste0('Mod3a/Result',k,'.csv'))

v_fit1 <- (as.matrix(unclass(fit1$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod3a',
                     mean = apply(v_fit1,2,mean),
                     median = apply(v_fit1,2,median)),file=paste0('Mod3a/Hyper/Var_Result',k,'.csv'))

hyper_fit1 <- as.matrix(unclass(fit1$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod3a',
                     param = colnames(hyper_fit1),
                     mean = colMeans(hyper_fit1),
                     median = apply(hyper_fit1,2,median)),
          file=paste0('Mod3a/Hyper/Hyper_Result',k,'.csv'))

params_fit1 <- data.frame(sim=k,
                          model='Mod3a',
                          area = 1:n_admin2,
                          scale = exp(colMeans(hyper_fit1)[str_detect(colnames(hyper_fit1),'log_sig2')])*scale,
                          df = df)
write.csv(params_fit1,
          file=paste0('Mod3a/Hyper/Chisq_params',k,'.csv'))

# diagnostics
summary_df <- fit1$summary()
diag_df    <- fit1$diagnostic_summary()

write.csv(data.frame(sim=k,
                     model = 'Mod3a',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail),
          file=paste0('Mod3a/Diagnostic/Diag_Summary_',k,'.csv'))

prob_params <- summary_df %>% filter((!is.na(rhat) & rhat > 1.01) | ess_bulk < 400 | ess_tail < 100)
if(nrow(prob_params)>0){
  write.csv(prob_params,file=paste0('Mod3a/Diagnostic/Diag_detect_',k,'.csv'))
}


rm(fit1,theta_fit1,v_fit1,hyper_fit1)


