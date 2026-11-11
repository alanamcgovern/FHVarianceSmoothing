
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

mod <- cmdstan_model(stan_file = "Standard.stan")
##mod <- cmdstan_model(stan_file = "Stan/Standard.stan")

setwd(paste0('Sim',setting))

load(file='objects.rda')
load(file='cmat_admin2.rda')
load(file='sampled_clusters.rda')
load(file='direct.rda')

folders <- c('Mod0')
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

dir.dat <- direct[[k]]
dir.dat <- merge(dir.dat,cmat_admin2)

change_id <- which(dir.dat$variance < 1e-5)
if(length(change_id)>0){
  dir.dat[change_id,]$variance <- NA
  dir.dat[change_id,]$mean <- NA
}

# get input 
data_areas <- which(!is.na(dir.dat$mean))

data_list = list(m=n_admin2,
                 m_data=length(data_areas),
                 data_areas = data_areas,
                 y=dir.dat$mean[data_areas],
                 v_hat = dir.dat$variance[data_areas],
                 p_mean = ncol(mean_model_matrix_admin2),
                 X = mean_model_matrix_admin2,
                 bym2_mean = 1,
                 # all bym2 stuff
                 N_edges = nrow(objects$nodes2),
                 node_1 = objects$nodes2$node1,
                 node_2 = objects$nodes2$node2,
                 car_scale = objects$car_scale2) 

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
  refresh=0)

theta_fit1 <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod0',
                     mean = apply(theta_fit1,2,mean),
                     upper = apply(theta_fit1,2,quantile,probs=0.95),
                     lower = apply(theta_fit1,2,quantile,probs=0.05)),file=paste0('Mod0/Result',k,'.csv'))

hyper_fit1 <- as.matrix(unclass(fit1$draws(variables = c("beta",'sig_u','phi1'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod0',
                     param = colnames(hyper_fit1),
                     mean = colMeans(hyper_fit1),
                     median = apply(hyper_fit1,2,median)),
          file=paste0('Mod0/Hyper/Hyper_Result',k,'.csv'))

# diagnostics
summary_df <- fit1$summary()
diag_df    <- fit1$diagnostic_summary()

write.csv(data.frame(sim=k,
                     model = 'Mod0',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail),
          file=paste0('Mod0/Diagnostic/Diag_Summary_',k,'.csv'))

prob_params <- summary_df %>% filter((!is.na(rhat) & rhat > 1.01) | ess_bulk < 400 | ess_tail < 100)
if(nrow(prob_params)>0){
  write.csv(prob_params,file=paste0('Mod0/Diagnostic/Diag_detect_',k,'.csv'))
}


rm(fit1,theta_fit1,hyper_fit1)



