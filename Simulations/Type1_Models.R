## fit Type 1 model (naive): 
# does not account for sample weights and sizes or unplanned domain
# only account for stratification in the mean and variance by including the urban proportion as a covariate in the latent layer

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

mod <- cmdstan_model(stan_file = "NonStrat.stan")
##mod <- cmdstan_model(stan_file = "Stan/NonStrat.stan")

setwd(paste0('Sim',setting))

load(file='objects.rda')
load(file='cmat_admin2.rda')
load(file='sampled_clusters.rda')
load(file='direct.rda')

folders <- c('Mod1','Mod1a')
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

change_id <- which(dir.dat$variance < 1e-10)
if(length(change_id)>0){
  dir.dat[change_id,]$variance <- NA
  dir.dat[change_id,]$mean <- NA
}

# get input 
data_areas <- which(!is.na(dir.dat$mean))
Cons_naive <- scale_naive <- df_naive <- v_hat_scaled_naive <- rep(NA,n_admin2)

for(area in data_areas){

  tmp_D <- sampled_clusters[[k]] %>% filter(admin2 == area)
  N_D <- c(sum(1-tmp_D$urban),sum(tmp_D$urban))
  
 # if(sum(N_D>0)==2){
    df_naive[area] <- pmax(1,sum(N_D) - sum(N_D>0))
    v_hat_scaled_naive[area] <- dir.dat$variance[area]*pmax(1,sum(N_D) - sum(N_D>0))
    Cons_naive[area] <- 1/(sum(N_D)*mean(tmp_D$n))
    scale_naive[area] <- 1/(sum(N_D)*mean(tmp_D$n)*pmax(1,sum(N_D) - sum(N_D>0)))
 # }
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

fit1 <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
  show_messages = F,
  show_exceptions = F,
  refresh=0)

theta_fit1 <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod1',
                     mean = apply(theta_fit1,2,mean),
                     upper = apply(theta_fit1,2,quantile,probs=0.95),
                     lower = apply(theta_fit1,2,quantile,probs=0.05)),file=paste0('Mod1/Result',k,'.csv'))

v_fit1 <- (as.matrix(unclass(fit1$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod1',
                     mean = apply(v_fit1,2,mean)),file=paste0('Mod1/Hyper/Var_Result',k,'.csv'))

hyper_fit1 <- as.matrix(unclass(fit1$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod1',
                     param = colnames(hyper_fit1),
                     mean = colMeans(hyper_fit1)),
          file=paste0('Mod1/Hyper/Hyper_Result',k,'.csv'))

params_fit1 <- data.frame(sim=k,
                          model='Mod1',
                          area = 1:n_admin2,
                          scale = exp(colMeans(hyper_fit1)[str_detect(colnames(hyper_fit1),'log_sig2')])*scale_naive,
                          df = df_naive)
write.csv(params_fit1,
          file=paste0('Mod1/Hyper/Chisq_params',k,'.csv'))

# diagnostics
summary_df <- fit1$summary()
diag_df    <- fit1$diagnostic_summary()

write.csv(data.frame(sim=k,
                     model = 'Mod1',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail),
          file=paste0('Mod1/Diagnostic/Diag_Summary_',k,'.csv'))

prob_params <- summary_df %>% filter((!is.na(rhat) & rhat > 1.01) | ess_bulk < 400 | ess_tail < 100)
if(nrow(prob_params)>0){
  write.csv(prob_params,file=paste0('Mod1/Diagnostic/Diag_detect_',k,'.csv'))
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
  adapt_delta = 0.99,
  show_messages = F,
  show_exceptions = F,
  refresh=0)

theta_fit1 <- as.matrix(unclass(fit1$draws(variables = c("theta"), format = "matrix")))
write.csv(data.frame(sim=k,
                     area = 1:n_admin2,
                     model = 'Mod1a',
                     mean = apply(theta_fit1,2,mean),
                     upper = apply(theta_fit1,2,quantile,probs=0.95),
                     lower = apply(theta_fit1,2,quantile,probs=0.05)),file=paste0('Mod1a/Result',k,'.csv'))

v_fit1 <- (as.matrix(unclass(fit1$draws(variables = c("v_data"), format = "matrix"))))
write.csv(data.frame(sim=k,
                     area = data_areas,
                     model = 'Mod1a',
                     mean = apply(v_fit1,2,mean)),file=paste0('Mod1a/Hyper/Var_Result',k,'.csv'))

hyper_fit1 <- as.matrix(unclass(fit1$draws(variables = c("log_sig2","beta",'gamma','sig_u','phi1','phi2'), format = "matrix")))

write.csv(data.frame(sim=k,
                     model = 'Mod1a',
                     param = colnames(hyper_fit1),
                     mean = colMeans(hyper_fit1)),
          file=paste0('Mod1a/Hyper/Hyper_Result',k,'.csv'))

params_fit1 <- data.frame(sim=k,
                          model='Mod1a',
                          area = 1:n_admin2,
                          scale = exp(colMeans(hyper_fit1)[str_detect(colnames(hyper_fit1),'log_sig2')])*scale_naive,
                          df = df_naive)
write.csv(params_fit1,
          file=paste0('Mod1a/Hyper/Chisq_params',k,'.csv'))

# diagnostics
summary_df <- fit1$summary()
diag_df    <- fit1$diagnostic_summary()

write.csv(data.frame(sim=k,
                     model = 'Mod1a',
                     divergences = sum(diag_df$num_divergent),
                     max_treedepth = sum(diag_df$num_max_treedepth),
                     min_efbmi = min(diag_df$ebfmi),
                     max_rhat = max(summary_df$rhat,na.rm=T),
                     bulk_ess = summary_df[summary_df$variable=='lp__',]$ess_bulk,
                     tail_ess = summary_df[summary_df$variable=='lp__',]$ess_tail),
          file=paste0('Mod1a/Diagnostic/Diag_Summary_',k,'.csv'))

prob_params <- summary_df %>% filter((!is.na(rhat) & rhat > 1.01) | ess_bulk < 400 | ess_tail < 100)
if(nrow(prob_params)>0){
  write.csv(prob_params,file=paste0('Mod1a/Diagnostic/Diag_detect_',k,'.csv'))
}


rm(fit1,theta_fit1,v_fit1,hyper_fit1)



