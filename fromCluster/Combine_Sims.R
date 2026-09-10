library(readr)
library(stringr)


models <- c('Mod1','Mod1a','Mod2','Mod2a')#,'Mod4','Mod4a')#,'Mod3')

args <- commandArgs(trailingOnly=TRUE)
setting <- as.numeric(args[1])

setwd(paste0('Sim',setting))

results <- NULL
for(mod in models){
  tmp_list <- lapply(1:100,function(k){
    read.csv(file = paste0(mod,'/Result',k,'.csv'))
  })
  results <- rbind(results,do.call(rbind,tmp_list))
}

write_csv(results,file='Summary.csv')


results <- NULL
for(mod in models){
  tmp_list <- lapply(1:100,function(k){
    read.csv(file = paste0(mod,'/Hyper/Var_Result',k,'.csv'))
  })
  results <- rbind(results,do.call(rbind,tmp_list)[,c('sim','area','model','mean')])
}

write_csv(results,file='Summary_var.csv')

results <- NULL
for(mod in models){
  tmp_list <- lapply(1:100,function(k){
    read.csv(file = paste0(mod,'/Hyper/Chisq_params',k,'.csv'))
  })
  results <- rbind(results,do.call(rbind,tmp_list)[,c('sim','area','model','scale','df')])
}

write_csv(results,file='Summary_chisq_params.csv')


results <- NULL
for(mod in models){
  tmp_list <- lapply(1:100,function(k){
    read.csv(file = paste0(mod,'/Hyper/Hyper_Result',k,'.csv'))
  })
  results <- rbind(results,do.call(rbind,tmp_list))
}

write_csv(results,file='Summary_hyper.csv')

results <- NULL
for(mod in models){
  tmp_list <- lapply(1:100,function(k){
    read.csv(file = paste0(mod,'/Diagnostic/Diag_Summary_',k,'.csv'))
  })
  results <- rbind(results,do.call(rbind,tmp_list))
}

write_csv(results,file='Summary_diag.csv')



