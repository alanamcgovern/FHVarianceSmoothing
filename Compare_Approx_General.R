

H <- 2
strata_totals <- c(100*100,100*100)

## new ----------
set.seed(12626)
for(setting in c('V1','V2','V3')){
  
  if(setting == 'V1'){
    N <- c(5,5)
    prob_D <- 1 # proportion in target domain
  }else if(setting=='V2'){
    N <- c(15,15)
    prob_D <- 0.2
  }else if(setting == 'V3'){
    N <- c(2,5)
    prob_D <- 1 # proportion in target domain
  }else{ 
    stop('error')
  }
  
  dat <- data.frame(h = unlist(lapply(1:H,function(h){rep(h,N[h])})))
  dat$D <- rbinom(nrow(dat),1,prob_D)
  
  dat$mu <- ifelse(dat$D==1 & dat$h==1,-1,
                   ifelse(dat$D==1 & dat$h==2,1,2))
  
  mu_diff <- 2
  #mu <- c(rep(-1,N[1]),rep(1,N[2]))
  sig2_k <- 2
  
  if(setting %in% c('V1','V2')){
    dat$n <- 10
  }else if(setting %in% c('V3')){
    dat$n <- rpois(nrow(dat),10)
  }else{
    stop('error')
  }
  
  dat$wt <- strata_totals[dat$h]/(N[dat$h]*dat$n)/1e3
  
  S <- diag(1/dat$n)
  
  sig2_c <- sig2_k/mean(dat$n)
  
  # get quantities-----
  
  if(sum(1-dat$D)>0)
    dat[dat$D==0,]$wt <- 0
  omega <- dat$wt*dat$n
  
  B_comps <- lapply(1:H,function(h){
    if(N[h]>0){
      return(N[h]/(N[h]-1)*(diag(1,N[h]) - 1/N[h]*matrix(1,N[h],N[h])))
    }
  })
  
  which.not.null <- c(1:H)[!sapply(B_comps, is.null)]
  B <- bdiag(B_comps[which.not.null])
  W <- diag(omega)%*%(diag(1,sum(N)) - (rep(1,sum(N))%*%t(omega))/sum(omega))
  M <- t(W)%*%B%*%W
  
  L <- chol(S)              # S = L' L
  Li <- backsolve(L, diag(nrow(S)))  # L^{-1}
  
  # Form scale-free matrix A = S^{1/2} C S^{1/2}
  A <- t(L) %*% M %*% L
  
  out <- eigen(A, symmetric = TRUE)
  V <- Li %*% out$vectors
  
  keep_id <- which(out$values > 1e-10)
  q <- out$values[keep_id]
  r <- length(keep_id)
  
  V <- V[,keep_id]
  nu <- t(V)%*%(dat$h-1)
  
  delta <- (nu*mu_diff)^2/sig2_k
  
  ## calculate exact -----
  dist_exact <- replicate(2e4,{
    sig2_k/(sum(omega)^2)*sum(q*sapply(1:r,function(i){rchisq(1,1,ncp =  delta[i])}))
  })
  
  q_exact <- quantile(dist_exact,probs=seq(0.01,0.99,0.01))
  
  ## compare to Satterthwhaite -----
  
  Q1 <- sum(q*(1+delta))
  Q2 <- sum(q^2*(1+2*delta))
  
  scale_satt <- sig2_k/sum(omega)^2*Q2/Q1
  df_satt <- Q1^2/Q2
  
  dist_satt <- scale_satt*rchisq(2e4,df_satt)
  q_satt <- quantile(dist_satt,seq(0.01,0.99,0.01))
  
  mean_satt <- scale_satt*df_satt
  var_satt <- scale_satt^2*2*df_satt
  skew_satt <- sqrt(8/df_satt)
  
  ## compare to debiased Satterthwhaite -----
  
  scale_satt_debias <- as.numeric(sig2_k/sum(omega)^2*(t(omega)%*%S%*%omega)*Q2/Q1^2)
  mean_satt_debias <- scale_satt_debias*df_satt
  var_satt_debias <- scale_satt_debias^2*2*df_satt
  
  dist_satt_debias <- scale_satt_debias*rchisq(2e4,df_satt)
  q_satt_debias <- quantile(dist_satt_debias,seq(0.01,0.99,0.01))
  
  
  ## compare to naive -----
  
  df_naive <- sum(dat$D) - length(unique(dat[dat$D==1,]$h))
  scale_naive <- sig2_c/(sum(dat$D)*df_naive)
  dist_naive <- sig2_c/(sum(dat$D)*df_naive)*rchisq(2e4,df_naive)
  
  q_naive <- quantile(dist_naive,probs=seq(0.01,0.99,0.01))
  
  mean_naive <- scale_naive*df_naive
  var_naive <- scale_naive^2*2*df_naive
  skew_naive <- sqrt(8/df_naive)
  
  ## plot -----
  quants <- data.frame(exact = q_exact,
                       Naive = q_naive,
                       Satterthwhaite = q_satt,
                       Debiased_Satterthwhaite = q_satt_debias) %>%
    pivot_longer(
      cols = c(Naive, Satterthwhaite,Debiased_Satterthwhaite),
      names_to = "approx",
      values_to = "value"
    )
  quants[quants$approx=='Debiased_Satterthwhaite',]$approx <- 'Debiased Satterthwhaite'
  
  params <- data.frame(
    approx = c("Naive", "Satterthwhaite","Debiased Satterthwhaite"),
    label  = c(paste0('scale = ',round(scale_naive,5),'<br> df = ',df_naive,'<br> mean = ',round(mean_naive,3),'<br> var = ',round(var_naive,6)), 
               paste0('scale = ',round(scale_satt,5),'<br> df = ',round(df_satt,2),'<br> mean = ',round(mean_satt,3),'<br> var = ',round(var_satt,6)),
               paste0('scale = ',round(scale_satt_debias,5),'<br> df = ',round(df_satt,2),'<br> mean = ',round(mean_satt_debias,3),'<br> var = ',round(var_satt_debias,6)))
  )
  
  lims <- range(quants$exact,quants$value)
  
  g_tmp <- ggplot(quants, aes(x = exact, y = value)) +
    geom_point(size=0.75) + geom_abline(intercept = 0,slope=1,lwd=0.5,col='gray40') +
    facet_wrap(~ factor(approx,
                        levels = c("Naive", "Satterthwhaite","Debiased Satterthwhaite")), 
               ncol = 3) +
    xlim(lims) + ylim(lims) +
    ggtext::geom_richtext(
      data = params,
      size=3,
      aes(label = label),
      inherit.aes = FALSE,
      x = -Inf, y = Inf,
      hjust = -0.05, 
      vjust = 1.1
    ) + theme_bw() + theme(strip.text = element_text(size=12)) +
    labs(x = "Exact", y = 'Approximation') 
  
  if(setting == 'V1'){
    g0 <- g_tmp + ggtitle('A. All design assumptions met')
  }else if(setting=='V2'){
    g1 <- g_tmp + ggtitle('B. Unplanned domain')
  }else if(setting=='V3'){
    g2 <- g_tmp + ggtitle('C. Unequal sample size')
  }else{stop('error')
  }
  
}

ggarrange(plotlist = list(g0,g1,g2),nrow=3)




## correct specification and estimation -----
N <- c(5,5)
prob_D <- 1 # proportion in target domain

mu <- c(rep(-1,N[1]),rep(1,N[2]))
sig2_k <- 2 #individual level variance

dat <- data.frame(h = unlist(lapply(1:H,function(h){rep(h,N[h])})),
                  n = rep(10,sum(N)))
dat$wt <- strata_totals[dat$h]/(N[dat$h]*dat$n)/1e3

#dat <- as.data.table(dat)

S <- diag(1/dat$n)
sig2_c <- sig2_k/mean(dat$n)

# get quantities-

dat$D <- rbinom(nrow(dat),1,prob_D)

if(sum(1-dat$D)>0)
  dat[dat$D==0,]$wt <- 0
omega <- dat$wt*dat$n

B_comps <- lapply(1:H,function(h){
  if(N[h]>0){
    return(N[h]/(N[h]-1)*(diag(1,N[h]) - 1/N[h]*matrix(1,N[h],N[h])))
  }
})

which.not.null <- c(1:H)[!sapply(B_comps, is.null)]
B <- bdiag(B_comps[which.not.null])
W <- diag(omega)%*%(diag(1,sum(N)) - (rep(1,sum(N))%*%t(omega))/sum(omega))
M <- t(W)%*%B%*%W

L <- chol(S)              # S = L' L
Li <- backsolve(L, diag(nrow(S)))  # L^{-1}
A <- t(L) %*% M %*% L
out <- eigen(A, symmetric = TRUE)
u <- Li %*% out$vectors

keep_id <- which(out$values > 1e-10)
q <- out$values[keep_id]
r <- length(q)
u <- u[,keep_id]
delta <- (t(u)%*%mu)^2/sig2_k

## calculate exact 

dist_exact <- replicate(2e4,{
  sig2_k/(sum(omega)^2)*sum(q*sapply(1:r,function(i){rchisq(1,1,ncp =  delta[i])}))
})

q_exact <- quantile(dist_exact,probs=seq(0.01,0.99,0.01))

## compare to naive 

df_naive <- sum(dat$D) - length(unique(dat[dat$D==1,]$h))
scale_naive <- sig2_c/(sum(dat$D)*df_naive)
dist_naive <- scale_naive*rchisq(2e4,df_naive)

q_naive <- quantile(dist_naive,probs=seq(0.01,0.99,0.01))

## compare to Satterthwhaite 

Q1 <- sum(q*(1+delta))
Q2 <- sum(q^2*(1+2*delta))

scale_satt <- sig2_k/(sum(omega)^2)*Q2/Q1
df_satt <- Q1^2/Q2
dist_satt <- scale_satt*rchisq(2e4,df_satt)
q_satt <- quantile(dist_satt,seq(0.01,0.99,0.01))

## satt debiased
scale_satt_debias <- as.numeric(sig2_k/(sum(omega)^2)*(t(omega) %*%S %*% omega)*Q2/Q1^2)
dist_satt_debias <- scale_satt_debias*rchisq(2e4,df_satt)
q_satt_debias <- quantile(dist_satt_debias,seq(0.01,0.99,0.01))

## plot 
quants <- data.frame(exact = q_exact,
                     Naive = q_naive,
                     Satterthwhaite = q_satt,
                     Debiased_Satterthwhaite = q_satt_debias) %>%
  pivot_longer(
    cols = c(Naive, Satterthwhaite,Debiased_Satterthwhaite),
    names_to = "approx",
    values_to = "value"
  )

params <- data.frame(
  approx = c("Naive", "Satterthwhaite","Debiased_Satterthwhaite"),
  label  = c(paste0('scale = ',round(scale_naive,4),'\n df = ',df_naive), 
             paste0('scale = ',round(scale_satt,4),'\n df = ',round(df_satt,2)),
             paste0('scale = ',round(scale_satt_debias,4),'\n df = ',round(df_satt,2)))
)
lims <- range(quants$exact,quants$value)
g <- ggplot(quants, aes(x = exact, y = value)) +
  xlim(lims) + ylim(lims) +
  geom_point(size=0.75) + geom_abline(intercept = 0,slope=1,lwd=0.5,col='gray40') +
  facet_wrap(~ approx, ncol = 2) +
  geom_label(
    data = params,
    aes(label = label),
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE
  ) + theme_bw() + theme(strip.text = element_text(size=12)) +
  labs(x = "Exact", y = 'Approximation') 

print(g)

## misspecification -----
set.seed(12626)
for(setting in c('V1','V2','V3','V4')){

if(setting=='V1'){
  N <- c(15,15)
  prob_D <- 0.2
}else if(setting == 'V2'){
  N <- c(2,5)
  prob_D <- 1 # proportion in target domain
}else if(setting == 'V3'){
  N <- c(5,5)
  prob_D <- 1
}else if(setting == 'V4'){
   N <- c(15,20)
  prob_D <- 0.2 # proportion in target domain
}else{ 
  stop('error')
}
 
dat <- data.frame(h = unlist(lapply(1:H,function(h){rep(h,N[h])})))
dat$D <- rbinom(nrow(dat),1,prob_D)

dat$mu <- ifelse(dat$D==1 & dat$h==1,-1,
                 ifelse(dat$D==1 & dat$h==2,1,2))

mu_diff <- 2
#mu <- c(rep(-1,N[1]),rep(1,N[2]))
sig2_k <- 2

if(setting %in% c('V1','V3')){
  dat$n <- 10
}else if(setting %in% c('V2','V4')){
  dat$n <- rpois(nrow(dat),10)
}else{
  stop('error')
}

dat$wt <- strata_totals[dat$h]/(N[dat$h]*dat$n)/1e3

S <- diag(1/dat$n)

sig2_c <- sig2_k/mean(dat$n)

# get quantities-----

if(sum(1-dat$D)>0)
  dat[dat$D==0,]$wt <- 0
omega <- dat$wt*dat$n

B_comps <- lapply(1:H,function(h){
  if(N[h]>0){
    return(N[h]/(N[h]-1)*(diag(1,N[h]) - 1/N[h]*matrix(1,N[h],N[h])))
  }
})

which.not.null <- c(1:H)[!sapply(B_comps, is.null)]
B <- bdiag(B_comps[which.not.null])
W <- diag(omega)%*%(diag(1,sum(N)) - (rep(1,sum(N))%*%t(omega))/sum(omega))
M <- t(W)%*%B%*%W

L <- chol(S)              # S = L' L
Li <- backsolve(L, diag(nrow(S)))  # L^{-1}

# Form scale-free matrix A = S^{1/2} C S^{1/2}
A <- t(L) %*% M %*% L

out <- eigen(A, symmetric = TRUE)
V <- Li %*% out$vectors

keep_id <- which(out$values > 1e-10)
q <- out$values[keep_id]
r <- length(keep_id)

V <- V[,keep_id]
nu <- t(V)%*%(dat$h-1)

delta <- (nu*mu_diff)^2/sig2_k

# generate empirical -------

df_t <- 5
s = sqrt(sig2_k/dat$n*(df_t-2)/df_t)

dist_emp <- replicate(2e4,{
  if(setting %in% c('V1','V2')){
    y <- as.vector(Rfast::rmvnorm(1,dat$mu,sig2_k*S))
  }else if(setting %in% c('V3','V4')){
    y <- dat$mu + s*rt(sum(N),df_t)
  }else{
    return(NA)
  }
  return(as.numeric(t(y)%*%M%*%y/sum(omega)^2))
})

q_emp <- quantile(dist_emp,probs=seq(0.01,0.99,0.01))

## calculate exact -----
dist_exact <- replicate(2e4,{
  sig2_k/(sum(omega)^2)*sum(q*sapply(1:r,function(i){rchisq(1,1,ncp =  delta[i])}))
})

q_exact <- quantile(dist_exact,probs=seq(0.01,0.99,0.01))

## compare to Satterthwhaite -----

Q1 <- sum(q*(1+delta))
Q2 <- sum(q^2*(1+2*delta))

scale_satt <- sig2_k/sum(omega)^2*Q2/Q1
df_satt <- Q1^2/Q2

dist_satt <- scale_satt*rchisq(2e4,df_satt)
q_satt <- quantile(dist_satt,seq(0.01,0.99,0.01))

mean_satt <- scale_satt*df_satt
var_satt <- scale_satt^2*2*df_satt
skew_satt <- sqrt(8/df_satt)

## compare to debiased Satterthwhaite -----

scale_satt_debias <- as.numeric(sig2_k/sum(omega)^2*(t(omega)%*%S%*%omega)*Q2/Q1^2)
mean_satt_debias <- scale_satt_debias*df_satt
var_satt_debias <- scale_satt_debias^2*2*df_satt

dist_satt_debias <- scale_satt_debias*rchisq(2e4,df_satt)
q_satt_debias <- quantile(dist_satt_debias,seq(0.01,0.99,0.01))


## compare to naive -----

df_naive <- sum(dat$D) - length(unique(dat[dat$D==1,]$h))
scale_naive <- sig2_c/(sum(dat$D)*df_naive)
dist_naive <- sig2_c/(sum(dat$D)*df_naive)*rchisq(2e4,df_naive)

q_naive <- quantile(dist_naive,probs=seq(0.01,0.99,0.01))

mean_naive <- scale_naive*df_naive
var_naive <- scale_naive^2*2*df_naive
skew_naive <- sqrt(8/df_naive)

## plot -----
quants <- data.frame(exact = q_exact,
                     Naive = q_naive,
                     Satterthwhaite = q_satt,
                     Debiased_Satterthwhaite = q_satt_debias) %>%
  pivot_longer(
    cols = c(Naive, Satterthwhaite,Debiased_Satterthwhaite),
    names_to = "approx",
    values_to = "value"
  )

params <- data.frame(
  approx = c("Naive", "Satterthwhaite","Debiased_Satterthwhaite"),
  label  = c(paste0('scale = ',round(scale_naive,5),'<br> df = ',df_naive,'<br> mean = ',round(mean_naive,3),'<br> var = ',round(var_naive,6)), 
             paste0('scale = ',round(scale_satt,5),'<br> df = ',round(df_satt,2),'<br> mean = ',round(mean_satt,3),'<br> var = ',round(var_satt,6)),
             paste0('scale = ',round(scale_satt_debias,5),'<br> df = ',round(df_satt,2),'<br> mean = ',round(mean_satt_debias,3),'<br> var = ',round(var_satt_debias,6)))
)

lims <- range(quants$exact,quants$value)

g_tmp <- ggplot(quants, aes(x = exact, y = value)) +
  geom_point(size=0.75) + geom_abline(intercept = 0,slope=1,lwd=0.5,col='gray40') +
  facet_wrap(~ approx, ncol = 3) +
  xlim(lims) + ylim(lims) +
  ggtext::geom_richtext(
    data = params,
    size=3,
    aes(label = label),
    inherit.aes = FALSE,
    x = -Inf, y = Inf,
    hjust = -0.05, 
   vjust = 1.1
  ) + theme_bw() + theme(strip.text = element_text(size=12)) +
  labs(x = "Exact", y = 'Approximation') 

if(setting=='V1'){
  g1 <- g_tmp + ggtitle('A. Unplanned domain')
}else if(setting=='V2'){
  g2 <- g_tmp + ggtitle('B. Unequal sample size')
}else if(setting=='V4'){
  g4 <- g_tmp + ggtitle('D. Unplanned domain, unequal sample size, & non-Gaussian observations')
}else if(setting=='V3'){
  g3 <- g_tmp + ggtitle('C. Non-Gaussian observations, t-distributed with df=5')
}else{stop('error')
}

}

ggarrange(plotlist = list(g1,g2,g3,g4),nrow=4)


## wrong mean estimation ------
set.seed(12726)
  N <- c(15,20)
  prob_D <- 0.2 # proportion in target domain
  
#  N <- c(5,15)
#  prob_D <- 1# proportion in target domain
  
  mu <- c(rep(-1,N[1]),rep(1,N[2]))
  sig2_k <- 2
  
  dat <- data.frame(h = unlist(lapply(1:H,function(h){rep(h,N[h])})),
                    n = rpois(sum(N),10)) 
  dat$wt <- 10*strata_totals[dat$h]/(N[dat$h]*dat$n)/1e3
  
  dat <- as.data.table(dat)
  
  S <- diag(1/dat$n)
  sig2_c <- sig2_k/mean(dat$n)
  
  # get quantities-----
  
  dat$D <- rbinom(nrow(dat),1,prob_D)
  
  if(sum(1-dat$D)>0)
    dat[dat$D==0,]$wt <- 0
  omega <- dat$wt*dat$n
  
  B_comps <- lapply(1:H,function(h){
    if(N[h]>0){
      return(N[h]/(N[h]-1)*(diag(1,N[h]) - 1/N[h]*matrix(1,N[h],N[h])))
    }
  })
  
  which.not.null <- c(1:H)[!sapply(B_comps, is.null)]
  B <- bdiag(B_comps[which.not.null])
  W <- diag(omega)%*%(diag(1,sum(N)) - (rep(1,sum(N))%*%t(omega))/sum(omega))
  M <- t(W)%*%B%*%W
  
  L <- chol(S)              # S = L' L
  Li <- backsolve(L, diag(nrow(S)))  # L^{-1}
  A <- t(L) %*% M %*% L
  out <- eigen(A, symmetric = TRUE)
  u <- Li %*% out$vectors
  
  keep_id <- which(out$values > 1e-10)
  q <- out$values[keep_id]
  r <- length(q)
  u <- u[,keep_id]
  delta <- (t(u)%*%mu)^2/sig2_k
  
  ## calculate dists under correct means -----
  # exact
  dist_exact <- replicate(2e4,{
    sig2_k/(sum(omega)^2)*sum(q*sapply(1:r,function(i){rchisq(1,1,ncp =  delta[i])}))
  })
  
  q_exact <- quantile(dist_exact,probs=seq(0.01,0.99,0.01))
  
  #Satterthwhaite 
  Q1 <- sum(q*(1+delta))
  Q2 <- sum(q^2*(1+2*delta))
  
  scale_correct <- sig2_k/sum(omega)^2*Q2/Q1
  df_correct <- Q1^2/Q2
  mean_correct <- scale_correct*df_correct
  var_correct <- 2*scale_correct^2*df_correct
  
  dist_satt <- scale_correct*rchisq(2e4,df_correct)
  q_satt_correct <- quantile(dist_satt,seq(0.01,0.99,0.01))
  
  ## naive 
  
  df_naive <- sum(dat$D) - length(unique(dat[dat$D==1,]$h))
  scale_naive <- sig2_c/(sum(dat$D)*df_naive)
  mean_naive <- df_naive*scale_naive
  var_naive <- 2*scale_naive^2*df_naive
  
  dist_naive <- scale_naive*rchisq(2e4,df_naive)
  q_naive <- quantile(dist_naive,probs=seq(0.01,0.99,0.01))
  
  ## calculate dists under incorrect means ------
  q_satt <- list() 
  scale_satt <- df_satt <- mean_satt <- var_satt <- NULL
  for(setting in c('M1','M2','M3')){
  ## add estimation error ----
  if(setting == 'M1'){
    mu <- c(rep(-0.5,N[1]),rep(0.5,N[2]))
  }else if(setting == 'M2'){
    mu <- c(rep(-1.5,N[1]),rep(1.5,N[2]))
  }else if(setting == 'M3'){
    mu <- c(rep(0,N[1]),rep(2,N[2]))
  }else{
    stop('error')
  }
  
  delta <- (t(u)%*%mu)^2/sig2_k
  
  ## compare to Satterthwhaite -----
  
  Q1 <- sum(q*(1+delta))
  Q2 <- sum(q^2*(1+2*delta))
  
  dist_satt <- sig2_k/sum(omega)^2*Q2/Q1*rchisq(2e4,Q1^2/Q2)
  q_satt[[length(q_satt)+1]] <- quantile(dist_satt,seq(0.01,0.99,0.01))
  scale_satt[length(scale_satt)+1] <- sig2_k/sum(omega)^2*Q2/Q1
  df_satt[length(df_satt)+1] <- Q1^2/Q2
  mean_satt[length(mean_satt)+1] <- Q1^2/Q2*sig2_k/sum(omega)^2*Q2/Q1
  var_satt[length(var_satt)+1] <- 2*(sig2_k/sum(omega)^2*Q2/Q1)^2*Q1^2/Q2
  
  }
  
  ## plot -----
  quants <- data.frame(exact = q_exact,
                       Naive = q_naive,
                       Satt_M1 = q_satt[[1]],
                       Satt_M2 = q_satt[[2]],
                #       Satt_M3 = q_satt[[3]],
                       Satt_correct = q_satt_correct) %>%
    pivot_longer(
      cols = c(Naive, Satt_M1,Satt_M2,Satt_correct),
      names_to = "approx",
      values_to = "value"
    )
  
  params <- data.frame(
    approx = c("Naive", "Satt_correct","Satt_M1","Satt_M2"),
    label  = c(paste0('scale = ',round(scale_naive,5),'<br> df = ',df_naive,'<br> mean = ',round(mean_naive,3),'<br> var = ',round(var_naive,6)), 
               paste0('scale = ',round(scale_correct,5),'<br> df = ',round(df_correct,2),'<br> mean = ',round(mean_correct,3),'<br> var = ',round(var_correct,6)),
               paste0('scale = ',round(scale_satt[1],5),'<br> df = ',round(df_satt[1],2),'<br> mean = ',round(mean_satt[1],3),'<br> var = ',round(var_satt[1],6)),
               paste0('scale = ',round(scale_satt[2],5),'<br> df = ',round(df_satt[2],2),'<br> mean = ',round(mean_satt[2],3),'<br> var = ',round(var_satt[2],6)))
  )
  
  setting_labs <- c(
    Naive = "Naive",
    Satt_correct = "Satterthwaite~(hat(theta)[1]==-1~','~hat(theta)[2]==1)",
    Satt_M1 = "Satterthwaite~(hat(theta)[1]==-0.5~','~hat(theta)[2]==0.5)",
    Satt_M2 = "Satterthwaite~(hat(theta)[1]==-1.5~','~hat(theta)[2]==1.5)")
  
  lims <- range(quants$exact,quants$value)
  
  g <- ggplot(quants, aes(x = exact, y = value)) +
    geom_point(size=0.75) + geom_abline(intercept = 0,slope=1,lwd=0.5,col='gray40') +
    facet_wrap(~ approx, ncol = 2,labeller = as_labeller(setting_labs,label_parsed)) +
    xlab(lims) + ylab(lims) +
    ggtext::geom_richtext(
      data = params,
      size=3,
      aes(label = label),
      inherit.aes = FALSE,
      x = -Inf, y = Inf,
      hjust = -0.05, 
      vjust = 1.1
    )  + theme_bw() + theme(strip.text = element_text(size=12)) +
    labs(x = "Exact", y = 'Approximation') 
  
 print(g)


## wrong variance estimation ------
 set.seed(12526)
 N <- c(15,20)
 prob_D <- 0.2 # proportion in target domain

 mu <- c(rep(-1,N[1]),rep(1,N[2]))
 sig2_k <- 2
 
 dat <- data.frame(h = unlist(lapply(1:H,function(h){rep(h,N[h])})),
                   n = rpois(sum(N),10)) 
 dat$wt <- strata_totals[dat$h]/(N[dat$h]*dat$n)/1e3
 
 dat <- as.data.table(dat)
 
 S <- diag(1/dat$n)
 sig2_c <- sig2_k/mean(dat$n)
 
 # get quantities-----
 
 dat$D <- rbinom(nrow(dat),1,prob_D)
 
 if(sum(1-dat$D)>0)
   dat[dat$D==0,]$wt <- 0
 omega <- dat$wt*dat$n
 
 B_comps <- lapply(1:H,function(h){
   if(N[h]>0){
     return(N[h]/(N[h]-1)*(diag(1,N[h]) - 1/N[h]*matrix(1,N[h],N[h])))
   }
 })
 
 which.not.null <- c(1:H)[!sapply(B_comps, is.null)]
 B <- bdiag(B_comps[which.not.null])
 W <- diag(omega)%*%(diag(1,sum(N)) - (rep(1,sum(N))%*%t(omega))/sum(omega))
 M <- t(W)%*%B%*%W
 
 L <- chol(S)              # S = L' L
 Li <- backsolve(L, diag(nrow(S)))  # L^{-1}
 A <- t(L) %*% M %*% L
 out <- eigen(A, symmetric = TRUE)
 u <- Li %*% out$vectors
 
 keep_id <- which(out$values > 1e-10)
 q <- out$values[keep_id]
 r <- length(q)
 u <- u[,keep_id]
 delta <- (t(u)%*%mu)^2/sig2_k
 
 ## calculate dists under correct variances -----
 # exact
 dist_exact <- replicate(2e4,{
   sig2_k/(sum(omega)^2)*sum(q*sapply(1:r,function(i){rchisq(1,1,ncp =  delta[i])}))
 })
 
 q_exact <- quantile(dist_exact,probs=seq(0.01,0.99,0.01))
 
 #Satterthwhaite 
 Q1 <- sum(q*(1+delta))
 Q2 <- sum(q^2*(1+2*delta))
 
 scale_satt_correct <- sig2_k/sum(omega)^2*Q2/Q1
 df_satt_correct <- Q1^2/Q2
 mean_satt_correct <- scale_satt_correct*df_satt_correct
 var_satt_correct <- 2*scale_satt_correct^2*df_satt_correct
 
 dist_satt <- scale_satt_correct*rchisq(2e4,df_satt_correct)
 q_satt_correct <- quantile(dist_satt,seq(0.01,0.99,0.01))
 
 ## naive 
 
 df_naive_correct <- sum(dat$D) - length(unique(dat[dat$D==1,]$h))
 scale_naive_correct <- sig2_c/(sum(dat$D)*df_naive_correct)
 mean_naive_correct <- scale_naive_correct*df_naive_correct
 var_naive_correct <- 2*scale_naive_correct^2*df_naive_correct
 
 dist_naive <- scale_naive_correct*rchisq(2e4,df_naive_correct)
 q_naive_correct <- quantile(dist_naive,probs=seq(0.01,0.99,0.01))
 
 ## calculate dists under incorrect variances ------
 q_satt <- q_naive <- list() 
 scale_satt <- df_satt <- mean_satt <- var_satt <- scale_naive <- mean_naive <- var_naive <- NULL
 for(setting in c('V1','V2')){
   ## add estimation error ----
   if(setting == 'V1'){
     sig2_k <- 4
   }else if(setting == 'V2'){
     sig2_k <- 8
   }else{
     stop('error')
   }
   
   sig2_c <- sig2_k/10
   delta <- (t(u)%*%mu)^2/sig2_k
   
   ## Satterthwhaite ----
   
   Q1 <- sum(q*(1+delta))
   Q2 <- sum(q^2*(1+2*delta))
   
   dist_satt <- sig2_k/sum(omega)^2*Q2/Q1*rchisq(2e4,Q1^2/Q2)
   q_satt[[length(q_satt)+1]] <- quantile(dist_satt,seq(0.01,0.99,0.01))
   scale_satt[length(scale_satt)+1] <- sig2_k/sum(omega)^2*Q2/Q1
   df_satt[length(df_satt)+1] <- Q1^2/Q2
   mean_satt[length(mean_satt)+1] <- Q1^2/Q2*sig2_k/sum(omega)^2*Q2/Q1
   var_satt[length(var_satt)+1] <- 2*(sig2_k/sum(omega)^2*Q2/Q1)^2*Q1^2/Q2
   
   ## naive ----
   scale_naive[length(scale_naive) + 1] <- sig2_c/(sum(dat$D)*df_naive_correct)
   mean_naive[length(mean_naive) + 1] <- sig2_c/(sum(dat$D)*df_naive_correct)*df_naive_correct
   var_naive[length(var_naive) + 1] <- 2*(sig2_c/(sum(dat$D)*df_naive_correct))^2*df_naive_correct
   
   dist_naive <-  sig2_c/(sum(dat$D)*df_naive_correct)*rchisq(2e4,df_naive_correct)
   q_naive[[length(q_naive) + 1]] <- quantile(dist_naive,probs=seq(0.01,0.99,0.01))
   
   
 }
 
 
 ## plot -----
 quants <- data.frame(exact = q_exact,
                      Naive_correct = q_naive_correct,
                      Naive_V1 = q_naive[[1]],
                      Naive_V2 = q_naive[[2]],
                      Satt_V1 = q_satt[[1]],
                      Satt_V2 = q_satt[[2]],
                      Satt_correct = q_satt_correct) %>%
   pivot_longer(
     cols = c(Naive_correct,Naive_V1,Naive_V2, Satt_V1,Satt_V2,Satt_correct),
     names_to = "approx",
     values_to = "value"
   )
 
 params <- data.frame(
   approx = c("Naive_correct","Naive_V1","Naive_V2", "Satt_correct","Satt_V1","Satt_V2"),
   label  = c(paste0('scale = ',round(scale_naive_correct,5),'<br> df = ',df_naive_correct,'<br> mean = ',round(mean_naive_correct,3),'<br> var = ',round(var_naive_correct,6)),
              paste0('scale = ',round(scale_naive[1],5),'<br> df = ',df_naive_correct,'<br> mean = ',round(mean_naive[1],3),'<br> var = ',round(var_naive[1],6)),
              paste0('scale = ',round(scale_naive[2],5),'<br> df = ',df_naive_correct,'<br> mean = ',round(mean_naive[2],3),'<br> var = ',round(var_naive[2],6)),
              paste0('scale = ',round(scale_satt_correct,5),'<br> df = ',round(df_satt_correct,2),'<br> mean = ',round(mean_satt_correct,3),'<br> var = ',round(var_satt_correct,6)),
              paste0('scale = ',round(scale_satt[1],5),'<br> df = ',round(df_satt[1],2),'<br> mean = ',round(mean_satt[1],3),'<br> var = ',round(var_satt[1],6)),
              paste0('scale = ',round(scale_satt[2],5),'<br> df = ',round(df_satt[2],2),'<br> mean = ',round(mean_satt[2],3),'<br> var = ',round(var_satt[2],6)))
 )
 
 setting_labs <- c(Naive_correct = "Naive~(hat(sigma)^2==2)",
                   Naive_V1 = "Naive~(hat(sigma)^2==4)",
                   Naive_V2 = "Naive~(hat(sigma)^2==8)",
                   Satt_V1 = "Satterthwaite~(hat(sigma)^2==4)",
                   Satt_V2 = "Satterthwaite~(hat(sigma)^2==8)",
                   Satt_correct = "Satterthwaite~(hat(sigma)^2==2)")

 lims <- range(quants$exact,quants$value)
 g <- ggplot(quants, aes(x = exact, y = value)) +
   geom_point(size=0.75) + geom_abline(intercept = 0,slope=1,lwd=0.5,col='gray40') +
   facet_wrap(~ approx, ncol = 3,labeller = as_labeller(setting_labs,label_parsed)) +
   xlim(lims) + ylim(lims) +
   ggtext::geom_richtext(
     data = params,
     size=3,
     aes(label = label),
     inherit.aes = FALSE,
     x = -Inf, y = Inf,
     hjust = -0.05, 
     vjust = 1.1
   )  + theme_bw() + theme(strip.text = element_text(size=12)) +
   labs(x = "Exact", y = 'Approximation') 
 
 print(g)
 
 
 
 