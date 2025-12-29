# A = all assumptions met
# V1 = unplanned domain
# V2 = different strata sample size
# V3 = not identically distributed

# all assumptions of Naive approximation are met -------

## settings 
H <- 2
strata_totals <- c(500,500)

## correct specification and estimation -----
N <- c(10,10)
prob_D <- 1 # proportion in target domain

mu <- c(rep(-1,N[1]),rep(1,N[2]))
sig2_k <- 2 #individual level variance
sig2_c <- sig2_k/10

dat <- data.frame(h = unlist(lapply(1:H,function(h){rep(h,N[h])})),
                  n = rpois(sum(N),10)) 
dat$wt <- strata_totals/(N[dat$h]*dat$n)/10

dat <- as.data.table(dat)

S <- diag(1/10,sum(N))

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

dist_exact <- replicate(5e4,{
  sig2_k/(sum(omega)^2)*sum(q*sapply(1:r,function(i){rchisq(1,1,ncp =  delta[i])}))
})

q_exact <- quantile(dist_exact,probs=seq(0.01,0.99,0.01))

## compare to naive 

df_naive <- sum(dat$D) - length(unique(dat[dat$D==1,]$h))
dist_naive <- sig2_c/(sum(dat$D)*df_naive)*rchisq(5e4,df_naive)

q_naive <- quantile(dist_naive,probs=seq(0.01,0.99,0.01))

## compare to satterwhaite 

Q1 <- sum(q*(1+delta))
Q2 <- sum(q^2*(1+2*delta))

dist_satt <- sig2_k/(sum(omega)^2)*Q2/Q1*rchisq(5e4,Q1^2/Q2)
q_satt <- quantile(dist_satt,seq(0.01,0.99,0.01))

## plot 
quants <- data.frame(exact = q_exact,
                     Naive = q_naive,
                     Satterwhaite = q_satt) %>%
  pivot_longer(
    cols = c(Naive, Satterwhaite),
    names_to = "approx",
    values_to = "value"
  )

params <- data.frame(
  approx = c("Naive", "Satterwhaite"),
  label  = c(paste0('scale = ',round(sig2_c/(sum(dat$D)*df_naive),4),'\n df = ',df_naive), 
             paste0('scale = ',round(sig2_k/(sum(omega)^2)*Q2/Q1,4),'\n df = ',round(Q1^2/Q2,2)))
)

g <- ggplot(quants, aes(x = exact, y = value)) +
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

for(setting in c('V1','V2','V3','V4')){

if(setting =='V3'){
  N <- c(10,10)
  prob_D <- 1 # proportion in target domain
}else if(setting == 'V1'){
  N <- rep(50,H) # in planned
  prob_D <- 0.2 # proportion in target domain
}else if(setting == 'V2'){
  N <- c(5,15)
  prob_D <- 1 # proportion in target domain
}else if(setting =='V4'){
 # N <- c(5,15)
#  prob_D <- 1 # proportion in target domain
   N <- c(5*5,5*15)
  prob_D <- 0.2 # proportion in target domain
}else{ 
  stop('error')
}
 
mu <- c(rep(-1,N[1]),rep(1,N[2]))
sig2_k <- 2
sig2_c <- sig2_k/10

dat <- data.frame(h = unlist(lapply(1:H,function(h){rep(h,N[h])})),
                  n = rpois(sum(N),10)) 
dat$wt <- strata_totals/(N[dat$h]*dat$n)

dat <- as.data.table(dat)

# setting A, V1, V2
if(setting %in% c('V1','V2')){
  S <- diag(1/10,sum(N))
}else if(setting %in% c('V3','V4')){
  S <- diag(1/dat$n)
}else{
  stop('error')
}

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

# generate empirical -------

dist_emp <- replicate(5e4,{
  y <- Rfast::rmvnorm(1,mu,sig2_k*S)
  return(as.numeric(y%*%M%*%t(y)/sum(omega)^2))
})

q_emp <- quantile(dist_emp,probs=seq(0.01,0.99,0.01))

## calculate exact -----
dist_exact <- replicate(5e4,{
  sig2_k/(sum(omega)^2)*sum(q*sapply(1:r,function(i){rchisq(1,1,ncp =  delta[i])}))
})

q_exact <- quantile(dist_exact,probs=seq(0.01,0.99,0.01))

## compare to satterwhaite -----

Q1 <- sum(q*(1+delta))
Q2 <- sum(q^2*(1+2*delta))

dist_satt <- sig2_k/sum(omega)^2*Q2/Q1*rchisq(5e4,Q1^2/Q2)
q_satt <- quantile(dist_satt,seq(0.01,0.99,0.01))

## compare to naive -----

df_naive <- sum(dat$D) - length(unique(dat[dat$D==1,]$h))
dist_naive <- sig2_c/(sum(dat$D)*df_naive)*rchisq(5e4,df_naive)

q_naive <- quantile(dist_naive,probs=seq(0.01,0.99,0.01))

## plot -----
quants <- data.frame(exact = q_exact,
                     Naive = q_naive,
                     Satterwhaite = q_satt) %>%
  pivot_longer(
    cols = c(Naive, Satterwhaite),
    names_to = "approx",
    values_to = "value"
  )

params <- data.frame(
  approx = c("Naive", "Satterwhaite"),
  label  = c(paste0('scale = ',round(sig2_c/(sum(dat$D)*df_naive),4),'\n df = ',df_naive), 
             paste0('scale = ',round(sig2_k/sum(omega)^2*Q2/Q1,4),'\n df = ',round(Q1^2/Q2,2)))
)

g_tmp <- ggplot(quants, aes(x = exact, y = value)) +
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

if(setting=='V1'){
  g1 <- g_tmp + ggtitle('A. Unplanned domain')
}else if(setting=='V2'){
  g2 <- g_tmp + ggtitle('B. Unequal sample size and weights')
}else if(setting=='V3'){
  g3 <- g_tmp + ggtitle('C. Unequal variance')
}else if(setting=='V4'){
  g4 <- g_tmp + ggtitle('A + B + C')
}else{
  stop('error')
}

}

ggarrange(plotlist = list(g1,g2,g3,g4),nrow=4)


