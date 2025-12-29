
#load("~/Desktop/Research/Project 2/FHVariance_Smoothing/Kenya_Example/KEN_2022_cont_outcomes.rda")

# cluster.sample <- dat %>% filter(!is.na(haz)) %>% dplyr::select(admin1,admin2,v025,v005,cluster) %>% group_by(admin1,admin2,v005,v025,cluster) %>% summarise(n=n()) %>% arrange((v025))
# cluster.sample$v005 <- cluster.sample$v005/1e6


setwd("/Users/alanamcgovern/Desktop/Research/Project 2/FHVariance_Smoothing")

# load simulation objects -------------
setting <- 3
load(file=paste0('Simulations/Sim',setting,'/sim',setting,'_sampled_clusters.rda'))

cluster.sample <- sampled_clusters[[1]] %>% rename(v025 = urban, v005 = wt)

#cluster.sample <- cluster.sample[1:3,]

# unstratified admin1 -----------

tmp <- sampled_clusters[[k]] %>% filter(admin1==1) %>% mutate(wt=wt/1e3)
tmp$wt <- pmin(tmp$wt, quantile(tmp$wt,0.95))
N = nrow(tmp)

omega <- tmp$n*tmp$wt ## Nx1
D = diag(omega^2)
M = diag(1,N) - (rep(1,N) %*% t(omega))/sum(omega)
S = diag(1/tmp$n) # structure of the covariance matrix
q <- eigen(sqrt(S)%*%t(M)%*%D%*%M%*%sqrt(S))$values

#pdf('/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/Variance_Approx_N=34.pdf')
par(mfrow=c(2,2))
# empirical
mu <- -1
sigma <- 20

W_emp <- replicate(5e4,{
  y <- as.vector(Rfast::rmvnorm(1,rep(mu,N),sigma^2*S))
  y_bar <- sum(omega*y)/sum(omega)
  N/(sum(omega)^2*(N-1))*sum(omega^2*(y-y_bar)^2)
})

# theoretical (exact)
W_exact <- replicate(5e4,{
  sigma^2*N/(sum(omega)^2*(N-1))*sum(q*rchisq(N,1))
})

plot(quantile(W_emp,seq(0.01,0.99,0.01)),quantile(W_exact,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Exact")
abline(0,1)

# theoretical (roughest approximation, df = N-1) -- very wrong!
W_rough_approx <- replicate(5e4,{
  1/(N*(N-1))*rchisq(1,N-1)
})

plot(quantile(W_emp,seq(0.01,0.99,0.01)),quantile(W_rough_approx,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical',main = 'df=N-#clusters-1')
abline(0,1)


# theoretical (satterwhaite approximation) -- MUCH better, a bit skewed in the tails as number of clusters increases
W_satt_approx <- replicate(5e4, {
  N/(sum(omega)^2*(N-1))*(sum(q^2)/sum(q))*rchisq(1,sum(q)^2/sum(q^2))
})

plot(quantile(W_emp,seq(0.01,0.99,0.01)),quantile(W_satt_approx,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main='Satterhwaite approx')
abline(0,1)

# theoretical (pearson approximation) --- better than satterwhaite
gamma <- c(sum(q),
           sqrt(2*sum(q^2)),
           sqrt(8)*sum(q^3)/(sum(q^2)^(3/2)))

W_pears_approx <- replicate(5e4,{
  N/(sum(omega)^2*(N-1))*(gamma[1] - 2*gamma[2]/gamma[3] + rgamma(1,4/gamma[3]^2,scale=0.5*gamma[2]*gamma[3]))
})

plot(quantile(W_emp,seq(0.01,0.99,0.01)),quantile(W_pears_approx,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main='Pearson approx')
abline(0,1)
#dev.off()

# stratified admin1 ----------
area <- 6
#tmp <- sampled_clusters[[k]] %>% filter(admin1 == admin.key[admin.key$admin2 == area,]$admin1) %>% arrange(urban) #%>% mutate(wt=wt/1e2)

tmp <- sim_sample %>% filter(admin1 == admin.key[admin.key$admin2 == area,]$admin1) %>% group_by(admin2,urban,cluster) %>% reframe(yc = mean(value),n=n(),wt=unique(wt))
#  tmp$wt <- pmin(tmp$wt,quantile(tmp$wt,0.95))
tmp[tmp$admin2 != area,]$wt <- 0
N <- c(sum(1-tmp$urban),sum(tmp$urban))

#v_hat_scaled <- dir.dat$variance[area]*sum(tmp$wt^2*tmp$n)

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
S_inv <- diag(1/diag(S))

r <- rankMatrix(sqrt(S)%*%C%*%sqrt(S))
out <- eigen(sqrt(S)%*%C%*%sqrt(S))
q <- out$values[1:r]
U <- out$vectors[,1:r]

## choose mean and variance
 mu <- rep(1,sum(N))
 sig2 <- 2
# #mu <- c(rep(-1,N[1]), rep(1,N[2]))
# V <- 0.05
# sig2 <- V*sum(omega)^2/sum(tmp$wt^2*tmp$n)

# mu <- rep(admin2_strata_mean[area],sum(N))
# #mu <- admin2_strata_mean[tmp$admin2]
# V <- sum(tmp$wt^2*tmp$n)/sum(omega)^2*sig2[area]
 
 V <- sum(tmp$wt^2*tmp$n)/sum(omega)^2*sig2

nu <- sqrt(sum(tmp$wt^2*tmp$n))/sum(omega)*t(U)%*%sqrt(S_inv)%*%mu
delta <- nu^2/V


#(1/sum(omega)^2*t(tmp$yc)%*%C%*%tmp$yc)

# tmp <- data.table::as.data.table(tmp)
# shell_dat <- tmp[rep(seq_len(sum(N)),times = tmp$n)]
# 
dist_emp <- replicate(5e3,{
  shell_dat$y <- rnorm(nrow(shell_dat),mu,sqrt(sig2))
  y_c <- shell_dat[,mean(y),by=cluster]$V1
  as.numeric(1/sum(omega)^2*t(y_c)%*%C%*%y_c)
#  as.numeric(1/sum(omega)^2*t(y)%*%C%*%y/sig2)
})

dist_emp <- replicate(5e3,{
  y_c <- rnorm(sum(N),mu,sqrt(sig2))
  y_c <- shell_dat[,mean(y),by=cluster]$V1
  as.numeric(1/sum(omega)^2*t(y_c)%*%C%*%y_c)
  #  as.numeric(1/sum(omega)^2*t(y)%*%C%*%y/sig2)
})

#pdf('/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/Variance_Approx_Strat.pdf')

#par(mfrow=c(2,2))

#dist_exact <- replicate(5e4, V/sum(tmp$wt^2*tmp$n)*sum(sapply(1:length(q),function(j){q[j]*rchisq(1,1,delta[j])})))
dist_exact <- replicate(5e4, sig2/sum(omega)^2*sum(sapply(1:length(q),function(j){q[j]*rchisq(1,1,delta[j])})))

hist(dist_exact)

plot(quantile(dist_emp,seq(0.01,0.99,0.01)),quantile(dist_exact,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Exact")
abline(0,1)

dist_rough <- (sum(N)-length(N))/(sum(cluster.sample$n))*rchisq(5e4,sum(N)-length(N))
#dist_rough <- V/(sum(tmp$wt>0)-2)*rchisq(5e4,sum(tmp$wt>0)-2)


plot(quantile(dist_exact,seq(0.01,0.99,0.01)),quantile(dist_rough,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Rough (#clusters - #strata)")
abline(0,1)


scale <- sum(q^2*(1 + 2*delta))/sum(q*(1+delta))
df <- (sum(q*(1+delta)))^2/sum(q^2*(1 + 2*delta))

dist_satt <- V/sum(tmp$wt^2*tmp$n)*scale*rchisq(5e4,df)

plot(quantile(dist_exact,seq(0.01,0.99,0.01)),quantile(dist_satt,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Satterwhaite")
abline(0,1)


## not accounting for the non-centrality (difference in means between strata)
scale <- sum(q^2)/sum(q)
df <- (sum(q))^2/sum(q^2)

dist_satt_central <- V/sum(tmp$wt^2*tmp$n)*scale*rchisq(5e4,df)

plot(quantile(dist_exact,seq(0.01,0.99,0.01)),quantile(dist_satt_central,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Satterwhaite (Central)")
abline(0,1)


#dev.off()

# unstratified admin2 -----------

tmp <- sampled_clusters[[1]] %>% filter(admin1==40) %>% mutate(wt=wt/1e3)
N = nrow(tmp)

tmp[tmp$admin2 != 261,]$wt <- 0

omega <- tmp$n*tmp$wt ## Nx1
S = diag(1/tmp$n) # structure of the covariance matrix

#pdf('/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/Variance_Approx_N=34.pdf')
par(mfrow=c(2,2))
# empirical
mu <- -1
sigma <- 2

W_emp <- replicate(5e4,{
  y <- as.vector(Rfast::rmvnorm(1,rep(mu,N),sigma^2*S))
  y_bar <- sum(omega*y)/sum(omega)
  N/(sum(omega)^2*(N-1))*sum(omega^2*(y-y_bar)^2)/sigma^2
})

tmp$wt <- pmin(tmp$wt, quantile(tmp$wt,0.95))
omega <- tmp$n*tmp$wt ## Nx1
D = diag(omega^2)
M = diag(1,N) - (rep(1,N) %*% t(omega))/sum(omega)
q <- eigen(sqrt(S)%*%t(M)%*%D%*%M%*%sqrt(S))$values

# theoretical (exact)
W_exact <- replicate(5e4,{
  N/(sum(omega)^2*(N-1))*sum(q*rchisq(N,1))
})

plot(quantile(W_emp,seq(0.01,0.99,0.01)),quantile(W_exact,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Exact")
abline(0,1)

# theoretical (roughest approximation, df = N-1) -- very wrong!
W_rough_approx <- replicate(5e4,{
  1/(N*(N-1))*rchisq(1,N-1)
})

plot(quantile(W_emp,seq(0.01,0.99,0.01)),quantile(W_rough_approx,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical',main = 'df=N-#clusters-1')
abline(0,1)


# theoretical (satterwhaite approximation) -- MUCH better, a bit skewed in the tails as number of clusters increases
W_satt_approx <- replicate(5e4, {
  N/(sum(omega)^2*(N-1))*(sum(q^2)/sum(q))*rchisq(1,sum(q)^2/sum(q^2))
})

plot(quantile(W_emp,seq(0.01,0.99,0.01)),quantile(W_satt_approx,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main='Satterhwaite approx')
abline(0,1)



# stratified admin2 ------------
area <- 2

# tmp <- dir.dat %>% filter(admin1==admin.key[admin.key$admin2==area,]$admin1) %>% group_by(admin1,admin2,v005,v025,cluster) %>% 
#   summarise(n=n()) %>% arrange(desc(v025))
# tmp$wt <- tmp$v005/1e6
# tmp$urban <- ifelse(tmp$v025==2,0,tmp$v025)

tmp <- sampled_clusters[[k]] %>% filter(admin1 == admin.key[admin.key$admin2 == area,]$admin1) %>% arrange(urban)

tmp[tmp$admin2 != area,]$wt <- 0
# tmp$wt <- pmin(tmp$wt,quantile(tmp[tmp$wt>0,]$wt,0.95))
N <- c(sum(1-tmp$urban),sum(tmp$urban))

tmp_D <- dir.dat %>% filter(admin2==area) %>% group_by(cluster,v025) %>% summarise(n=n())
N_D <- c(sum(tmp_D$v025==2),sum(tmp_D$v025==1))

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

mu <- c(0,1)
mu_vec <- c(rep(mu[1],N[1]),rep(mu[2],N[2]))
sig2 <- 2
kappa <- 3
V <- as.numeric(sig2*t(omega)%*%S%*%omega/(t(omega)%*%omega))
S <- diag((1+kappa*tmp$urban)/tmp$n)
S_inv <- diag(1/diag(S))


## generate empirical distribution
dist_emp <- replicate(1e4,{
  y_c <- as.vector(Rfast::rmvnorm(1,mu_vec,sig2*S))
  as.numeric(1/t(omega)%*%omega*t(y_c)%*%C%*%(y_c))
})


r <- rankMatrix(sqrt(S)%*%C%*%sqrt(S))
q <-  eigen(sqrt(S)%*%C%*%sqrt(S))$values[1:r]
U <-  eigen(sqrt(S)%*%C%*%sqrt(S))$vector[,1:r]
nu <- as.numeric(sqrt(t(omega)%*%S%*%omega)/sqrt(t(omega)%*%omega))*(t(U)%*%sqrt(S_inv)%*%tmp$urban)*(mu[2]-mu[1])
delta <- nu^2/V

## exact distribution
dist_exact <- replicate(5e4,{
  V/(t(omega)%*%S%*%omega)*sum(q*sapply(1:r,function(i){rchisq(1,1,ncp =  delta[i])}))
})

# plot(quantile(dist_emp,seq(0.01,0.99,0.01)),quantile(dist_exact,seq(0.01,0.99,0.01)),
#      xlab = 'Empirical',ylab='Theoretical', main = "Exact")
# abline(0,1)
par(mfrow=c(1,2))
## naive approx
dist_naive <- V/(sum(N_D) - sum(N_D>0))*rchisq(5e4,sum(N_D) - sum(N_D>0))
plot(quantile(dist_exact,seq(0.01,0.99,0.01)),quantile(dist_naive,seq(0.01,0.99,0.01)),
     xlab = 'Exact',ylab='Approximation', main = "Naive (scale = 0.2, df=5)")
abline(0,1)

## saterwhaite approx
Q1 <- sum(q*(1+delta))
Q2 <- sum(q^2*(1+2*delta))

dist_satt <- V/as.numeric(t(omega)%*%S%*%omega)*Q2/Q1*rchisq(5e4,Q1^2/Q2)
plot(quantile(dist_exact,seq(0.01,0.99,0.01)),quantile(dist_satt,seq(0.01,0.99,0.01)),
     xlab = 'Exact',ylab='Approximation', main = "Satterwhaite (scale=0.29, df=5)")
abline(0,1)

## saterwhaite approx without accounting for anything strata specific ------------
S <- diag(1/tmp$n)
S_inv <- diag(1/diag(S))
r <- rankMatrix(sqrt(S)%*%C%*%sqrt(S))
q <-  eigen(sqrt(S)%*%C%*%sqrt(S))$values[1:r]

Q1 <- sum(q)
Q2 <- sum(q^2)

dist_satt_nokappa <- V/as.numeric(t(omega)%*%S%*%omega)*Q2/Q1*rchisq(5e4,Q1^2/Q2)
plot(quantile(dist_exact,seq(0.01,0.99,0.01)),quantile(dist_satt_nokappa,seq(0.01,0.99,0.01)),
     xlab = 'Exact',ylab='Approximation', main = "Satterwhaite")
abline(0,1)



plot(quantile(rchisq(5e4,Q1^2/Q2),seq(0.01,0.99,0.01)),
     quantile(rchisq(5e4,sum(N_D) - sum(N_D>0)),seq(0.01,0.99,0.01)))
     

