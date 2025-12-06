
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
cluster.sample$urban <- ifelse(cluster.sample$v025==2,0,cluster.sample$v025)

N <- c(sum(cluster.sample$urban==1), #urban
        sum(cluster.sample$urban==0)) #rural 

omega <- cluster.sample$v005*cluster.sample$n ## Nx1

M <- diag(omega) %*% (diag(1,sum(N)) - (matrix(1,sum(N),ncol=1) %*% t(omega))/sum(omega))

# if both strata have samples
if(sum(N==0)==0){
  B_comps <- lapply(1:2,function(i){
    if(N[i]>1){
      (N[i]/(N[i]-1))*(diag(1,N[i]) - matrix(1,N[i],1)%*%matrix(1,1,N[i])/N[i])
    }else if(N[i]==1){
      matrix(1,1,1)
    }else{NULL} # shouldnt get here
  })
  
  B <- bdiag(B_comps)
}else{
  which.strat <- which(N!=0) # which strata has clusters?
  B <- (N[which.strat]/(N[which.strat]-1))*(diag(1,N[which.strat]) - matrix(1,N[which.strat],1)%*%matrix(1,1,N[which.strat])/N[which.strat])
}

C <- t(M)%*%B%*%M

## choose mean and variance
#mu <- rep(1,sum(N))
mu <- c(rep(-1,N[1]), rep(1,N[2]))
sig2 <- 2 
S <- diag(1/cluster.sample$n) # difference in variance between strata can be accounted for here
S_inv <- diag(1/diag(S))

out <- eigen((S)^(1/2)%*%C%*%(S)^(1/2))

q <- out$values[abs(out$values)>1e-10]
U <- out$vectors[,abs(out$values)>1e-10]


nu <- t(U) %*% sqrt(S_inv) %*%mu
delta <- nu^2/sig2


dist_emp <- replicate(5e4,{
  y <- as.vector(Rfast::rmvnorm(1,mu,sig2*S))
  as.numeric(1/sum(omega)^2*t(y)%*%C%*%y/sig2)
})

#pdf('/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/Variance_Approx_Strat.pdf')

par(mfrow=c(2,2))

dist_exact <- replicate(5e4, 1/sum(omega)^2*sum(sapply(1:length(q),function(j){q[j]*rchisq(1,1,delta[j])})))

plot(quantile(dist_emp,seq(0.01,0.99,0.01)),quantile(dist_exact,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Exact")
abline(0,1)

#dist_rough <- (sum(N)-length(N))/(sum(cluster.sample$n))*rchisq(5e4,sum(N)-length(N))
dist_rough <- 1/(sum(cluster.sample$n)*(sum(N)-length(N)))*rchisq(5e4,sum(N)-length(N))


plot(quantile(dist_emp,seq(0.01,0.99,0.01)),quantile(dist_rough,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Rough (#clusters - #strata)")
abline(0,1)


scale <- sum(q^2*(1 + 2*delta))/sum(q*(1+delta))
df <- (sum(q*(1+delta)))^2/sum(q^2*(1 + 2*delta))

dist_satt <- 1/sum(omega)^2*scale*rchisq(5e4,df)

plot(quantile(dist_emp,seq(0.01,0.99,0.01)),quantile(dist_satt,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Satterwhaite")
abline(0,1)


## not accounting for the non-centrality (difference in means between strata)
scale <- sum(q^2)/sum(q)
df <- (sum(q))^2/sum(q^2)

dist_satt_central <- 1/sum(omega)^2*scale*rchisq(5e4,df)

plot(quantile(dist_emp,seq(0.01,0.99,0.01)),quantile(dist_satt_central,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Satterwhaite (Central)")
abline(0,1)


#dev.off()

# unstratified admin2 -----------

tmp <- sampled_clusters %>% filter(admin1==40) %>% mutate(wt=wt/1e3)
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

