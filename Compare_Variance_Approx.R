
load(file='/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/KEN_2022_cont_outcomes.rda')

cluster.sample <- dat %>% filter(!is.na(haz)) %>% dplyr::select(admin1,admin2,v025,v005,cluster) %>% group_by(admin1,admin2,v005,v025,cluster) %>% summarise(n=n())
cluster.sample$v005 <- cluster.sample$v005/1e6

cluster.sample <- cluster.sample %>% filter(admin1==1)
#cluster.sample <- cluster.sample[1:3,]

# unstratified admin1 -----------
## Let W = V/sigma^2
# number of clusters
N = nrow(cluster.sample)

omega <- cluster.sample$v005*cluster.sample$n ## Nx1
D = diag(omega^2)
M = diag(1,N) - (rep(1,N) %*% t(omega))/as.numeric(t(rep(1,N)) %*% t(t(omega)))
S = diag(1/cluster.sample$n) # structure of the covariance matrix
q <- eigen(S^(1/2)%*%t(M)%*%D%*%M%*%S^(1/2))$values

#pdf('/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/Variance_Approx_N=34.pdf')
par(mfrow=c(2,2))
# empirical
mu <- 1
sigma <- 2

W_emp <- replicate(5e4,{
  y <- as.vector(Rfast::rmvnorm(1,rep(mu,N),sigma^2*S))
  y_bar <- sum(omega*y)/sum(omega)
  N/(sum(omega)^2*(N-1))*sum(omega^2*(y-y_bar)^2)/sigma^2
})

# theoretical (exact)
W_exact <- replicate(5e4,{
  N/(sum(omega)^2*(N-1))*sum(q*rchisq(N,1))
})

plot(quantile(W_emp,seq(0.01,0.99,0.01)),quantile(W_exact,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Exact")
abline(0,1)

# theoretical (roughest approximation, df = N-1) -- very wrong!
W_rough_approx <- replicate(5e4,{
  sum(omega^2/cluster.sample$n)/(sum(omega)^2*(N-1))*rchisq(1,N-1)
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

omega <- matrix(cluster.sample$v005*cluster.sample$n,sum(N),1) ## Nx1

M <- diag(omega[,1]) %*% (diag(1,sum(N)) - ((omega %*% matrix(1,1,sum(N)))/(sum(omega))))

B_comps <- lapply(1:length(N),function(i){
  if(N[i]>1){
    (N[i]/(N[i]-1))*(diag(1,N[i]) - matrix(1,N[i],1)%*%matrix(1,1,N[i])/N[i])
  }else{
    matrix(1,1,1)
  }
  
})
B <- bdiag(B_comps)
C <- t(M)%*%B%*%M

## choose mean and variance
mu <- c(rep(-1,N[1]),
        rep(1,N[2]))
S <- diag(1,sum(N)) # difference in variance between strata can be accounted for here
sig2 <- 1


tmp <- eigen(S^(1/2)%*%C%*%S^(1/2))
q <- tmp$values
q <- q[abs(q)>1e-10]
U <- tmp$vectors
nu <- t(U) %*% S^(1/2)%*%C%*% mu
nu <- nu[abs(nu)>1e-10]
delta <- nu^2/q^2


dist_emp <- replicate(5e4,{
  y <- as.vector(Rfast::rmvnorm(1,mu,sqrt(sig2*S)))
  as.numeric(1/sum(omega)^2*t(y)%*%C%*%y/sig2)
})

#pdf('/Users/alanamcgovern/Desktop/Research/Project 2/FH Variance/Variance_Approx_Strat.pdf')

par(mfrow=c(2,2))

dist_exact <- replicate(5e4, 1/sum(omega)^2*sum(sapply(1:length(q),function(j){q[j]*rchisq(1,1,delta[j])})))

plot(quantile(dist_emp,seq(0.01,0.99,0.01)),quantile(dist_exact,seq(0.01,0.99,0.01)),
     xlab = 'Empirical',ylab='Theoretical', main = "Exact")
abline(0,1)

dist_rough <- 1/(sum(N)*(sum(N)-length(N)))*rchisq(5e4,sum(N)-length(N))

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
