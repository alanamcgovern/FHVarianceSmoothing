
## TWO STAGE PLANNED (unstratified) ------------
which.area <- 1
dat.tmp <- dir.dat %>% filter(admin1==which.area)
dat.tmp$v005 <- dat.tmp$v005/1e6

my.svydesign <- survey::svydesign(ids = ~cluster, 
                                  strata = ~admin1,
                                  weights = ~v005, data = dat.tmp)
lm.ob <- survey::svymean(value ~ 1, design = my.svydesign)
lm.ob

vcov(lm.ob)

# get mean
y_bar <- sum(dat.tmp$value*dat.tmp$v005)/sum(dat.tmp$v005)

# get cluster totals of linearized variable
cluster_dt <- dat.tmp %>% group_by(cluster) %>% reframe(wt=unique(v005),yc=mean(value),n=n()) %>% mutate(u = wt*n*(yc-y_bar))

cluster_dt %>% summarise(n()/(n()-1)*(1/sum(n*wt)^2)*sum((u-mean(u))^2))

## TWO STAGE PLANNED (stratified) ------------
which.area <- 1
dat.tmp <- dir.dat %>% filter(admin1==which.area)

my.svydesign <- survey::svydesign(ids = ~cluster, 
                                  strata = ~v023,
                                  weights = ~v005, data = dat.tmp)
lm.ob <- survey::svymean(value ~ 1, design = my.svydesign)
lm.ob

vcov(lm.ob)

# get mean
y_bar <- sum(dat.tmp$value*dat.tmp$v005)/sum(dat.tmp$v005)

wt_sum <- sum(dat.tmp$v005)

cluster_dt <- dat.tmp %>% group_by(cluster,v025) %>% reframe(wt=unique(v005),yc=mean(value),n=n()) %>% mutate(u = wt*n*(yc-y_bar))
cluster_dt %>% group_by(v025) %>% summarise(t = n()/(n()-1)*(1/wt_sum^2)*sum((u-mean(u))^2)) %>% summarise(sum(t))

## TWO STAGE UNPLANNED (unstratified) -------
dat.tmp <- dir.dat %>% filter(admin1==1)

my.svydesign <- survey::svydesign(ids = ~cluster, 
                                  strata = ~admin1,
                                  weights = ~v005, data = dat.tmp)

tmp <- subset(my.svydesign, (admin2 == as.character(2)))

lm.ob <- survey::svymean(value ~ 1, design = tmp)
lm.ob
vcov(lm.ob)

dat.tmp[dat.tmp$admin2 != 2,]$v005 <- 0

# get mean
y_bar <- sum(dat.tmp$value*dat.tmp$v005)/sum(dat.tmp$v005)


cluster_dt <- dat.tmp %>% group_by(cluster) %>% reframe(wt=unique(v005),yc=mean(value),n=n())
N <- nrow(cluster_dt)

N <- sum(cluster_dt$wt>0)

cluster_dt %>% summarise(N/(N-1)*(1/sum(n*wt)^2)*sum(wt^2*n^2*(yc-y_bar)^2))

## TWO STAGE UNPLANNED (stratified) ----------
dat.tmp <- dir.dat %>% filter(admin1==1)

my.svydesign <- survey::svydesign(ids = ~cluster, 
                                  strata = ~v023,
                                  weights = ~v005, data = dat.tmp)

tmp <- subset(my.svydesign, (admin2 == as.character(3)))

lm.ob <- survey::svymean(value ~ 1, design = tmp)
lm.ob
vcov(lm.ob)

dat.tmp[dat.tmp$admin2 != 3,]$v005 <- 0

# get mean
y_bar <- sum(dat.tmp$value*dat.tmp$v005)/sum(dat.tmp$v005)

wt_sum <- sum(dat.tmp$v005)

cluster_dt <- dat.tmp %>% group_by(cluster,v025) %>% reframe(wt=unique(v005),yc=mean(value),n=n()) %>% mutate(u = wt*n*(yc-y_bar))
cluster_dt %>% group_by(v025) %>% summarise(t = n()/(n()-1)*(1/wt_sum^2)*sum((u-mean(u))^2)) %>% summarise(sum(t))

S <- diag(1/cluster_dt$n)
omega <- cluster_dt$wt*cluster_dt$n

N <- c(sum(cluster_dt$v025==2),sum(cluster_dt$v025==1))
B_comps <- lapply(1:2,function(h){
  if(N[h]>0){
    return(N[h]/(N[h]-1)*(diag(1,N[h]) - 1/N[h]*matrix(1,N[h],N[h])))
  }
})

which.not.null <- c(1:2)[!sapply(B_comps, is.null)]
B <- bdiag(B_comps[which.not.null])
W <- diag(omega)%*%(diag(1,sum(N)) - (rep(1,sum(N))%*%t(omega))/sum(omega))
M <- t(W)%*%B%*%W

as.numeric(t(cluster_dt$yc)%*%M%*%cluster_dt$yc/sum(omega)^2)

