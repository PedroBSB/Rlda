rm(list=ls())
setwd('C:\\Users\\Pedro Albuquerque\\Desktop\\Denis')
set.seed(2929)


source('VB functions.R')
tmp=data.matrix(read.csv('fake data 7.csv',as.is=T))

delta_elbo<-Inf
thresh<-0.0001
gamma<-0.1
a_phi<-0.01
b_phi<-0.99
n_obs<-5
maxit<-1000
ncommun<-10
nobs=5
ncommun=10


dat=aggregate.data(tmp)
nloc=nrow(dat)
nspp=ncol(dat)

m0=m1=array(abs(rnorm(nloc*nspp*ncommun)),dim=c(nloc,nspp,ncommun),
            dimnames=list(paste('loc',1:nloc,sep=''),
                          paste('spp',1:nspp,sep=''),
                          paste('comm',1:ncommun,sep='')))
tm1=tm0=m1=m0

a0=b0=a=b=matrix(1,nloc,ncommun)
c0=d0=c=d=matrix(1,ncommun,nspp)


##R: Get A and B
tmp=get.ab(ncommun,nloc,nspp,nobs,gamma,m1,m0,dat)
a=tmp$a
b=tmp$b
##RcppArmadillo: Get A and B
tmp0 = generateABmatrix(nloc,ncommun,nspp,nobs, gamma, tm1, tm0, dat)
a0<-tmp0$matA
b0<-tmp0$matB
sum(b0-b)

##R: Get M1 and M0
tmp=get.m1m0(nloc,ncommun,nspp,a,b,c,d)
m1=tmp$m1
m0=tmp$m0
##RcppArmadillo: Get M1 and M0
tmp0 = generateM1M0matrix(nloc,ncommun,nspp, a0, b0, c0, d0)
tm1<-tmp0$M1
tm0<-tmp0$M0



##R: Get C and D
tmp=get.cd(ncommun,nspp,a_phi,b_phi,nobs,dat,m1,m0)
c=tmp$c
d=tmp$d

##RcppArmadillo: Get C and D
tmp0 = generateCDmatrix(ncommun, nspp, nobs, a_phi, b_phi,  tm1,  tm0, dat)
c0<-tmp0$matC
d0<-tmp0$matD


##R: ELBO
get.elbo(a,b,c,d,m1,m0,dat,nobs,nloc,nspp,ncommun,gamma,a_phi,b_phi)

##RcppArmadillo: ELBO
generateELBO(a0, b0, c0, d0, tm1, tm0, dat, ncommun, nspp, nobs, nloc, gamma, a_phi, b_phi)





#-288674.3 -287268.2 -285876.1 -284498.4 -283135.2 -281786.8
teste<- lda_binomial_var(dat, ncommun, 100, n_obs, gamma, a_phi, b_phi, thresh, delta_elbo, m1, m0,  TRUE)
plot(teste$Elbo,type="l")
Theta<-matrix(teste$Theta[100,],nrow=nloc,ncol=ncommun)





