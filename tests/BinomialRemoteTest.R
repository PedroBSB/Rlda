rm(list=ls(all=T))
library('MCMCpack')
set.seed(1)

nbands=7
nloc=1000
ncommun=5
ndig.values=256 #number of digital values per pixel for remote sensing data

#theta: nloc * ncommun
if (ncommun==3) base=floor(nloc/(ncommun-1))
if (ncommun==5) base=floor(nloc/(ncommun-2))
if (ncommun==10) base=floor(nloc/(ncommun-4))

x=seq(from=-1,to=1,length.out=base)
y=sqrt(1-(x^2))*0.1
min1=0.0001
y[y<min1]=min1
# plot(x,y)

init=floor(nloc/ncommun)
seq1=c(seq(from=1,to=nloc,by=init),nloc)

theta=matrix(min1,nloc,ncommun)
for (i in 1:ncommun){
  seq2=seq1[i]:(seq1[i]+base-1)
  seq3=seq2[seq2<=nloc]
  theta[seq3,i]=y[1:length(seq3)]
}
theta=theta/matrix(apply(theta,1,sum),nloc,ncommun)
theta.true=theta

plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(1:nloc,theta[,i],col=i)

#omega matrix
# seq1=seq(from=0.95,to=0.1,length.out=ncommun)
omega.true=omega=matrix(runif(ncommun*nbands),ncommun,nbands)

#generate remote sensing data
remote=matrix(NA,nloc,nbands)
for (i in 1:nloc){
  for (j in 1:nbands){
    prob=theta[i,]%*%omega[,j]
    remote[i,j]=rbinom(1,size=ndig.values,prob=prob)
  }
}
image(remote)

pop<-as.data.frame(matrix(256,ncol=ncol(remote),nrow=nrow(remote)))
remote<-as.data.frame(remote)

data=remote
pop=pop
n_community=10
alpha0=1
alpha1=1
gamma=0.01
n_gibbs=1000
ll_prior = TRUE
display_progress = TRUE
Rcpp::sourceCpp("src/Utils.cpp")

#teste<- rlda.binomialRemote(data=remote, pop=pop, n_community=10, alpha0=1 , alpha1=1, gamma=0.01,
#                              n_gibbs=1000, ll_prior = TRUE, display_progress = TRUE)
