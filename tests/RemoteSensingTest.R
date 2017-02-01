rm(list=ls(all=T))
library('MCMCpack')
set.seed(1)

nbands=7
nspp=40
nloc=1000
ncommun=5
nind=1000 #number of individuals per plot for forestry data
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

#phi matrix
phi=matrix(NA,ncommun,nspp)
for (i in 1:ncommun) phi[i,]=rdirichlet(1,rep(1,nspp))
phi.true=phi

#generate forestry data
forest=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  prob=theta[i,]%*%phi
  forest[i,]=rmultinom(1,size=nind,prob=prob)
}
image(forest)

#setwd('U:\\pedro\\remote sensing_FIA')
#write.csv(forest,'fake data forest.csv',row.names=F)
#write.csv(remote,'fake data remote.csv',row.names=F)
