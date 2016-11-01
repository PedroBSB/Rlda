#Clean ethe Working Directory
rm(list=ls())
set.seed(17)
library(gtools)

#Working directory
setwd("U:\\pedro\\binomial\\original r gibbs sampler")
source('gibbs functions.R')
dat=data.matrix(read.csv('fake data.csv',as.is=T)[,-1])

#Number of communities
community <- 10

#Number of bands
bands<-ncol(dat)

#Number of locations
locations<-nrow(dat)

#priors
a1=1; a2=1
gamma=1

#initial values
Theta=matrix(1/community,locations,community)
Omega=matrix(runif(community*bands),bands,community)

#assign everybody to community 1
max.dat=max(dat)
zmat=expand.grid(comm=1,bands=1:bands,reflect=0:1,loc=1:locations,n=NA)
for (i in 1:bands){
  cond=zmat$bands==i & zmat$reflect==1
  zmat[cond,'n']=dat[,i]
  cond=zmat$bands==i & zmat$reflect==0
  zmat[cond,'n']=max.dat-dat[,i]
}

param=list(Theta=Theta,Omega=Omega,zmat=zmat)
ngibbs=1000
vec.omega=matrix(NA,ngibbs,community*bands)
vec.theta=matrix(NA,ngibbs,community*locations)
vec.loglikel=matrix(NA,ngibbs,1)

for (i in 1:ngibbs){
  print(i)
  print(aggregate(n~comm,data=param$zmat,sum)$n)
  tmp=update.Theta(param)
  param$Theta=tmp$Theta
  param$Vjl=tmp$Vjl
  
  param$Omega=update.Omega(param)
  param$zmat=update.zmat(param)
  
  vec.omega[i,]=param$Omega
  vec.theta[i,]=param$Theta
  vec.loglikel[i]=calc.loglikel(param)
}

boxplot(param$Theta)

plot(NA,xlim=c(0,locations),ylim=c(0,1))
for (i in c(1,2,4)) lines(param$Theta[,i],col=i)
