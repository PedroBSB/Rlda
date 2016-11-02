#Clean ethe Working Directory
rm(list=ls())
set.seed(17)
library(gtools)
library(Rlda)

#Working directory
#source("C:\\Users\\p.albuquerque\\Dropbox\\Rlda\\tests\\Binomial Denis\\gibbs functions.R")
dat=data.matrix(read.csv('C:\\Users\\p.albuquerque\\Dropbox\\Rlda\\tests\\Binomial Denis\\fake data.csv',as.is=T)[,-1])

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

pop<-matrix(40,nrow=nrow(dat),ncol=ncol(dat))
res<- lda_binomial(dat,pop,community,a1,a2,gamma,100,TRUE,TRUE)


test<-generateZBinomial(dat,pop, Theta, Omega)

sum<-0
for(i in 1:length(test)){
  z<-test[[i]]
  which(z$Community==0)
}

#Rcpp::sourceCpp("src/Binomial.cpp")
#res<- lda_binomial(dat,pop,community,a1,a2,gamma,100,TRUE,FALSE)


#param=list(Theta=Theta,Omega=Omega,zmat=zmat)
#ngibbs=1000
#vec.omega=matrix(NA,ngibbs,community*bands)
#vec.theta=matrix(NA,ngibbs,community*locations)
#vec.loglikel=matrix(NA,ngibbs,1)

#for (i in 1:ngibbs){
#  print(i)
#  print(aggregate(n~comm,data=param$zmat,sum)$n)
#  tmp=update.Theta(param)
#  param$Theta=tmp$Theta
#  param$Vjl=tmp$Vjl

#  param$Omega=update.Omega(param)
#  param$zmat=update.zmat(param)

#  vec.omega[i,]=param$Omega
#  vec.theta[i,]=param$Theta
#  vec.loglikel[i]=calc.loglikel(param)
#}
Theta<-matrix(res$Theta[100,],ncol=community,nrow=1000)
#boxplot(param$Theta)
boxplot(Theta)

plot(NA,xlim=c(0,locations),ylim=c(0,1))
for (i in c(1,2,4)) lines(param$Theta[,i],col=i)
