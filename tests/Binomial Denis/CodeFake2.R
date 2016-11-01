#Clean ethe Working Directory
rm(list=ls())
#Working directory
setwd("U:\\pedro\\binomial\\original r gibbs sampler")

#Generate the fake data
set.seed(5587)
library(gtools)

#Number of communities
community <- 3

#Number of Species
bands<-10

#Number of locations
locations<-1000

#Generate the Omega
alpha0<-0.1
alpha1<-0.1
Omega<-matrix(rbeta(bands*community,alpha0,alpha1),nrow=bands,ncol=community)
Theta=data.matrix(read.csv('theta.csv',as.is=T))

plot(NA,NA,xlim=c(1,locations),ylim=c(0,1))
for (i in 1:3) lines(1:locations,Theta[,i],col=i)

#Generate POP matrix
pop<-matrix(40,nrow=locations,ncol=bands)

#Generate Data
options(warn=2)
tmp<-rep(0,bands)
DATA<-matrix(NA,locations,bands)
for(l in 1:locations){
  for(b in 1:bands){
    prob=sum(Theta[l,]*Omega[b,])
    tmp[b]<-rbinom(1,pop[l,b],prob)
  }
  DATA[l,]<-tmp
}
image(DATA)

#Create rownames colnames
rownames(DATA)<-paste0("Location ",seq(1,nrow(DATA)))
colnames(DATA)<-paste0("Bands ",seq(1,ncol(DATA)))
POP<-as.data.frame(pop)
rownames(POP)<-paste0("Location ",seq(1,nrow(POP)))
colnames(POP)<-paste0("Bands ",seq(1,ncol(POP)))

write.csv(DATA,'fake data.csv')