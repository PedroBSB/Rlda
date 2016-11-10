
#Clean ethe Working Directory
rm(list=ls())

#Generate the fake data
set.seed(5587)
library(gtools)

#Number of communities
community <- 4

#Number of Species
species<-100

#Number of locations
locations<-1000

#Number of individuals in each location
size<-rpois(locations,40)

#Define the Hyperparameters
beta<-rep(1,species)

#Generate Phi
Phi<-rdirichlet(community,beta)

#Number of covariates
nCovariates <- 5

#Design matrix
designMat<-MASS::mvrnorm(locations,rep(0,nCovariates),diag(1,nCovariates,nCovariates))

Beta<-list()
Theta<-matrix(NA,locations,community)
#Generate Beta
for(l in 1:locations){
  betaMat<-MASS::mvrnorm(community,seq(1,nCovariates),diag(1,nCovariates,nCovariates))
  Beta[[l]]<-betaMat
  for(c in 1:community){
    Theta[l,c]<-exp(sum(designMat[l,]*betaMat[c,]))
  }
}

#Normalize Theta
Theta<-t(apply(Theta,1,function(x)x/sum(x)))

#Generate Z
tmp<-cbind(size,Theta)
Z<-t(apply(tmp,1,function(x)rmultinom(1, x[1], x[2:(community+1)])))

#Generate S
FIA<-matrix(NA,locations,species)
for(l in 1:locations){
  tmp<-rep(0,species)
  for(c in 1:community){
    tmp<-tmp+t(as.data.frame(rmultinom(1, Z[l,c], Phi[c,])))
  }
  FIA[l,]<-tmp
}

#Create rownames colnames
rownames(FIA)<-paste0("Location ",seq(1,nrow(FIA)))
colnames(FIA)<-paste0("Specie ",seq(1,ncol(FIA)))


res<-lda_covariate(FIA,as.data.frame(designMat), community, beta, 100,TRUE,TRUE )

gg<-res[[100]]
