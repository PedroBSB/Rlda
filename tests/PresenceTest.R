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

#Generate the Theta
alpha0<-0.01
alpha1<-0.01
Theta<-matrix(rbeta(species*community,alpha0,alpha1),nrow=species,ncol=community)

#Generate V
gamma<-0.01
vMat<-matrix(rbeta(locations*community,1,gamma),nrow=locations,ncol=community)
vMat[,community]<-1

#Generate Phi
Phi<-apply(vMat,1,function(x){
  prod<-1;
  Phi<-rep(NA,community)
  for(c in 1:community){
    vNumber <- x[c];
    if (c == 1) prod<-1;
    if (c >  1) prod<-prod*(1.0-x[c-1]);
    Phi[c]<-vNumber*prod;
  }
  Phi
})
Phi<-t(Phi)

#Generate Z
Z<-t(apply(Phi,1,function(x)rmultinom(1,1,x)))

#Generate Data
tmp<-as.data.frame(t(rep(0,species)))
DATA<-data.frame()
for(l in 1:locations){
  for(s in 1:species){
    tmp[s]<-rbinom(1,1,sum(Theta[s,]*Phi[l,]))
  }
  DATA<-rbind(DATA,tmp)
}

#Create rownames colnames
rownames(DATA)<-paste0("Location ",seq(1,nrow(DATA)))
colnames(DATA)<-paste0("Specie ",seq(1,ncol(DATA)))


###########################################################################################
#############################   GIBBS SAMPLING PROCEDURE    ###############################
###########################################################################################
alpha0<-0.01
alpha1<-0.01
gamma<-0.01
set.seed(9292)

#Estimate the parameters by Gibbs Sampling (Time difference of 18.73626 secs)
start.time <- Sys.time()

res<-lda_bernoulli(DATA, 4, alpha0, alpha1, gamma, 1000, TRUE, TRUE)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#Get the logLikelihood with Priors
logLikeli<-res$logLikelihood

#Plot the logLikelihood
plot(logLikeli,type="l")

tt<-res$Theta
