#Clean ethe Working Directory
rm(list=ls())

#Generate the fake data
set.seed(5587)
library(gtools)
library(Rlda)
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

#Generate V
gamma<-0.01
vMat<-matrix(rbeta(locations*community,1,gamma),nrow=locations,ncol=community)
vMat[,community]<-1

#Generate Theta
Theta<-apply(vMat,1,function(x){
  prod<-1;
  Theta<-rep(NA,community)
  for(c in 1:community){
    vNumber <- x[c];
    if (c == 1) prod<-1;
    if (c >  1) prod<-prod*(1.0-x[c-1]);
    Theta[c]<-vNumber*prod;
  }
  Theta
})
Theta<-t(Theta)

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


###########################################################################################
#############################   GIBBS SAMPLING PROCEDURE    ###############################
###########################################################################################
beta<-rep(1,species)
gamma<-0.01
set.seed(9292)


#Estimate the parameters by Gibbs Sampling (Time difference of 20.28974 secs)
start.time <- Sys.time()

res<-rlda.multinomial(as.data.frame(FIA), 4, beta, NA, 100,  FALSE, TRUE)
#res<-lda_multinomial_burn(FIA, 4, beta, NA, 1000, 50,  FALSE, TRUE)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#Get the logLikelihood with Priors
logLikeli<-res$logLikelihood

#Plot the logLikelihood
plot(logLikeli,type="l")

