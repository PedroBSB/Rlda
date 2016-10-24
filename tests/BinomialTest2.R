rm(list=ls())
library(Rlda)
#Read the SP500 data
data(sp500)
#Create size
spSize<-as.data.frame(matrix(100,
                             ncol=ncol(sp500),
                             nrow=nrow(sp500)))
sp5001<-na.omit(sp500)
#Set seed
set.seed(5874)
#Hyperparameters for each prior distribution
gamma <-0.01
alpha0<-0.01
alpha1<-0.01

nLocations<-nrow(sp500)
nBands<-ncol(sp500)
nCommunity<-5

#Intialize Theta
Theta<-matrix(rep(0.1,nLocations*nCommunity),nrow=nLocations,ncol=nCommunity)
vMat<-matrix(NA,nrow=nLocations,ncol=nCommunity)
Phi<-matrix(rep(0.1,nBands*nCommunity),nrow=nBands,ncol=nCommunity)


zMat  = generateZBinomial(as.matrix(sp500), as.matrix(spSize), Theta, Phi)
Theta = generateThetaBinomial(zMat,vMat, nLocations, nCommunity, gamma)
Phi   = generatePhiBinomial(zMat, nBands, nCommunity, alpha0, alpha1)


