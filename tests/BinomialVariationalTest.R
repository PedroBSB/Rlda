rm(list=ls())
setwd('C:\\Users\\Pedro Albuquerque\\Desktop\\Denis')
tmp=data.matrix(read.csv('fake data 7.csv',as.is=T))
res<- rlda.binomialVR(data=tmp, loc.id='loc.id', n_community=10, alpha0=0.01, alpha1=0.99, gamma=0.1, maxit=1000, thresh=0.0001)
