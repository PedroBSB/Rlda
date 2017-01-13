#######################################################################################################################
############################################   GENERATE FAKE DATA        ##############################################
#######################################################################################################################

rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(1)
nind=100
nloc=1000
ncomm=4
np=3
nspp=30

#design matrix
xmat=matrix(runif(nloc*np,min=-1,max=1),nloc,np)
xmat[,1]=1 #for intercept
colnames(xmat)=c('interc',paste('cov',1:(np-1),sep=''))

#regression parameters
betas=matrix(0,ncomm,np)
betas[1:ncomm,1]=c(-0.2,0.4,0.2,0)
betas[1:2,2]=c(-1,1)
betas[3:4,3]=c(-1,1)
betas.true=betas

#expand design matrix for each individual in each location
ind=rep(1:nloc,each=nind)
xmat1=xmat[ind,]

#generate y
y=expand.grid(ind=1:nind,loc=1:nloc)
y$y1=y$y2=y$y3=y$y4=NA

for (i in 1:ncomm){
  media=xmat1%*%betas[i,]
  nome=paste('y',i,sep='')
  y[,nome]=rnorm(nloc*nind,mean=media,sd=1)
}

#get z
nomes=paste('y',1:ncomm,sep='')
max1=apply(y[,nomes],1,max)
max2=matrix(max1,nloc*nind,ncomm)
z=apply(y[,nomes]==max2,1,which)
y$z=z

#generate phi
phi=matrix(NA,ncomm,nspp)
gamma=0.1
for (i in 1:ncomm){
  phi[i,]=rdirichlet(1,rep(gamma,nspp))
}
phi.true=phi

#generate w
y$w=NA
for (i in 1:ncomm){
  cond=y$z==i
  tmp=rmultinom(sum(cond),size=1,prob=phi[i,])
  ind=apply(tmp==1,2,which)
  y$w[cond]=ind
}
y.true=y

#generate aggregated dataset
agg=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  tmp=table(y[y$loc==i,'w'])
  aux=rep(0,nspp)
  aux[as.numeric(names(tmp))]=tmp
  agg[i,]=aux
}

#######################################################################################################################
############################################     GIBBS SAMPLING          ##############################################
#######################################################################################################################

yData<-as.data.frame(y[,c("y1","y2","y3","y4")])
xData<-as.data.frame(xmat1)
specVector<-y$w
gamma<-0.1
varBetas<-c(10,rep(1,np-1))
n_community<-ncomm
n_specie<-nspp

listYZW<-list(Y=as.matrix(yData),Z=rep(1,nrow(yData)),W=specVector)
xMat<-as.matrix(xData)

###Betas
betasMat<-generateBetas(listYZW, xMat, varBetas)

###Lista
phiMat<-phi
listYZW<-generateZ(listYZW, xMat, betasMat, phiMat)

###Phi
phiMat<-generatePhi(listYZW, xMat, gamma, n_specie)

###Gibbs
teste <- lda_multinomial_cov(yData, xData, specVector, varBetas, n_community, n_specie, gamma, 10000, TRUE, TRUE)
