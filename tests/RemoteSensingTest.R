rm(list=ls(all=T))
library('MCMCpack')
librray('Rcpp')
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
#image(forest)

#### Gibbs functions
fix.MH=function(lo,hi,old1,new1,jump){
  jold=pnorm(hi,mean=old1,sd=jump)-pnorm(lo,mean=old1,sd=jump)
  jnew=pnorm(hi,mean=new1,sd=jump)-pnorm(lo,mean=new1,sd=jump)
  log(jold)-log(jnew) #add this to pnew
}
#----------------------------------------------------------------------------------------------
tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi

  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))

  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution

  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually

  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1

  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)

  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]
  accept <- length(keep)

  list(x = x0, accept = accept)
}
#-------------------------------
print.adapt = function(accept1z,jump1z){
  accept1=accept1z; jump1=jump1z;

  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }

  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
    jump1[[k]][cond] = jump1[[k]][cond]*2
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }

  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------------------------------
update.omega=function(param,jump){
  omega.orig=omega.old=param$omega
  tmp=tnorm(nbands*ncommun,lo=0,hi=1,mu=omega.old,sig=jump)
  novos=matrix(tmp,ncommun,nbands)

  tmp=fix.MH(lo=0,hi=1,omega.old,novos,jump)
  ajuste=matrix(tmp,ncommun,nbands)

  prior.old=matrix(dbeta(omega.old,a.omega,b.omega,log=T),ncommun,nbands)
  prior.new=matrix(dbeta(novos    ,a.omega,b.omega,log=T),ncommun,nbands)

  for (i in 1:ncommun){
    omega.new=omega.old
    omega.new[i,]=novos[i,]

    prob.old=param$theta%*%omega.old
    llk.old=colSums(dbinom(remote,size=ndig.values,prob=prob.old,log=T))

    prob.new=param$theta%*%omega.new
    llk.new=colSums(dbinom(remote,size=ndig.values,prob=prob.new,log=T))

    k=acceptMH(llk.old+prior.old[i,],
               llk.new+prior.new[i,]+ajuste[i,],
               omega.old[i,],omega.new[i,],F)
    omega.old[i,]=k$x
  }
  list(omega=omega.old,accept=omega.old!=omega.orig)
}
#------------------------------------------------
# make.from.SB.mat=function(vmat){
#   n=ncol(vmat)
#   mat=matrix(NA,nrow(vmat),n)
#   mat[,1]=vmat[,1]
#   mat[,2]=vmat[,2]*(1-vmat[,1])
#   for (i in 3:n){
#     mat[,i]=vmat[,i]*apply(1-vmat[,1:(i-1)],1,prod)
#   }
#   mat
# }

# vmat=matrix(rnorm(nloc*ncommun),nloc,ncommun)
# res=make.from.SB.mat(param$v)
# res1=convertSBtoNormal(vmat=param$v,
#                        ncol=ncol(param$v),nrow=nrow(param$v),
#                        prod=rep(1,nloc))
# range(res1-res)
#------------------------------------------------
update.theta=function(param,jump){
  v.orig=v.old=param$v
  tmp=tnorm(nloc*(ncommun-1),lo=0,hi=1,mu=v.old[,-ncommun],sig=jump[,-ncommun])
  novos=cbind(matrix(tmp,nloc,ncommun-1),1)
  ajuste=matrix(fix.MH(lo=0,hi=1,v.old,novos,jump),nloc,ncommun)

  prior.old=matrix(dbeta(v.old,1,gamma,log=T),nloc,ncommun)
  prior.new=matrix(dbeta(novos,1,gamma,log=T),nloc,ncommun)

  for (j in 1:(ncommun-1)){ #last column has to be 1
    v.new=v.old
    v.new[,j]=novos[,j]

    theta.old=convertSBtoNormal(vmat=v.old,ncol=ncommun,nrow=nloc,prod=rep(1,nloc))
    theta.new=convertSBtoNormal(vmat=v.new,ncol=ncommun,nrow=nloc,prod=rep(1,nloc))

    #contribution from reflectance data
    pold=theta.old%*%param$omega
    pnew=theta.new%*%param$omega
    p1.old=0#rowSums(dbinom(remote,size=ndig.values,pold,log=T))
    p1.new=0#rowSums(dbinom(remote,size=ndig.values,pnew,log=T))

    #contribution from forestry data
    p2.old=rowSums(forest*log(theta.old%*%param$phi))
    p2.new=rowSums(forest*log(theta.new%*%param$phi))

    k=acceptMH(p1.old+p2.old+prior.old[,j],
               p1.new+p2.new+prior.new[,j]+ajuste[,j],
               v.old[,j],v.new[,j],F)
    v.old[,j]=k$x
  }
  theta=convertSBtoNormal(vmat=v.old,ncol=ncommun,nrow=nloc,prod=rep(1,nloc))
  list(theta=theta,v=v.old,accept=v.old!=v.orig)
}
#------------------------------------------------
update.phi=function(param,jump){
  x.orig=x.old=param$x
  tmp=tnorm(ncommun*(nspp-1),lo=0,hi=1,mu=x.old[,-nspp],sig=jump[,-nspp])
  novos=cbind(matrix(tmp,ncommun,nspp-1),1)
  ajuste=matrix(fix.MH(lo=0,hi=1,x.old,novos,jump),ncommun,nspp)

  b.phi1=matrix(b.phi,ncommun,nspp,byrow=T)
  prior.old=dbeta(x.old,a.phi,b.phi1,log=T)
  prior.new=dbeta(novos,a.phi,b.phi1,log=T)

  unifs=matrix(runif(ncommun*nspp),ncommun,nspp)

  for (j in 1:(nspp-1)){ #last column has to be 1
    for (k in 1:ncommun){
      x.new=x.old
      x.new[k,j]=novos[k,j]

      phi.old=convertSBtoNormal(vmat=x.old,ncol=nspp,nrow=ncommun,prod=rep(1,ncommun))
      phi.new=convertSBtoNormal(vmat=x.new,ncol=nspp,nrow=ncommun,prod=rep(1,ncommun))

      #       ind=j:nspp
      #       pold=param$theta%*%phi.old[,ind]
      #       pnew=param$theta%*%phi.new[,ind]
      #       p1.old=sum(forest[,ind]*log(pold))
      #       p1.new=sum(forest[,ind]*log(pnew))

      pold=param$theta%*%phi.old
      pnew=param$theta%*%phi.new
      p1.old=sum(forest*log(pold))
      p1.new=sum(forest*log(pnew))

      k1=acceptMH(p1.old+prior.old[k,j],
                  p1.new+prior.new[k,j]+ajuste[k,j],
                  x.old[k,j],x.new[k,j],F)
      x.old[k,j]=k1$x
    }
  }
  #   z=updatePhi(x.old,x.old,novos,ajuste,prior.old,prior.new,unifs,
  #               nspp,ncommun,nloc,param$theta,forest)

  phi=convertSBtoNormal(vmat=x.old,ncol=nspp,nrow=ncommun,prod=rep(1,ncommun))
  list(phi=phi,x=x.old,accept=x.old!=x.orig)
}

###### Gibbs
#useful stuff
nind=rowSums(forest)
nloc=nrow(remote)
nbands=ncol(remote)
ndig.values=256 #number of digital values per pixel for remote sensing data
ncommun=10
nspp=ncol(forest)

#priors
gamma=0.01
a.omega=b.omega=1
psi=1
a.phi=psi
i=1:nspp
b.phi=psi*(nspp-i)

#initial values
phi=matrix(1/nspp,ncommun,nspp)
omega=matrix(runif(ncommun*nbands),ncommun,nbands)
theta=matrix(1/ncommun,nloc,ncommun)
v=theta
v[,ncommun]=1
x=phi
x[,nspp]=1

#stuff for gibbs sampling
ngibbs=1000
param=list(theta=theta,omega=omega,v=v,phi=phi,x=x)
vec.theta=matrix(NA,ngibbs,nloc*ncommun)
vec.omega=matrix(NA,ngibbs,ncommun*nbands)
vec.phi=matrix(NA,ngibbs,ncommun*nspp)

#stuff for MH algorithm
jump1=list(omega=matrix(1,ncommun,nbands),
           v=matrix(0.3,nloc,ncommun),
           x=matrix(0.05,ncommun,nspp))
accept1=list(omega=matrix(0,ncommun,nbands),
             v=matrix(0,nloc,ncommun),
             x=matrix(0,ncommun,nspp))
accept.output=50

Rcpp::sourceCpp('aux1.cpp')

i<-1

tmp=update.theta(param,jump1$v)
param$theta=tmp$theta #theta.true#
param$v=tmp$v
accept1$v=accept1$v+tmp$accept

tmp=update.omega(param,jump1$omega)
param$omega=tmp$omega #omega.true#
accept1$omega=accept1$omega+tmp$accept

tmp=update.phi(param,jump1$x)
param$phi=tmp$phi#phi.true#
param$x=tmp$x
accept1$x=accept1$x+tmp$accept

################################################################################################################################
################################################################################################################################
################################################################################################################################

thetaMat <- theta
Omega <- param$omega
jumpList <- list(omega=matrix(1,ncommun,nbands),
                v=matrix(0.3,nloc,ncommun),
                x=matrix(0.05,ncommun,nspp))
jumpAcceptance <- list(Omega=matrix(0,ncommun,nbands),
                 Theta=matrix(0,nloc,ncommun),
                 Phi=matrix(0,ncommun,nspp))
remoteMat <- as.matrix(remote)
maxBand <- 256
a <- 1.0
b<- 1.0
gamma <- 0.1

matX<-as.matrix(x)
forestMat<-as.matrix(forest)
bPhi <- b.phi
aPhi <- 1.0

vMatrix<-param$v

##Generate Omega
omega<-generateOmegaRemote(thetaMat, Omega, jumpList,jumpAcceptance, remoteMat, maxBand, a, b)

##Generate Phi
phi <- generatePhiRemote(thetaMat, matX, forestMat, jumpList,jumpAcceptance, bPhi, aPhi)


##Generate Theta
theta <- generateThetaRemote(vMatrix, omega,phi, forestMat,jumpAcceptance, jumpList, gamma)


################################################################################################################################
################################################################################################################################
################################################################################################################################

n_community<-ncommun
aOmega<- a
bOmega <-b
psi <- 1
accept_output <- 50
teste <- lda_remote(remoteMat, forestMat, jumpList, n_community, maxBand, gamma, aOmega, bOmega, psi, accept_output, 2000, TRUE)
save.image("C:\\Users\\p.albuquerque\\Desktop\\Resultados.Rdata")

################################################################################################################################
################################################################################################################################
################################################################################################################################

theta<-matrix(teste$Theta[2000,],nrow =n_community ,ncol=nrow(remoteMat))
phi<-matrix(teste$Phi[2000,],nrow = n_community,ncol=ncol(forestMat))
omega<-matrix(teste$Omega[2000,],nrow = n_community,ncol=ncol(remoteMat))

boxplot(theta)

rango=c(0,0.2)
plot(phi.true,phi[1:5,],xlim=rango,ylim=rango)
lines(rango,rango)

rango=c(0,1)
plot(omega.true,omega[1:5,],xlim=rango,ylim=rango)
lines(rango,rango)

plot(rep(1,5),1:5,col=1:5)
