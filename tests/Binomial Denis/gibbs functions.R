update.Theta=function(param){
  Theta=Vjl=matrix(NA,locations,community)
  zmat=param$zmat
  
  #calculate the number of individuals in each community and locations
  res=aggregate(n~comm+loc,data=zmat,sum)
  #make sure that empty communities have zeroes for 'n'
  aux=expand.grid(loc=1:locations,comm=1:community)
  res1=merge(res,aux,all=T)
  cond=is.na(res1$n)
  res1[cond,'n']=0
  
  for (i in 1:(community-1)){
    cond=res1$comm==i
    njl=res1[cond,'n']
    
    #calculate sum of individuals from communities greater than i
    cond=res1$comm>i
    tmp=aggregate(n~loc,data=res1[cond,],sum)
    n.maiorjl=tmp$n
    
    Vjl[,i]=rbeta(locations,njl+1,n.maiorjl+gamma)
  }
  Vjl[,community]=1
  
  #calculate Theta
  Theta[,1]=Vjl[,1]
  Theta[,2]=Vjl[,2]*(1-Vjl[,1])
  for (i in 3:community){
    Theta[,i]=Vjl[,i]*apply(1-Vjl[,1:(i-1)],1,prod)
  }
  list(Theta=Theta,Vjl=Vjl)
}
#------------------------------------------------------
update.Omega=function(param){
  Omega=matrix(NA,bands,community)
  zmat=param$zmat
  
  #calculate number of individuals in each community and band with reflectance = 1
  cond=zmat$reflect==1
  n1=aggregate(n~comm+bands,data=zmat[cond,],sum)

  #calculate number of individuals in each community and band with reflectance = 0
  cond=zmat$reflect==0
  n0=aggregate(n~comm+bands,data=zmat[cond,],sum)
  
  #make sure that we have 0's for empty communities 
  aux=expand.grid(comm=1:community,bands=1:bands)
  
  n1a=merge(n1,aux,all=T)
  cond=is.na(n1a$n)
  n1a[cond,'n']=0
  
  n0a=merge(n0,aux,all=T)
  cond=is.na(n0a$n)
  n0a[cond,'n']=0
  
  #draw from the corresponding beta distributions
  for (i in 1:community){
    tmp1=n1a[n1a$comm==i,'n']
    tmp0=n0a[n0a$comm==i,'n']
    
    Omega[,i]=rbeta(bands,tmp1+a1,tmp0+a2)
  }
  Omega
}
#----------------------------------------------
update.zmat=function(param){
  res=matrix(NA,bands*community*locations*2,5)
  oo=1
  for (b in 1:bands){
    for (l in 1:locations){
      #for individuals that reflect, calculate prob of belonging to each community
      prob=param$Omega[b,]*param$Theta[l,]
      prob1=prob/sum(prob)
      
      #draw community membership and store results
      if (dat[l,b]!=0){
        ztmp=rmultinom(1,size=dat[l,b],prob=prob1)
        ind=which(ztmp!=0)
        fim1=cbind(b,l,1,ind,ztmp[ind])
      }
      
      #for individuals that don't reflect, calculate prob of belonging to each community
      prob=(1-param$Omega[b,])*param$Theta[l,]
      prob1=prob/sum(prob)

      #draw community membership and store results
      if (max.dat-dat[l,b]!=0){
        ztmp=rmultinom(1,size=max.dat-dat[l,b],prob=prob1)
        ind=which(ztmp!=0)
        fim0=cbind(b,l,0,ind,ztmp[ind])
      }
      
      #store results
      if (dat[l,b]!=0 & max.dat-dat[l,b]!=0) fim=rbind(fim1,fim0)       
      if (dat[l,b]==0 & max.dat-dat[l,b]!=0) fim=fim0
      if (dat[l,b]!=0 & max.dat-dat[l,b]==0) fim=fim1
      
      res[oo:(oo+nrow(fim)-1),]=fim
      oo=oo+nrow(fim)
    }
  }
  res1=data.frame(res)
  colnames(res1)=c('bands','loc','reflect','comm','n')
  res1
}
#---------------------------------------------
calc.loglikel=function(param){ 
  probs=param$Theta%*%t(param$Omega)
  loglikel=dbinom(dat,size=max.dat,prob=probs,log=T)
  
  prior.omega=dbeta(param$Omega,a1,a2,log=T)
  prior.Vjl=dbeta(param$Vjl,1,gamma,log=T)
  sum(loglikel)+sum(prior.omega)+sum(prior.Vjl)
}