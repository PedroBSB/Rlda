set.seed(2929)

a=b=matrix(abs(rnorm(10000)),1000,10)
c=d=matrix(abs(rnorm(1000)),10,100)

#teste<- generateThetaBinomialVariational(1000, 10 , a, b)
teste<- generateM1M0matrix(1000,10, 100, a, b, c, d)


nloc=1000
ncommun=10
nspp=100

get.m1m0=function(nloc,ncommun,nspp,a,b,c,d){
  #get stb prior piece
  stb=digamma(a)-digamma(a+b)
  stb[,ncommun]=0 #because V_{lK}=1

  res=rep(0,nloc)
  for (i in 2:ncommun){
    res=res+digamma(b[,i-1])-digamma(a[,i-1]+b[,i-1])
    stb[,i]=stb[,i]+res
  }

  #get data part
  pdat1=digamma(c)-digamma(c+d)
  pdat0=digamma(d)-digamma(c+d)

  m0=m1=array(NA,dim=c(nloc,nspp,ncommun),
              dimnames=list(paste('loc',1:nloc,sep=''),
                            paste('spp',1:nspp,sep=''),
                            paste('comm',1:ncommun,sep='')))
  for (i in 1:nloc){
    for (j in 1:nspp){
      m1.tmp=exp(pdat1[,j]+stb[i,])
      m1[i,j,]=m1.tmp/sum(m1.tmp)

      m0.tmp=exp(pdat0[,j]+stb[i,])
      m0[i,j,]=m0.tmp/sum(m0.tmp)
    }
  }
  list(m1=m1,m0=m0)
}


teste2<-get.m1m0(1000,10, 100, a, b, c, d)

head(teste$M1)
head(teste2$m1)


all(teste$M1==teste2$m1)
