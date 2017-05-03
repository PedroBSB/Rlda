#' @name rlda.bernoulli
#' @title Gibbs Sampling for LDA Presence and Absence
#' @description Compute the Gibbs Sampling for LDA Presence and Absence
#' @param data - DataFrame with Presence and Absecence (Zeros and Ones)
#' @param n_community - Number of communities
#' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
#' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.bernoulli<-function(data, n_community, alpha0, alpha1, gamma,
                         n_gibbs, ll_prior = TRUE, display_progress = TRUE){
  #Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(alpha0, "numeric"))
  stopifnot(inherits(alpha1, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))

  #Use a function not exported
  # Execute the LDA for the Bernoulli entry
  res <- lda_bernoulli(data, n_community,
            alpha0, alpha1, gamma,
            n_gibbs, ll_prior, display_progress)
  #Type distribution
  res$type<- "Bernoulli"
  #Number of communities
  res$n_community<- n_community
  #Sample size
  res$N<- nrow(data)
  #Covariates
  res$Species<- colnames(data)
  #Alpha0
  res$alpha0<- alpha0
  #Alpha1
  res$alpha1<- alpha1
  #Gamma
  res$gamma<- gamma
  #Number of gibbs
  res$n_gibbs<- n_gibbs
  #Species
  res$colnames<-colnames(data)
  #Locations
  res$rownames<-rownames(data)
  #Create the class
  class(res) <- c("list", "rlda")
  return(res)
}

#' @name rlda.bernoulliSB
#' @title Gibbs Sampling for LDA Presence and Absence with  Stick-Breaking Priors
#' @description Compute the Gibbs Sampling for LDA Presence and Absence with Stick-Breaking Priors
#' @param data - DataFrame with Presence and Absecence (Zeros and Ones)
#' @param n_community - Number of communities
#' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
#' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.bernoulliSB<-function(data, n_community, alpha0, alpha1, gamma,
                         n_gibbs, ll_prior = TRUE, display_progress = TRUE){
  #Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(alpha0, "numeric"))
  stopifnot(inherits(alpha1, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))

  #Dictionary
  a.phi<-alpha0
  b.phi<-alpha1
  loc.id<-seq(1,nrow(data))
  nspp<-ncol(data)
  ncomm<-n_community
  nloc <- nrow(data)
  y<-as.matrix(data)
  ngibbs<- n_gibbs

  ############################################################################################

  print.adapt = function(accept1z,jump1z){
    accept1=accept1z; jump1=jump1z;

    for (k in 1:length(accept1)){
      z=accept1[[k]]/accept.output
    }

    for (k in 1:length(jump1)){
      cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<10000
      jump1[[k]][cond] = jump1[[k]][cond]*2
      cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
      jump1[[k]][cond] = jump1[[k]][cond]*0.5
      accept1[[k]][]=0
    }

    return(list(jump1=jump1,accept1=accept1))
  }
  #----------------------------
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
  rmvnorm1=function (n, sigma, pre0.9_9994 = FALSE)
  {
    #   retval <- chol(sigma, pivot = TRUE)
    #   o <- order(attr(retval, "pivot"))
    #   retval <- retval[, o]
    s. <- svd(sigma)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    R = t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
    retval
  }
  #----------------------------
  fix.probs=function(probs){
    cond=probs<0.00001
    probs[cond]=0.00001
    cond=probs>0.99999
    probs[cond]=0.99999
    probs
  }
  #----------------------------
  update.phi=function(param,jump){
    phi.orig=phi.old=param$phi
    proposed=matrix(tnorm(nspp*ncomm,lo=0,hi=1,mu=phi.old,sig=jump),ncomm,nspp)

    for (i in 1:ncomm){
      phi.new=phi.old
      phi.new[i,]=proposed[i,]
      adj=fix.MH(lo=0,hi=1,phi.old[i,],phi.new[i,],jump[i,])

      prob.old=fix.probs(param$theta%*%phi.old)[loc.id,]
      prob.new=fix.probs(param$theta%*%phi.new)[loc.id,]
      pold=colSums(dbinom(y,size=1,prob=prob.old,log=T))+dbeta(phi.old[i,],a.phi,b.phi,log=T)
      pnew=colSums(dbinom(y,size=1,prob=prob.new,log=T))+dbeta(phi.new[i,],a.phi,b.phi,log=T)
      k=acceptMH(pold,pnew+adj,phi.old[i,],phi.new[i,],F)
      phi.old[i,]=k$x
    }
    list(phi=phi.old,accept=phi.orig!=phi.old)
  }
  #-----------------------------
  update.theta=function(param,jump){
    v.orig=v.old=param$vmat
    tmp=tnorm(nloc*(ncomm-1),lo=0,hi=1,mu=v.old[,-ncomm],sig=jump[,-ncomm])
    novos=cbind(matrix(tmp,nloc,ncomm-1),1)
    ajuste=matrix(fix.MH(lo=0,hi=1,v.old,novos,jump),nloc,ncomm)

    prior.old=matrix(dbeta(v.old,1,gamma,log=T),nloc,ncomm)
    prior.new=matrix(dbeta(novos,1,gamma,log=T),nloc,ncomm)

    for (j in 1:(ncomm-1)){ #last column has to be 1
      v.new=v.old
      v.new[,j]=novos[,j]

      theta.old=convertSBtoNormal(vmat=v.old,ncol=ncomm,nrow=nloc,prod=rep(1,nloc))
      theta.new=convertSBtoNormal(vmat=v.new,ncol=ncomm,nrow=nloc,prod=rep(1,nloc))

      #contribution from reflectance data
      pold=fix.probs(theta.old%*%param$phi)[loc.id,]
      pnew=fix.probs(theta.new%*%param$phi)[loc.id,]

      k=rowSums(dbinom(y,size=1,prob=pold,log=T))
      p1.old=aggregatesum(Tobesum=k, nind=nloc, nobs=length(k), ind=loc.id)

      k=rowSums(dbinom(y,size=1,prob=pnew,log=T))
      p1.new=aggregatesum(Tobesum=k, nind=nloc, nobs=length(k), ind=loc.id)

      k=acceptMH(p1.old+prior.old[,j],
                 p1.new+prior.new[,j]+ajuste[,j],
                 v.old[,j],v.new[,j],F)
      v.old[,j]=k$x
    }
    theta=convertSBtoNormal(vmat=v.old,ncol=ncomm,nrow=nloc,prod=rep(1,nloc))
    list(theta=theta,v=v.old,accept=v.old!=v.orig)
  }
  ############################################################################################
  #initial values
  theta=matrix(1/ncomm,nloc,ncomm)
  vmat=theta
  vmat[,ncomm]=1
  phi=matrix(colMeans(y),ncomm,nspp,byrow=T)
  phi[phi==1]=0.99999

  #Initial Values Gibbs
  vec.theta=matrix(0,ngibbs,nloc*ncomm)
  vec.phi=matrix(0,ngibbs,ncomm*nspp)
  vec.logl=matrix(NA,ngibbs,1)
  param=list(theta=theta,phi=phi,vmat=vmat)

  #for MH algorithm
  jump1=list(vmat=matrix(0.1,nloc,ncomm),
             phi=matrix(0.1,ncomm,nspp))
  accept1=list(vmat=matrix(0,nloc,ncomm),
               phi=matrix(0,ncomm,nspp))
  accept.output=50

  # create progress bar
  if(display_progress) pb <- txtProgressBar(min = 0, max = n_gibbs, style = 3)
  count=0
  for (i in 1:ngibbs){
    tmp=update.phi(param,jump1$phi)
    param$phi=tmp$phi #phi.true #
    accept1$phi=accept1$phi+tmp$accept

    tmp=update.theta(param,jump1$vmat)
    param$theta=tmp$theta
    param$vmat=tmp$v
    accept1$vmat=accept1$vmat+tmp$accept

    if (i%%accept.output==0 & i<1000){
      k=print.adapt(accept1,jump1)
      accept1=k$accept1
      jump1=k$jump1
    }

    #to assess convergence, examine logl
    prob=fix.probs(param$theta%*%param$phi)[loc.id,]
    loglikel <- sum(dbinom(y,size=1,prob=prob,log=T))

    #Numerical Correction
    param$phi[param$phi==0]<- 1e-10
    param$vmat[param$vmat==0]<- 1e-10


    if(ll_prior){
      loglikel <- sum(dbeta(param$phi,alpha0,alpha1,log=T))+
        sum(dbeta(param$vmat[,-n_community],1,gamma,log=T))
    }

    vec.logl[i]=loglikel
    vec.theta[i,]=param$theta
    vec.phi[i,]=param$phi

    #Progress Bar
    if(display_progress) setTxtProgressBar(pb, i)
  }
  if(display_progress) close(pb)

  #Use a function not exported
  # Execute the LDA for the Bernoulli entry
  res <- list("Theta"=vec.theta,
              "Phi"=vec.phi,
              "logLikelihood"=vec.logl)

  #Type distribution
  res$type<- "Bernoulli"
  #Number of communities
  res$n_community<- n_community
  #Sample size
  res$N<- nrow(data)
  #Covariates
  res$Species<- colnames(data)
  #Alpha0
  res$alpha0<- alpha0
  #Alpha1
  res$alpha1<- alpha1
  #Gamma
  res$gamma<- gamma
  #Number of gibbs
  res$n_gibbs<- n_gibbs
  #Species
  res$colnames<-colnames(data)
  #Locations
  res$rownames<-rownames(data)
  #Create the class
  class(res) <- c("list", "rlda")
  return(res)
}






#' @name rlda.multinomial
#' @title Gibbs Sampling for LDA Abundance with Stick-Breaking
#' @description Compute the Gibbs Sampling for LDA Abundance with Stick-Breaking
#' @param data - dataFrame with Abundance
#' @param int n_community - Number of communities
#' @param beta - NumericVector for beta (Sx1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.multinomial<-function(data, n_community, beta, gamma,
                           n_gibbs, ll_prior = TRUE, display_progress = TRUE){
  #Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(beta, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))

  #Use a function not exported
  # Execute the LDA for the Multinomial entry
  res <- lda_multinomial(data, n_community, beta, gamma,
           n_gibbs, ll_prior, display_progress)
  #Type distribution
  res$type<- "Multinomial"
  #Number of communities
  res$n_community<- n_community
  #Sample size
  res$N<- nrow(data)
  #Covariates
  res$Species<- colnames(data)
  #Beta
  res$beta<- beta
  #Gamma
  res$gamma<- gamma
  #Number of gibbs
  res$n_gibbs<- n_gibbs
  #Species
  res$colnames<-colnames(data)
  #Locations
  res$rownames<-rownames(data)
  #Create the class
  class(res) <- c("list", "rlda")
  return(res)
}

#' @name rlda.binomial
#' @title Compute the Gibbs Sampling for LDA Binomial
#' @description Compute the Gibbs Sampling for LDA Binomial
#' @param DATA - DataFrame with Presence and Absecence (Binomial)
#' @param POP - DataFrame with Population Size (Binomial)
#' @param int n_community - Number of communities
#' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
#' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.binomial<-function(data, pop, n_community, alpha0 , alpha1, gamma,
                           n_gibbs, ll_prior = TRUE, display_progress = TRUE){
  #Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(pop, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(alpha0, "numeric"))
  stopifnot(inherits(alpha1, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))
  if(nrow(data)!=nrow(pop)){
    stop('Both "data" and "pop" must have the same number of rows.')
  }
  # Execute the LDA for the Binomial entry
  res <- lda_binomial(data, pop, n_community,  alpha0 , alpha1, gamma,
                      n_gibbs, ll_prior, display_progress)
  #Type distribution
  res$type<- "Binomial"
  #Maximum value
  res$max<-max(pop)
  #Number of communities
  res$n_community<- n_community
  #Sample size
  res$N<- nrow(data)
  #Covariates
  res$Species<- colnames(data)
  #Alpha0
  res$alpha0<- alpha0
  #Alpha1
  res$alpha1<- alpha1
  #Gamma
  res$gamma<- gamma
  #Number of gibbs
  res$n_gibbs<- n_gibbs
  #Locations
  res$rownames<-rownames(data)
  #Create the class
  class(res) <- c("list", "rlda")
  return(res)
}

#' Plot the Rlda object.
#'
#' @param x rlda object
#' @param ... ignored
#' @export
plot.rlda <- function(x, burnin=0.1, ...){
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin>1 || burnin<0))
  #Burn-in
  i<- ceiling(x$n_gibbs*burnin)
  #Plot the log-likelihood
  plot(x$logLikelihood[i:x$n_gibbs],type="l",xlab="Gibbs iteration", ylab="Log-Likelihood",main="Log-Likelihood")
  par(ask=T)
  #Plot the box-plot Theta
  tmp<- colMeans(x$Theta[i:x$n_gibbs,])
  theta<- matrix(tmp,x$N,x$n_community)
  colnames(theta)=paste("Cluster ", 1:x$n_community,sep='')
  rownames(theta)=x$rownames
  boxplot(theta,main="Theta matrix",  ylab="Probability")
  par(ask=T)
  #Plot the box-plot Phi
  tmp<- colMeans(x$Phi[i:x$n_gibbs,])
  phi<- matrix(tmp,x$n_community,length(x$Species))
  rownames(phi)=paste("Cluster ", 1:x$n_community,sep='')
  colnames(phi)=x$Species
  barplot(phi, main="Phi matrix", ylab="Probability", legend=rownames(phi))
  invisible(x)
}

#' Summarize the Bayesian LDA.
#'
#' @param object rlda object
#' @param ... ignored
#' @export
summary.rlda <-function(object, burnin=0.1, silent=FALSE, ...){
  stopifnot(inherits(object, "rlda"))
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin>1 || burnin<0))
  stopifnot(inherits(silent, "logical"))

  #Burn-in
  i<- ceiling(object$n_gibbs*burnin)
  seq<-i:object$n_gibbs
  if(!silent){
    cat(paste("Total number of gibbs sampling:", object$n_gibbs,
                "\nNumber of clusters:", object$n_community,
                "\nNumber of variables:", length(object$Species)))
  }
  #Summary Theta
  tmp<- colMeans(object$Theta[i:object$n_gibbs,])
  theta<- matrix(tmp,object$N,object$n_community)
  colnames(theta)=paste("Cluster ", 1:object$n_community,sep='')
  rownames(theta)=object$rownames
  #Summary Phi
  tmp<- colMeans(object$Phi[i:object$n_gibbs,])
  phi<- matrix(tmp,object$n_community,length(object$Species))
  rownames(phi)=paste("Cluster ", 1:object$n_community,sep='')
  colnames(phi)=object$Species

  return(list("Theta"=theta,"Phi"=phi))
}


#' Predict the Bayesian LDA.
#'
#' @param object rlda object
#' @param ... ignored
#' @export
predict.rlda <-function(object, data, nclus=NA, burnin=0.1, places.round=0, ...){
  stopifnot(inherits(object, "rlda"))
  stopifnot(inherits(nclus, "numeric"))
  stopifnot(inherits(places.round, "numeric"))
  stopifnot(nclus>0 || is.na(nclus))
  stopifnot(places.round>=0)
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin>1 || burnin<0))

  #Matrix
  summ<-summary.rlda(object,burnin,T)
  phi<-summ$Phi

  if(is.na(nclus)){
    nclus<-nrow(phi)
  }

  #Create a matrix with all possible combinations of proportions
  seq1<- seq(from=0, to=1, by=0.05)
  combo<- expand.grid(p1=seq1)
  for(i in 2:(nclus-1)){
    temp<- expand.grid(var=seq1)
    colnames(temp)<-paste0("p",i)
    combo<-merge(combo,temp)
  }
  cond<- apply(combo, 1, sum)<= 1
  combo1<- combo[cond, ]
  combo1[ ,paste0("p",nclus)]<- 1-apply(combo1, 1, sum)


  #Calculate implied binomial probabilities
  probs<- data.matrix(combo1)%*%data.matrix(phi)

  #Import the data for the desired region
  dat1<- data[,object$Species]
  nbands<- length(object$Species)

  #Let's change the range of our data to start at zero.
  tmp<- apply(dat1,2,range)
  dat2<- dat1-matrix(tmp[1,],nrow(dat1),nbands,byrow=T)
  tmp<- apply(dat2,2,range)
  max1<- tmp[2,]
  max2<- matrix(max1,nrow(probs),length(max1),byrow=T)

  #Divisor
  div<-10^(places.round)
  max2<-floor(max2/div)
  dat2<-floor(dat2/div)
  df_args <- c(as.data.frame(dat2), sep="")
  dat2full<-as.data.frame(dat2)
  dat2full$ID<-as.character(do.call(paste, df_args))
  dat2full$Sort<-seq(1,nrow(dat2full))
  dat2<-unique(dat2)

  #Keep the same scale
  max2<-max2*div
  dat2<-dat2*div

  ncl<- detectCores()
  cl <- makeCluster(ncl)
  registerDoParallel(cl)
  #find which proportion of endmembers that yields the highest likelihood
  res2 <- foreach(i=1:nrow(dat2), .combine=rbind) %dopar% {
    rasc=matrix(as.integer(dat2[i,]),nrow=nrow(probs),ncol=ncol(dat2),byrow=T)
    llikel=dbinom(rasc,size=max2,prob=probs,log=T)
    fim=apply(llikel,1,sum)
    ind=which(fim==max(fim))
    as.numeric(combo1[ind,])
  }
  colnames(res2)=paste('prop',1:nclus,sep='')
  rownames(res2)=NULL
  #Stop clusters
  stopCluster(cl)

  #Convert to data.frame
  dat2<-as.data.frame(dat2)
  df_args <- c(dat2, sep="")
  res2<-as.data.frame(res2)
  res2$ID<-as.character(do.call(paste, df_args))
  final<-merge(dat2full,res2,by="ID",all=T)
  final<- final[ , -which(names(final) %in% c("ID"))]
  final <- final[order(final$Sort),]
  final<- final[, paste('prop',1:nclus,sep='')]
  return(final)
}





