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

#' @name rlda.bernoulliMH
#' @title Gibbs Sampling for LDA Presence and Absence with Metropolis-Hasting
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
rlda.bernoulliMH<-function(data, loc.id, n_community, alpha0, alpha1, gamma,
                         n_gibbs, nadapt, ll_prior = TRUE, display_progress = TRUE){
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
  dat<-data
  a.phi<-alpha0
  b.phi<-alpha1
  nspp<-ncol(data)
  ncomm<-n_community
  nloc <- nrow(data)
  y<-as.matrix(data)
  ngibbs<- n_gibbs

  #initial values
  #convert from a bunch of bernoulli to a single binomial per location
  tmp=aggregate.data(dat)
  y=tmp$dat
  loc.id=tmp$loc.id
  nspp=ncol(y)
  nloc=length(unique(loc.id))
  n=tmp$n
  nmat=matrix(n,nloc,nspp)

  #initial values
  theta=matrix(1/ncomm,nloc,ncomm)
  vmat=theta
  vmat[,ncomm]=1
  phi=matrix(0.2,ncomm,nspp,byrow=T)

  #gibbs stuff
  vec.theta=matrix(0,ngibbs,nloc*ncomm)

  #Remove the loc.id
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
    tmp=update.phiAbundanceSB(param=param,jump=jump1$phi,ncomm=ncomm,nspp=nspp,y=y,nmat=nmat,a.phi=a.phi,b.phi=b.phi)
    param$phi=tmp$phi
    accept1$phi=accept1$phi+tmp$accept

    tmp=update.thetaAbundanceSB(param=param,jump=jump1$vmat,nloc=nloc,ncomm=ncomm,y=y,nmat=nmat,gamma=gamma)
    param$theta=tmp$theta
    param$vmat=tmp$v
    accept1$vmat=accept1$vmat+tmp$accept

    if (i%%accept.output==0 & i<nadapt){
      k=print.adapt(parmAccept=accept1,parmJump=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }

    #to assess convergence, examine logl
    prob=get.logl(theta=param$theta,phi=param$phi,y=y,nmat=nmat)
    loglikel=sum(prob)
    if(ll_prior){
      loglikel=loglikel+sum(dbeta(param$phi,a.phi,b.phi,log=T))+
        sum(dbeta(param$vmat[,-ncomm],1,gamma,log=T))
    }

    vec.logl[i]=loglikel
    vec.theta[i,]=param$theta

    #Remove the loc.id
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
  res$N<- (length(vec.theta[n_gibbs,])/ncomm)
  #Covariates
  res$Species<- seq(1,nspp)
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
  res$rownames<-unique(loc.id)
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


#' @name rlda.binomialMH
#' @title Compute the Gibbs Sampling for LDA Binomial for Remote Sensing
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
rlda.binomialMH<-function(data, pop, n_community, alpha0 , alpha1, gamma,
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

  #Dictionary
  a.omega<-alpha0
  b.omega<-alpha1
  nbands<-ncol(data)
  ncommun<-n_community
  nloc <- nrow(data)
  ngibbs<- n_gibbs
  ndig.values<-as.matrix(pop)
  remote<-as.matrix(data)

  #initial values
  omega=matrix(runif(ncommun*nbands),ncommun,nbands)
  theta=matrix(1/ncommun,nloc,ncommun)
  v=theta
  v[,ncommun]=1

  #stuff for gibbs sampling
  param=list(theta=theta,omega=omega,v=v,gamma=gamma)
  vec.theta=matrix(NA,ngibbs,nloc*ncommun)
  vec.omega=matrix(NA,ngibbs,ncommun*nbands)

  #stuff for MH algorithm
  jump1=list(omega=matrix(1,ncommun,nbands),
             v=matrix(0.3,nloc,ncommun))
  accept1=list(omega=matrix(0,ncommun,nbands),
               v=matrix(0,nloc,ncommun))
  accept.output=50

  # create progress bar
  if(display_progress) pb <- txtProgressBar(min = 0, max = ngibbs, style = 3)

  for (i in 1:ngibbs){
    tmp=update.thetaRemote(remote, param,jump1$v,ncommun,nloc,ndig.values)
    param$theta=tmp$theta
    param$v=tmp$v
    accept1$v=accept1$v+tmp$accept

    tmp=update.omegaRemote(remote, param,jump1$omega,ncommun,nbands,ndig.values,a.omega,b.omega)
    param$omega=tmp$omega
    accept1$omega=accept1$omega+tmp$accept

    if (i%%accept.output==0 & i<1000){
      k=print.adapt(parmAccept=accept1,parmJump=jump1,accept.output=accept.output,FALSE)
      accept1=k$accept1
      jump1=k$jump1
    }

    vec.theta[i,]=param$theta
    vec.omega[i,]=param$omega
    #Progress Bar
    if(display_progress) setTxtProgressBar(pb, i)
  }
  if(display_progress) close(pb)

  vec.logl<-rep(0,nloc)

  #Use a function not exported
  # Execute the LDA for the Bernoulli entry
  res <- list("Theta"=vec.theta,
              "Phi"=vec.omega,
              "logLikelihood"=vec.logl)

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
plot.rlda <- function(x, burnin=0.1, maxCluster=NA, ...){
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin>1 || burnin<0))
  #Burn-in
  i<- ceiling(x$n_gibbs*burnin)
  #Plot the log-likelihood
  plot(x$logLikelihood[i:x$n_gibbs],type="l",xlab="Gibbs iteration", ylab="Log-Likelihood",main="Log-Likelihood")
  par(ask=T)
  #Plot the box-plot Theta
  if(is.na(maxCluster)) maxCluster=x$n_community
  tmp<- colMeans(x$Theta[i:x$n_gibbs,])
  theta<- matrix(tmp,x$N,maxCluster)
  colnames(theta)=paste("Cluster ", 1:maxCluster,sep='')
  rownames(theta)=x$rownames
  boxplot(theta,main="Theta matrix",  ylab="Probability")
  par(ask=T)
  #Plot the box-plot Phi
  tmp<- colMeans(x$Phi[i:x$n_gibbs,])
  phi<- matrix(tmp,maxCluster,length(x$Species))
  rownames(phi)=paste("Cluster ", 1:maxCluster,sep='')
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





