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
predict.rlda <-function(object, data, nclus=5, burnin=0.1, ...){
  stopifnot(inherits(object, "rlda"))
  stopifnot(inherits(nclus, "numeric"))
  stopifnot(nclus>0)
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin>1 || burnin<0))

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

  #Matrix
  summ<-summary(object,burnin,T)
  theta<-summ$Theta
  #Calculate implied binomial probabilities
  probs<- data.matrix(combo1)%*%data.matrix(t(theta))
  #Import the data for the desired region
  dat1<- data
  nbands<- length(object$Species)
  #Let's change the range of our data to start at zero.
  tmp<- apply(dat1,2,range)
  dat2<- dat1-matrix(tmp[1,],nrow(dat1),nbands,byrow=T)
  tmp<- apply(dat2,2,range)
  max1<- tmp[2,]
  max2<- matrix(max1,nrow(probs),length(max1),byrow=T)
  #Find which proportion of communities  yield the highest probability
  res<- matrix(NA,nrow(dat2),nclus)
  for (i in 1:nrow(dat2)){
   rasc<- matrix(dat2[i,],nrow(probs),ncol=ncol(dat2),byrow=T)
   llikel<- dbinom(rasc,size=max2,prob=probs,log=T)
   fim<- apply(llikel,1,sum)
   ind<- which(fim==max(fim))
   res[i,]<- as.numeric(combo1[ind,])
  }
  colnames(res)<- paste('prop',1:nclus,sep='')

  return(res)
}







