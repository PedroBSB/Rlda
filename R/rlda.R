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
  res <- Rlda:::lda_bernoulli(data, n_community,
                       alpha0, alpha1, gamma,
                       n_gibbs, ll_prior, display_progress)
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
  #Create the class
  class(res) <- c("list", "rlda")
  return(res)
}

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
  res <- Rlda:::lda_multinomial(data, n_community, beta, gamma,
                         n_gibbs, ll_prior, display_progress)
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
  #Create the class
  class(res) <- c("list", "rlda")
  return(res)
}

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
  res <- Rlda:::lda_binomial(data, pop, n_community,  alpha0 , alpha1, gamma,
                             n_gibbs, ll_prior, display_progress)
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
  #Create the class
  class(res) <- c("list", "rlda")
  return(res)
}

#' Plot the Rlda object.
#'
#' @param d rlda object
#' @param ... ignored
#' @export
plot.rlda <- function(d,...){
  #Burn-in
  i<- ceiling(d$n_gibbs*0.1)
  #Plot the log-likelihood
  plot(d$logLikelihood[i:d$n_gibbs],type="l", main="Log-Likelihood")
  par(ask=T)
  #Plot the box-plot Theta
  tmp<- colMeans(d$Theta[i:d$n_gibbs,])
  theta<- matrix(tmp,d$N,d$n_community)
  colnames(theta)=paste(1:d$n_community,sep='')
  boxplot(theta,main="Theta matrix")
  par(ask=T)
  #Plot the box-plot Phi
  tmp<- colMeans(d$Phi[i:d$n_gibbs,])
  phi<- matrix(tmp,d$n_community,length(d$Species))
  rownames(phi)=paste(1:d$n_community,sep='')
  colnames(phi)=d$Species
  boxplot(phi,main="Phi matrix")
  invisible(d)
}

#' Summarize the Bayesian LDA.
#'
#' @param d rlda object
#' @param ... ignored
#' @export
summary.rlda <-function(d,...){
  #Burn-in
  i<- ceiling(d$n_gibbs*0.1)
  seq<-i:d$n_gibbs
  print(paste("Total number of gibbs sampling:", d$n_gibbs))
  print(paste("Number of clusters:", d$n_community))
  print(paste("Number of variables:", length(d$Species)))
}







