#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
// [[Rcpp::depends(RcppProgress)]]
#include "progress.hpp" // Within MacOS, one needs quote.
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

int whichLessDVSparsenessPresence(double value, NumericVector prob) {
  int res = -1;

  //'Create the categorical table
  double probcum = 0;

  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}


NumericVector rmultinomialDVSparsenessPresence(int size, NumericVector prob) {
  //'Initialize the NumericMatrix result
  NumericVector res(prob.length());
  //'Initialize the values
  for(int j=0;j<prob.length();j++){
    res(j)=0;
  }

  //'Generate the sample
  int count=0;
  //'Generate until get the sample size
  while(count<size){
    //'Draw a uniform
    double random = R::runif(0,1);
    //'Which category was draw ?
    int iPos=whichLessDVSparsenessPresence(random,prob);
    //'Increment the matrix
    res(iPos)=res(iPos)+1;
    //'Increment the counter
    count=count+1;
  }

  return res;
}

double matchSparsenessPresence(double value,NumericVector vec){
  for(int i=0;i<vec.length();i++){
    if(value==vec(i)){
      return(i);
    }
  }
  return(-1);
}

//Generate the Truncated Normal
double tnormSparseness(double lo,double hi,double mu,double sig){
      double q1 = R::pnorm(lo,mu,sig,true,false);
      double q2 = R::pnorm(hi,mu,sig,true,false);
      double z0 = R::runif(q1,q2);
      double z1  = R::qnorm(z0,mu,sig,true,false);
  return(z1);
}

NumericVector rdirichletSparseness(Rcpp::NumericVector parms) {
  NumericVector res(parms.size());
  double sample_sum = 0;
  for(int j=0;j<parms.size();j++){
    res(j) = R::rgamma(parms[j], 1);
    sample_sum += res(j);
  }
  for(int j = 0; j<parms.size();j++) {
    res(j) = res(j) / sample_sum ;
  }
  return (res);
}

NumericMatrix rdirichletSparseness(int n, Rcpp::NumericVector parms) {
  NumericMatrix res(n, parms.size());
  for(int i=0;i<n;i++){
    double sample_sum = 0;
    for(int j=0;j<parms.size();j++){
      res(i, j) = R::rgamma(parms[j], 1);
      sample_sum += res(i, j);
    }
    for(int j = 0; j<parms.size();j++) {
      res(i, j) = res(i, j) / sample_sum ;
    }
  }
  return (res);
}

double fix_MH(double lo, double hi,double old1, double new1, double jump){
  double jold=R::pnorm(hi,old1,jump,true,false)-R::pnorm(lo,old1,jump,true,false);
  double jnew=R::pnorm(hi,new1,jump,true,false)-R::pnorm(lo,new1,jump,true,false);
  return(std::log(jold)-std::log(jnew));
}

// if BLOCK, then accept as a block, otherwise, accept individually
double acceptMH(double p0,double p1,double x0,double x1){
  //acceptance PR
  double a=std::exp(p1 - p0);
  double z=R::runif(0.0,1.0);
  if(z<a){
    return(x0);
  }
  else{
    return(x1);
  }
return(-1.0);
}

NumericVector meltSparsenessBernoulli(NumericMatrix mat){
  //'Initialize the NumericVector
  NumericVector vec(mat.nrow()*mat.ncol());
  //'Initialize the position
  int pos=0;
  for(int col=0;col<mat.ncol();col++){
    for(int row=0;row<mat.nrow();row++){
      //'meltPresence the matrix
      vec(pos)=mat(row,col);
      //'Increment the position
      pos=pos+1;
    }
  }
  return(vec);
}
void updateThetaAndPhiSparsenessBernoulli(NumericMatrix &ThetaGibbs,NumericMatrix Theta,NumericMatrix &PhiGibbs, NumericMatrix Phi,int gibbs){
  //'meltPresence the Theta and Phi matrix
  ThetaGibbs(gibbs,_)=meltSparsenessBernoulli(Theta);
  PhiGibbs(gibbs,_) =meltSparsenessBernoulli(Phi);
}


/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/


List generateZPresenceSparseness(NumericMatrix binaryMat, NumericMatrix Theta, NumericMatrix Phi) {
  //'Number of locations
  int nLocations=binaryMat.nrow();
  //'Number of Species
  int nSpecies=binaryMat.ncol();
  //'Number of Communities
  int n_community=Theta.ncol();
  //Create the Zmat
  NumericMatrix zMat(nLocations, nSpecies);
  //'Matrix N the number of locations where specie s comes from community c
  NumericMatrix nMat(nSpecies,n_community);
  //'Matrix R the number of these individuals thar are observed
  NumericMatrix rMat(nSpecies,n_community);
  //'Matrix R the number of these individuals thar are observed
  NumericMatrix lMat(nLocations,n_community);
  lMat.fill(0);
  zMat.fill(-1);
  for(int l=0;l<nLocations;l++){
     for(int s=0;s<nSpecies;s++){
      //'Calculate the probability vector
      NumericVector probability = pow(Theta(s,_),binaryMat(l,s))*pow(1.0-Theta(s,_),1.0-binaryMat(l,s))*Phi(l,_);
      //'Find the sum
      double sumVec=sum(probability);
      //'Normalize the probability
      probability=probability/sumVec;
      //'Presence in locatiol l and specie s (always will be a draw)
      int iSize=1;
      //'Store the results
      NumericVector tmp = rmultinomialDVSparsenessPresence(iSize, probability);
      //'Find the community
      int iCommunity = (int) matchSparsenessPresence(1.0,tmp);
      zMat(l,s)=iCommunity;
      //'N matrix is the number of individuals in plot L that come from community C
      nMat(s,iCommunity)=nMat(s,iCommunity)+1.0;
      //'R matrix is the number of individuals in plot L that come from community C that are observed
      if(binaryMat(l,s)==1.0) rMat(s,iCommunity)=rMat(s,iCommunity)+1.0;
      //Update lMat
      lMat(l,iCommunity)=lMat(l,iCommunity)+1;
    }
  }

  //'Store the zMat
  List resZ = Rcpp::List::create(Rcpp::Named("N") = nMat,
                                 Rcpp::Named("R") = rMat,
                                 Rcpp::Named("Z") = zMat,
                                 Rcpp::Named("L") = lMat);

  return resZ;
}


NumericMatrix generatePhi(List zList,NumericVector gammaVec, int nCommunity){
  //'Getting the L matrix
  NumericMatrix lMat = zList[3];
  //Create the Phi matrix
  NumericMatrix Phi(lMat.nrow(),nCommunity);
  //For each location
  for(int l=0;l<lMat.nrow();l++){
    NumericVector parms = lMat(l,_)+gammaVec;
    Phi(l,_) = rdirichletSparseness(parms);
  }
  return(Phi);
}


NumericMatrix generateThetaSparseness(List zList,double alpha0, double alpha1) {
  //'Getting the N matrix
  NumericMatrix nMat = zList[0];
  //'Getting the R matrix
  NumericMatrix rMat = zList[1];
  //'Total number of species
  int nSpecies = nMat.nrow();
  //'Total number of communities
  int n_community = nMat.ncol();
  //'Initialize the Phi matrix
  NumericMatrix thetaMat(nSpecies,n_community);
  //'For each Specie
  for(int s=0;s<nSpecies;s++){
    //'For each community:
    for(int c=0;c<n_community;c++){
      double aBeta = rMat(s,c)+alpha0;
      double bBeta = nMat(s,c)-rMat(s,c)+alpha1;
      thetaMat(s,c)=R::rbeta(aBeta,bBeta);
    }
  }
  return thetaMat;
}


double poposedGammaPDF(NumericMatrix lMat,NumericVector gammaVec){
  //'Number of locations
  int nLocations=lMat.nrow();
  //Calculate the constant
  double sumGamma = Rcpp::sum(gammaVec);
  //Part 1:
  double part1 = nLocations*R::lgammafn(sumGamma);
  double part2 = (-1.0)*nLocations*Rcpp::sum(Rcpp::lgamma(gammaVec));

  double part3 = 0;
  double part4 = 0;
  for(int l=0;l<lMat.nrow();l++){
    part3 = part3+ Rcpp::sum(lMat(l,_)+gammaVec);
    part4 = part4+ (-1.0)* R::lgammafn(Rcpp::sum(lMat(l,_))+sumGamma);
  }
  return(part1+part2+part3+part4);
}

NumericVector generateGamma(List zList, NumericVector jump, NumericVector& acceptVector) {
  //'Getting the L matrix
  NumericMatrix lMat = zList[3];
  //'Total number of communities
  int n_community = lMat.ncol();
  //'Number of locations
  int nLocations=lMat.nrow();
  //'Getting the gamma vector
  NumericVector gammaVec(n_community);
  gammaVec.fill(1.0);
  //Initialize the Metropolis-Hasting:
  for(int c=0;c<n_community;c++){
    //Store the old gamma
    NumericVector gammaNewVec = gammaVec;
    //Generate the Truncated Normal
    double res = tnormSparseness(0.0,1.0,gammaVec(c),jump(c));
    gammaNewVec(c) = res;

    //Calculate the probability
    double pold=poposedGammaPDF(lMat,gammaVec);
    double pnew=poposedGammaPDF(lMat,gammaNewVec);

    //Correction based on Metropolis-Hasting
    double probCorrection = fix_MH(0.0,1.0,gammaVec(c),gammaNewVec(c),jump(c));
    double k = acceptMH(pold,pnew+probCorrection,gammaVec(c),gammaNewVec(c));
    if(gammaVec(c) != k){
      gammaVec(c) = k;
      acceptVector(c) = 1.0;
    }
  }

  return gammaVec;
}


//TODO: Check the prior for Phi. How can we add this prior into the log-likelihood
double ll_priorFunctionSparsnessBernoulli(NumericMatrix matDATA,int nLocations, int nSpecies,int n_community, NumericMatrix Theta, NumericMatrix Phi, double alpha0,double alpha1, NumericVector gammaVec, bool ll_prior=true) {
  //'Initialize the logLikelihoodVec
  double logLikelihood=0;
  //'Calculate the Loglikelihood and Prior
  if(ll_prior){
    //'Initialize the Theta_{sc} prior
    double priorTheta=0.0;
    //'For each location
    for(int l=0;l<nLocations;l++){
      //'Initiate the Theta Counter
      int thetaGibbsCount=0;

      //'For each Specie
      for(int s=0;s<nSpecies;s++){
        //'Initiate the Gibbs Counter
        int phiGibbsCount=0;

        //'Compute just one time the prior for Theta_{sc}
        if(l==0){
          //'Compute the prior for Theta_{sc}
          for(int c=0;c<n_community;c++){
            if(Theta(s,c)<1)priorTheta=priorTheta+R::dbeta(Theta(s,c),alpha0,alpha1,1);
          }
        }

        //'Initialize the product between Theta and Phi
        double thetaPhi=0.0;
        //'For each community
        for(int c=0;c<n_community;c++){
          thetaPhi=thetaPhi+Theta(s,c)*Phi(l,c);  //'DV: mudei aqui
        }
        //'Getting the x_{ls} observation
        double xLS=matDATA(l,s);
        if(xLS==1.0){
          logLikelihood=logLikelihood+log(thetaPhi);
        }
        else{
          if(thetaPhi<1){//'Otherwise a large thetaPhi>1 represents low likelihood, then we don't compute this value.
            logLikelihood=logLikelihood+log(1.0-thetaPhi);
          }
        }
      }
    }

    //'Update the logLikelihood
    logLikelihood=logLikelihood+priorTheta;
  }
  else{
    //'For each location
    for(int l=0;l<nLocations;l++){
      //'Initiate the Theta Counter
      int thetaGibbsCount=0;
      //'For each Specie
      for(int s=0;s<nSpecies;s++){
        //'Initiate the Gibbs Counter
        int phiGibbsCount=0;
        //'Initialize the product between Theta and Phi
        double thetaPhi=0.0;
        //'For each community
        for(int c=0;c<n_community;c++){
          thetaPhi=thetaPhi+Theta(s,c)*Phi(l,c); //'DV: mudei aqui
        }
        //'Getting the x_{ls} observation
        double xLS=matDATA(l,s);
        if(xLS==1.0){
          logLikelihood=logLikelihood+log(thetaPhi);
        }
        else{
          if(thetaPhi<1){//'Otherwise a large thetaPhi>1 represents low likelihood, then we don't compute this value.
            logLikelihood=logLikelihood+log(1.0-thetaPhi);
          }
        }
      }
    }
  }

  return(logLikelihood);
}

/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING PROCEDURE                  ************************************/
/***************************************************************************************************************************/

//' @name GibbsSamplingPresence
//' @title Gibbs Sampling for LDA Presence and Absence
//' @description Compute the Gibbs Sampling for LDA Presence and Absence
//' @param DATA - DataFrame with Presence and Absecence (Zeros and Ones)
//' @param int n_community - Number of communities
//' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
//' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
//' @param n_gibbs - Total number of Gibbs Samples
//' @param ll_prior - Likelihood compute with Priors ?
//' @param bool display_progress=true - Should I Show the progressBar ?
//' @return List - With Theta(n_gibbs,n_community*nSpecies), Phi(n_gibbs,nLocations*n_community) and logLikelihood
// [[Rcpp::export]]
List lda_bernoulli_sparseness(DataFrame data, int n_community, double alpha0, double alpha1, int n_gibbs, bool ll_prior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix matDATA = internal::convert_using_rfunction(data, "as.matrix");

  //'Total number of locations
  int nLocations = matDATA.nrow();

  //'Total number of species
  int nSpecies = matDATA.ncol();

  //'Intialize Theta
  NumericVector hyperTheta(n_community);
  hyperTheta.fill(1);
  NumericMatrix Theta=rdirichletSparseness(nSpecies,hyperTheta);

  //'Intialize Phi
  NumericVector hyperPhi(n_community);
  hyperPhi.fill(1);
  NumericMatrix Phi=rdirichletSparseness(nLocations,hyperPhi);

  //'Intialize gammaVec
  NumericVector gammaVec(n_community);
  gammaVec.fill(1);

  //Initialize the acceptVector
  NumericVector acceptVector(n_community);
  acceptVector.fill(0);

  //Initialize the Jump
  NumericVector jumpVec(n_community);
  jumpVec.fill(0.5);

  //'Initialize the ThetaGibbs
  NumericMatrix ThetaGibbs(n_gibbs,nSpecies*n_community);

  //'Initialize the PhiGibbs
  NumericMatrix PhiGibbs(n_gibbs,nLocations*n_community);

  //'Initialize the logLikelihood vector
  NumericVector logLikelihoodVec(n_gibbs);

  //'Intialize the progressbar
  Progress p(n_gibbs, display_progress);
  for (int g = 0; g < n_gibbs; ++g) {
    //'Verify if everything is ok
    if (Progress::check_abort() )
    Rcpp::stop("Operation cancelled by interrupt.");

    //'Generate zList
    List zList  = generateZPresenceSparseness(matDATA, Theta, Phi);

    //'Generate Phi
    Phi = generatePhi(zList, gammaVec, n_community);

    //'Generate Theta
    Theta = generateThetaSparseness(zList,alpha0,alpha1);

    //'Getting the N matrix
    NumericMatrix nMat = zList[0];
    if(g%50==0 & g<300){
      NumericVector z = acceptVector/50;
      for(int t=0;t<z.size();t++){
        if(z(t)>0.4 & jumpVec(t)<100){
          jumpVec(t)=jumpVec(t)*2;
        }
        else if(z(t)<0.1 & jumpVec(t)>0.001){
          jumpVec(t)=jumpVec(t)*0.5;
        }
      }
      acceptVector.fill(0.0);
    }

    //Generate gamma
    gammaVec = generateGamma(zList, jumpVec, acceptVector);

    //'Create the final Theta (n_gibbs,n_community*nSpecies) and final Phi (PhiGibbs)
    updateThetaAndPhiSparsenessBernoulli(ThetaGibbs, Theta, PhiGibbs, Phi, g);

    //'Initialize the logLikelihood
    double logLikelihood=ll_priorFunctionSparsnessBernoulli(matDATA, nLocations,
                                                  nSpecies,n_community,
                                                  Theta, Phi,
                                                  alpha0, alpha1, gammaVec,
                                                  ll_prior);
    //'Store the logLikelihood
    logLikelihoodVec(g)=logLikelihood;

    //'Increment the progress bar
    p.increment();

  }

  //'Store the results - Order change to agree with the other functions
  List resTemp = Rcpp::List::create(Rcpp::Named("Theta") = PhiGibbs,
                                    Rcpp::Named("Phi")  = ThetaGibbs,
                                    Rcpp::Named("logLikelihood")  =logLikelihoodVec);

  return resTemp;
}



