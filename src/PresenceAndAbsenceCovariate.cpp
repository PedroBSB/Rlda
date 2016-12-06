#include <RcppArmadillo.h>
#include <iostream>
#include <ctime>
#include <fstream>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "progress.hpp" // Within MacOS, one needs quote.
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/



double matchPresenceCova(double value,NumericVector vec){
  for(int i=0;i<vec.length();i++){
    if(value==vec(i)){
      return(i);
    }
  }
  return(-1);
}

NumericVector rdirichletPresenceCova(Rcpp::NumericVector parms) {
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

NumericMatrix rdirichletPresenceCova(int n, Rcpp::NumericVector parms) {
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


int whichLessDVPresenceCova(double value, NumericVector prob) {
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

NumericVector rmultinomialDVPresenceCova(int size, NumericVector prob) {
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
    int iPos=whichLessDVPresenceCova(random,prob);
    //'Increment the matrix
    res(iPos)=res(iPos)+1;
    //'Increment the counter
    count=count+1;
  }

  return res;
}


NumericVector invertedCumsumPresenceCova(NumericVector n){
  NumericVector table(n.length());
  table(n.length()-1)=n(n.length()-1);
  for(int i=(n.length()-2);i>-1;i--){
    table(i)=table(i+1)+n(i);
  }
  return(table);
}

NumericVector meltPresenceCova(NumericMatrix mat){
  //'Initialize the NumericVector
  NumericVector vec(mat.nrow()*mat.ncol());
  //'Initialize the position
  int pos=0;
  for(int col=0;col<mat.ncol();col++){
    for(int row=0;row<mat.nrow();row++){
      //'meltPresenceCova the matrix
      vec(pos)=mat(row,col);
      //'Increment the position
      pos=pos+1;
    }
  }
  return(vec);
}

void updateThetaAndPhiPresenceCova(NumericMatrix &ThetaGibbs,NumericMatrix Theta,NumericMatrix &PhiGibbs, NumericMatrix Phi,NumericMatrix &BetaGibbs, NumericMatrix Beta,int gibbs){
    //'meltPresenceCova the Theta and Phi matrix
    ThetaGibbs(gibbs,_)=meltPresenceCova(Theta);
    PhiGibbs(gibbs,_) =meltPresenceCova(Phi);
    BetaGibbs(gibbs,_) =meltPresenceCova(Beta);
}


double tnormPresence(double lo, double hi,double mu, double sig){
  double q1 = R::pnorm5(lo,mu,sig,1,0);
  double q2 = R::pnorm5(hi,mu,sig,1,0);
  double z = R::runif(q1,q2);
  z = R::qnorm5(z,mu,sig,1,0);
  return(z);
}

double fixMHPresence(double lo, double hi,double old1,double new1,double jump){
  double jold=R::pnorm5(hi,old1,jump,1,0)-R::pnorm5(lo,old1,jump,1,0);
  double jnew=R::pnorm5(hi,new1,jump,1,0)-R::pnorm5(lo,new1,jump,1,0);
  return(std::log(jold)-std::log(jnew));
}

/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/

double gammaMHPresence(NumericMatrix vMat, double gamma, double jump){
  double newGamma = tnormBinomial(0.0,1.0,gamma,jump);
  double pold = 0.0;
  double pnew = 0.0;
  for(int r=0;r<vMat.nrow();r++){
    for(int c=0;c<(vMat.ncol()-1);c++){
      pold=pold+R::dbeta(vMat(r,c),1.0,gamma,1);
      pnew=pnew+R::dbeta(vMat(r,c),1.0,newGamma,1);
    }
  }
  double pcorrection = fixMHBinomial(0,1,gamma,newGamma,jump);
  //Acceptance
  double a = std::exp(pnew+pcorrection-pold);
  double z = R::unif_rand();
  if(z<a){
    return(newGamma);
  }
  else{
    return(gamma);
  }
  return(0);
}


List generateZPresenceCova(NumericMatrix binaryMat, NumericMatrix Theta, NumericMatrix vMat, NumericMatrix &PhiMat) {
  //'Number of locations
  int nLocations=binaryMat.nrow();
  //'Number of Species
  int nSpecies=binaryMat.ncol();
  //'Number of Communities
  int n_community=Theta.ncol();
  //'Matrix N the number of locations where specie s comes from community c
  NumericMatrix nMat(nSpecies,n_community);
  //'Matrix R the number of these individuals thar are observed
  NumericMatrix rMat(nSpecies,n_community);
  //'Matrix M the the number of species in plot l that come from community c
  NumericMatrix mMat(binaryMat.nrow(),n_community);
  //'For each location sample from a Multinomial
  for(int l=0;l<nLocations;l++){
    //'Create Phi (Size Cx1) based on vMat V_cl \prod_(k=1)^(c-1)(1-V_kl )
    NumericVector Phi(n_community);
    //'Update the Phi \prod_(k=1)^(c-1)(1-V_kl )
    double prod=1;
    for(int c=0;c<n_community;c++){
      double vNumber = vMat(l,c);
      if (c == 0) prod=1;
      if (c >  0) prod=prod*(1.0-vMat(l,c-1));
      Phi(c)=vNumber*prod;
    }
    //'Update the Phi Matrix by reference
    PhiMat(l,_)=Phi;
    //'For each Specie
    for(int s=0;s<nSpecies;s++){
      //'Calculate the probability vector
      double sumVec=0;
      NumericVector probability = pow(Theta(s,_),binaryMat(l,s))*pow(1.0-Theta(s,_),1.0-binaryMat(l,s))*Phi;
      //'Find the sum
      sumVec=sum(probability);
      //'Normalize the probability
      probability=probability/sumVec;
      //'PresenceCova in locatiol l and specie s (always will be a draw)
      int iSize=1;
      //'Store the results
      NumericVector tmp = rmultinomialDVPresenceCova(iSize, probability);
      //'Find the community
      int iCommunity = (int) matchPresenceCova(1.0,tmp);
      //'N matrix is the number of individuals in plot L that come from community C
      nMat(s,iCommunity)=nMat(s,iCommunity)+1.0;
      //'R matrix is the number of individuals in plot L that come from community C that are observed
      if(binaryMat(l,s)==1.0) rMat(s,iCommunity)=rMat(s,iCommunity)+1.0;
      //'M matrix is the number of species in plot L that come from community C
      mMat(l,iCommunity)=mMat(l,iCommunity)+1.0;
    }
  }

  //'Store the zMat
  List resZ = Rcpp::List::create(Rcpp::Named("N") = nMat,
                                 Rcpp::Named("R") = rMat,
                                 Rcpp::Named("M") = mMat);
  return resZ;
}

NumericMatrix generateThetaPresenceCova(List zList,double alpha0, double alpha1) {
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


NumericMatrix generateVPresenceCova(List zList,int nLocations, double gamma) {
  //'Getting the M matrix
  NumericMatrix mMat = zList[2];
  //'Total number of communities
  int n_community = mMat.ncol();
  //'Initialize the Phi matrix
  NumericMatrix vMat(nLocations,n_community);
  //'Foreach Specie
  for(int l=0;l<nLocations;l++){
    //'For each community:
    NumericVector nGreater = invertedCumsumPresenceCova(mMat(l,_));
    for(int c=0;c<n_community;c++){
      //'nLC is the number of species in plot l that come from community c
      double nLC = mMat(l,c);
      if(c<n_community-1){
        //'Generate stick-breaking probabilities
        vMat(l,c)=R::rbeta(1.0+nLC,gamma+nGreater(c+1));
      }
      else{
        //'Ensure that the last community has 1
        vMat(l,c)=1.0;
      }
    }
  }
  return vMat;
}

double ll_priorFunctionPresenceCova(NumericMatrix matDATA,int nLocations, int nSpecies,int n_community, NumericMatrix vMat,NumericMatrix Theta, NumericMatrix Phi, double alpha0,double alpha1, double gamma, bool ll_prior=true) {
  //'Initialize the logLikelihoodVec
  double logLikelihood=0;
  //'Calculate the Loglikelihood and Prior
  if(ll_prior){
    //'Initialize the V_{cl} and Theta_{sc} prior
    double priorV=0.0;
    double priorTheta=0.0;
    //'For each location
    for(int l=0;l<nLocations;l++){
      //'Initiate the Theta Counter
      int thetaGibbsCount=0;
      //'Compute the prior for V_{cl}
      for(int c=0;c<n_community;c++){
        if(vMat(l,c)<1)priorV=priorV+R::dbeta(vMat(l,c),1,gamma,1);
      }

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
    logLikelihood=logLikelihood+priorTheta+priorV;
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

//' @name GibbsSamplingPresenceCova
//' @title Gibbs Sampling for LDA PresenceCova and Absence
//' @description Compute the Gibbs Sampling for LDA PresenceCova and Absence
//' @param DATA - DataFrame with PresenceCova and Absecence (Zeros and Ones)
//' @param int n_community - Number of communities
//' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
//' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
//' @param gamma - Hyperparameter  Beta(1,gamma)
//' @param n_gibbs - Total number of Gibbs Samples
//' @param ll_prior - Likelihood compute with Priors ?
//' @param bool display_progress=true - Should I Show the progressBar ?
//' @return List - With Theta(n_gibbs,n_community*nSpecies), Phi(n_gibbs,nLocations*n_community) and logLikelihood
// [[Rcpp::export]]
List lda_bernoulli_cov(DataFrame data,DataFrame design,  int n_community, double alpha0, double alpha1, double gamma, int n_gibbs, bool ll_prior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix matDATA = internal::convert_using_rfunction(data, "as.matrix");

  //'Convert to matrix
  NumericMatrix matDesign = internal::convert_using_rfunction(design, "as.matrix");


  //'Total number of locations
  int nLocations = matDATA.nrow();

  //'Total number of species
  int nSpecies = matDATA.ncol();

  //'Intialize Theta
  NumericVector hyperTheta(n_community);
  hyperTheta.fill(1);
  NumericMatrix Theta=rdirichletPresenceCova(nSpecies,hyperTheta);

  //'Intialize vMat
  NumericVector hyperV(n_community);
  hyperV.fill(1);
  NumericMatrix vMat=rdirichletPresenceCova(nLocations,hyperV);

  //'Initialize the ThetaGibbs
  NumericMatrix ThetaGibbs(n_gibbs,n_community*nSpecies);

  //'Initialize the PhiGibbs
  NumericMatrix PhiGibbs(n_gibbs,nLocations*n_community);

  //'Initialize the BetaGibbs
  NumericMatrix BetaGibbs(n_gibbs,n_community*matDesign.ncol());

  //'Initialize the logLikelihood vector
  NumericVector logLikelihoodVec(n_gibbs);

  //Check if gamma existis
  bool bgamma = false;
  if(std::isnan(gamma)){
    bgamma=true;
  }

  //'Intialize the progressbar
  Progress p(n_gibbs, display_progress);
  for (int g = 0; g < n_gibbs; ++g) {
    //'Verify if everything is ok
    if (Progress::check_abort() )
      Rcpp::stop("Operation cancelled by interrupt.");
    //'Initialize the Phi matrix
    NumericMatrix PhiMat(nLocations, n_community);

    //'Generate zList
    List zList  = generateZPresenceCova(matDATA, Theta, vMat, PhiMat);

    //'Generate Theta
    Theta = generateThetaPresenceCova(zList,alpha0,alpha1);

    //Generate gamma MH
    if(bgamma){
      gamma = gammaMHPresence(vMat, gamma, 0.5);
    }

    //'Generate vMat
    vMat = generateVPresenceCova(zList,nLocations, gamma);

    //Beta coefficientes
    NumericMatrix betasCoef(n_community,matDesign.ncol());
    for(int c=0;c<n_community;c++){

      //Create the logistic constants
      NumericVector logistic = Rcpp::log(PhiMat(_,c))/Rcpp::log(1.0-PhiMat(_,c));


      int count =0;
      //Count number of noniformative
      for(int e=0;e<logistic.size();e++){
        if(!std::isinf(logistic(e))){
          count=count+1;
        }
      }

      //Create only informative data
      arma::mat Xmat(count,matDesign.ncol());
      arma::vec bVec(count);
      count =0;
      for(int e=0;e<logistic.size();e++){
        if(!std::isinf(logistic(e))){
          //Convert the row
          for(int cc=0;cc<matDesign.ncol();cc++){
            Xmat(count,cc) = matDesign(e,cc);
          }
          bVec(count) = logistic(e);
          count=count+1;
        }
      }

      arma::vec betasVec = arma::solve(Xmat,bVec);

      //Store the solution
      for(int e=0;e< betasVec.n_elem;e++){
        betasCoef(c,e) = betasVec(e);
      }
    }

    //'Create the final Theta (n_gibbs,n_community*nSpecies) and final Phi (PhiGibbs)
    //'Create the final BetaGibbs (n_gibbs,n_community*n_covar )
    updateThetaAndPhiPresenceCova(ThetaGibbs, Theta, PhiGibbs, PhiMat,BetaGibbs, betasCoef, g);

    //'Initialize the logLikelihood
    double logLikelihood=ll_priorFunctionPresenceCova(matDATA, nLocations,
                                                       nSpecies,n_community,
                                                       vMat, Theta, PhiMat,
                                                       alpha0, alpha1, gamma,
                                                       ll_prior);
    //'Store the logLikelihood
    logLikelihoodVec(g)=logLikelihood;

    //'Increment the progress bar
    p.increment();

  }

  //'Store the results - Order change to agree with the other functions
  List resTemp = Rcpp::List::create(Rcpp::Named("Theta") = PhiGibbs,
                                    Rcpp::Named("Phi")  = ThetaGibbs,
                                    Rcpp::Named("Beta")  = BetaGibbs,
                                    Rcpp::Named("logLikelihood")  =logLikelihoodVec);

  return resTemp;
}


