#include <RcppArmadillo.h>
#include <iostream>
#include <ctime>
#include <fstream>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "progress.hpp"
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/



int whichLessCovariate(double value, NumericVector vector) {
  int res = -1;
  for (int i = 0; i < vector.size(); i++) {
    if (value < vector(i)) {
      res = i;
      break;
    }
  }
  return res;
}

NumericVector rmultinomialCovariate(int size, NumericVector prob) {
  //'Initialize the NumericMatrix result
  NumericVector res(prob.length());
  //'Create the categorical table
  NumericVector table(prob.length());
  //'Populate the categorical table
  table(0)=prob(0);
  for(int i=1;i<prob.length();i++){
    table(i)=table(i-1)+prob(i);
  }
  //'Generate the sample
  int count=0;
  //'Generate until get the sample size
  while(count<size){
    //'Draw a uniform
    double random = R::runif(0,1);
    //'Which category was draw ?
    int iPos=whichLessCovariate(random,table);
    //'Increment the matrix
    res(iPos)=res(iPos)+1;
    //'Increment the counter
    count=count+1;
  }

  return res;
}

NumericVector rdirichletCovariate(Rcpp::NumericVector parms) {
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

NumericMatrix rdirichletCovariate(int n, Rcpp::NumericVector parms) {
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

NumericVector invertedCumsumCovariate(NumericVector n){
  NumericVector table(n.length());
  table(n.length()-1)=n(n.length()-1);
  for(int i=(n.length()-2);i>-1;i--){
    table(i)=table(i+1)+n(i);
  }
  return(table);
}


NumericVector countElementsCovariate(List zList,int c,int nSpecies){
  //'Size of the list
  int nSize=zList.length();
  //'Initialize the vector with count
  NumericVector vecSpecie(nSpecies);
  for(int i=0;i<nSize;i++){
    NumericMatrix zMat = zList[i];
    vecSpecie = vecSpecie + zMat(_,c);
  }
  return(vecSpecie);
}

NumericVector meltCovariate(NumericMatrix mat){
  //'Initialize the NumericVector
  NumericVector vec(mat.nrow()*mat.ncol());
  //'Initialize the position
  int pos=0;
  for(int col=0;col<mat.ncol();col++){
    for(int row=0;row<mat.nrow();row++){
      //'meltCovariate the matrix
      vec(pos)=mat(row,col);
      //'Increment the position
      pos=pos+1;
    }
  }
  return(vec);
}

void updateThetaAndPhiCovariate(NumericMatrix &ThetaGibbs,NumericMatrix Theta,NumericMatrix &PhiGibbs, NumericMatrix Phi,NumericMatrix &BetaGibbs, NumericMatrix Beta,int gibbs){
  //'meltCovariate the Theta and Phi matrix
  ThetaGibbs(gibbs,_)=meltCovariate(Theta);
  PhiGibbs(gibbs,_) =meltCovariate(Phi);
  BetaGibbs(gibbs,_) =meltCovariate(Beta);
}


NumericMatrix sumarizeCommunitiesCovariate(List zList, int n_community){
  //'Total number of locations
  int nLocations = zList.length();
  //'Intialize the mMat
  NumericMatrix mMat(nLocations,n_community);
  //'Foreach location
  for(int l=0;l<nLocations;l++){
    //'Get the zMat
    NumericMatrix zMat=zList[l];
    //'For each community
    for(int c=0;c<n_community;c++){
      //'Sum all elements in this location
      mMat(l,c)=sum(zMat(_,c));
    }
  }
  return(mMat);
}



/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/

NumericMatrix generateThetaCovariate(NumericMatrix vMat) {
  //'Total number of locations
  int nLocations = vMat.nrow();
  //'Total number of communities
  int n_community = vMat.ncol();
  //'Initialize the Theta matrix
  NumericMatrix thetaMat(nLocations,n_community);
  //'Foreach location
  for(int l=0;l<nLocations;l++){
    NumericVector Theta(n_community);
    //'Update the Theta \prod_(k=1)^(c-1)(1-V_kl )
    double prod=1;
    //'For each community
    for(int c=0;c<n_community;c++){
      double vNumber = vMat(l,c);
      if (c == 0) prod=1;
      if (c >  0) prod=prod*(1.0-vMat(l,c-1));
      Theta(c)=vNumber*prod;
    }
    //'Store each row
    thetaMat(l,_)=Theta;
  }
  return thetaMat;
}

List generateZCovariate(NumericMatrix size, NumericMatrix Theta, NumericMatrix Phi) {
  //'Initialize the Z List
  List zList(Theta.nrow());
  //'Number of locations
  int nLocations=Theta.nrow();
  //'Number of Species
  int nSpecies=Phi.ncol();
  //'Number of Communities
  int n_community=Phi.nrow();
  //'For each location sample from a Multinomial
  for(int l=0;l<nLocations;l++){
    //'Initialize the zMat:
    NumericMatrix zMat(nSpecies,n_community);
    //'For each Specie
    for(int s=0;s<nSpecies;s++){
      //'Calculate the probability vector
      double sumVec=0;
      NumericVector probability = Theta(l,_)*Phi(_,s);
      //'Find the sum
      sumVec=sum(probability);
      //'Normalize the probability
      probability=probability/sumVec;
      //'Number of elements in location l and specie s
      int iSize=size(l,s);
      //'Store the results
      NumericVector tmp = rmultinomialCovariate(iSize, probability);
      //'PrintObjectLine(probability);
      zMat(s,_)=tmp;
    }
    //'Store the zMat
    zList(l)=zMat;
  }
  return zList;
}

NumericMatrix generatePhiCovariate(int n_community, List zList, NumericVector beta) {
  //'Initialize the Phi matrix
  NumericMatrix phiMat(n_community,beta.length());
  //'Number of Species
  int nSpecies=beta.length();
  //'For each community count
  for(int c=0;c<n_community;c++){
    //'How many members are from community c and specie s (s=1,...,S)?
    NumericVector nSpeciesVec = countElementsCovariate(zList,c,nSpecies);
    //'Create the parameter vector Dirichlet
    NumericVector parms = nSpeciesVec+beta+1;
    //'Generate the c-th phi
    NumericVector res = rdirichletCovariate(parms);
    //'Store the results
    phiMat(c,_)=res;
  }
  return phiMat;
}

NumericMatrix generateVCovariate(List zList,int nLocations,int n_community, double gamma) {
  //'Initialize the Phi matrix
  NumericMatrix vMat(nLocations,n_community);
  //'Create the mMat
  NumericMatrix mMat = sumarizeCommunitiesCovariate(zList,n_community);
  //'Foreach Specie
  for(int l=0;l<nLocations;l++){
    //'For each community:
    NumericVector nGreater = invertedCumsumCovariate(mMat(l,_));
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


double ll_priorFunctionCovariate(List zList, int g, int nSpecies,int n_community, NumericMatrix vMat,NumericMatrix Theta, NumericMatrix Phi, double gamma, bool ll_prior=true) {
  double logLikelihood=0;
  //'Calculate the Loglikelihood and Prior
  if(ll_prior){
    //'Initialize
    double priorV=0.0;
    double priorPhi=0.0;
    //'Calculate the likelihood based on data
    //'For each location
    for(int l=0;l<zList.length();l++){
      //'Compute the prior for V_{cl}
      for(int c=0;c<n_community;c++){
        if(vMat(l,c)<1)priorV=priorV+R::dbeta(vMat(l,c),1,gamma,1);
      }
      //'Get the zMat
      NumericMatrix zMat = zList[l];
      //'For each specie
      for(int s=0;s<nSpecies;s++){
        //'For each community
        for(int c=0;c<n_community;c++){
          //'Acumulate the Prior Phi
          if(zMat(s,c)>0){
            //'All elements in specie S and Community C
            if(Phi(c,s)>0) logLikelihood=logLikelihood+zMat(s,c)*log(Phi(c,s));
            //'All elements in location L and community C
            if(Theta(l,c)>0)logLikelihood=logLikelihood+sum(zMat(_,c))*log(Theta(l,c));
          }
        }
      }
    }
    logLikelihood=logLikelihood+priorV+priorPhi;
  }
  else{
    //'Calculate the likelihood based on data
    //'For each location
    for(int l=0;l<zList.length();l++){
      //'Get the zMat
      NumericMatrix zMat = zList[l];
      //'For each specie
      for(int s=0;s<nSpecies;s++){
        //'For each community
        for(int c=0;c<n_community;c++){
          if(zMat(s,c)>0){
            //'All elements in specie S and Community C
            if(Phi(c,s)>0) logLikelihood=logLikelihood+zMat(s,c)*log(Phi(c,s));
            //'All elements in location L and community C
            if(Theta(l,c)>0)logLikelihood=logLikelihood+sum(zMat(_,c))*log(Theta(l,c));
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

//' @name GibbsSamplingCovariate
//' @title Gibbs Sampling for LDA Covariate with Stick-Breaking
//' @description Compute the Gibbs Sampling for LDA Covariate with Stick-Breaking
//' @param data - dataFrame with Covariate
//' @param int n_community - Number of communities
//' @param beta - NumericVector for beta (Sx1)
//' @param gamma - Hyperparameter  Beta(1,gamma)
//' @param n_gibbs - Total number of Gibbs Samples
//' @param ll_prior - Likelihood compute with Priors ?
//' @param bool display_progress=true - Should I Show the progressBar ?
//' @return List - With Theta(n_gibbs,nLocations*n_community), Phi(n_gibbs,n_community*nSpecies) and logLikelihood
// [[Rcpp::export]]
List lda_multinomial_cov(DataFrame data, DataFrame design, int n_community, NumericVector beta, double gamma, int n_gibbs, bool ll_prior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix matdata = internal::convert_using_rfunction(data, "as.matrix");

  //'Convert to matrix
  NumericMatrix matDesign = internal::convert_using_rfunction(design, "as.matrix");

  //'Total number of locations
  int nLocations = matdata.nrow();

  //'Total number of species
  int nSpecies = matdata.ncol();

  //'Initialize the ThetaGibbs
  NumericMatrix ThetaGibbs(n_gibbs,nLocations*n_community);

  //'Initialize the PhiGibbs
  NumericMatrix PhiGibbs(n_gibbs,n_community*nSpecies);

  //'Initialize the BetaGibbs
  NumericMatrix BetaGibbs(n_gibbs,n_community*matDesign.ncol());

  //'Intialize Theta
  NumericVector hyperTheta(n_community);
  hyperTheta.fill(1);
  NumericMatrix Theta=rdirichletCovariate(nSpecies,hyperTheta);

  //'Initialize Phi
  NumericVector hyperPhi(nSpecies);
  hyperPhi.fill(1);
  NumericMatrix Phi=rdirichletCovariate(n_community,hyperPhi);

  //'Intialize vMat
  NumericVector hyperV(n_community);
  hyperV.fill(1);
  NumericMatrix vMat=rdirichletCovariate(nLocations,hyperV);

  //'Initialize the logLikelihood vector
  NumericVector logLikelihoodVec(n_gibbs);

  //'Intialize the progressbar
  Progress p(n_gibbs, display_progress);
  for (int g = 0; g < n_gibbs; ++g) {
    //'Verify if everything is ok
    if (Progress::check_abort() )
      Rcpp::stop("Operation cancelled by interrupt.");

    //'Generate Theta
    Theta = generateThetaCovariate(vMat);

    //'Generate zList
    List zList  = generateZCovariate(matdata, Theta, Phi);

    //'Generate Phi
    Phi = generatePhiCovariate(n_community, zList, beta);

    //'Generate vMat
    vMat = generateVCovariate(zList,nLocations,n_community, gamma);

    //Beta coefficientes
    NumericMatrix betasCoef(n_community,matDesign.ncol());
    for(int c=0;c<n_community;c++){

      //Create the logistic constants
      NumericVector logistic = Rcpp::log(Theta(_,c))/Rcpp::log(1.0-Theta(_,c));


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

/*
      //Replace extreme values
      for(int e=0;e<logistic.size();e++){
        if(logistic(e)<0 && std::isinf(logistic(e))){
          logistic(e)=1e-10;
        }
        else if(std::isinf(logistic(e)>0 && logistic(e))){
          logistic(e)=1e+10;
        }
      }
     //Solve the logistic link function
     arma::mat Xmat = as<arma::mat>(wrap(matDesign));
     arma::vec bVec = as<arma::vec>(wrap(logistic));
*/

      arma::vec betasVec = arma::solve(Xmat,bVec);

      //Store the solution
      for(int e=0;e< betasVec.n_elem;e++){
        betasCoef(c,e) = betasVec(e);
      }
    }

    //'Create the final ThetaGibbs (n_gibbs,nLocations*n_community) and final PhiGibbs (n_gibbs,n_community*nSpecies)
    //'Create the final BetaGibbs (n_gibbs,n_community*n_covar )
    updateThetaAndPhiCovariate(ThetaGibbs, Theta, PhiGibbs, Phi,BetaGibbs, betasCoef, g);

    //'Initialize the logLikelihood
    double logLikelihood=ll_priorFunctionCovariate(zList, g, nSpecies, n_community,
                                                   vMat, Theta, Phi, gamma, ll_prior);

    //'Store the logLikelihood
    logLikelihoodVec(g)=logLikelihood;

    //'Increment the progress bar
    p.increment();

  }

  //'Store the results
  List resTemp = Rcpp::List::create(Rcpp::Named("Theta") = ThetaGibbs,
                                    Rcpp::Named("Phi")  = PhiGibbs,
                                    Rcpp::Named("Beta")  = BetaGibbs,
                                    Rcpp::Named("logLikelihood")  = logLikelihoodVec);

  return resTemp;
}
