#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <ctime>
#include <fstream>
#include <cmath>
#include "progress.hpp"
#include "rgens.h"
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;


/***************************************************************************************************************************/
/*********************************                      HEADER            **************************************************/
/***************************************************************************************************************************/

arma::mat rmvnorm(unsigned int n, const arma::vec& mu, const arma::mat& S);
arma::mat rwishart(unsigned int df, const arma::mat& S);
arma::mat riwishart(unsigned int df, const arma::mat& S);

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

void convertRcpptoARMA(NumericMatrix matIn, arma::mat& matOut){
  for(int r=0;r<matIn.nrow();r++){
    for(int c=0;c<matIn.ncol();c++){
      matOut(r,c)=matIn(r,c);
    }
  }
}

void convertRcpptoARMA(NumericVector vetIn, arma::vec& vetOut){
  for(int r=0;r<vetIn.size();r++){
    vetOut(r)=vetIn(r);
  }
}

void convertARMAtoRcpp(arma::mat matIn, NumericMatrix& matOut){
  for(int r=0;r<matIn.n_rows;r++){
    for(int c=0;c<matIn.n_cols;c++){
      matOut(r,c)=matIn(r,c);
    }
  }
}

void convertARMAtoRcpp(arma::vec vetIn, NumericVector& vetOut){
  for(int r=0;r<vetIn.n_elem ;r++){
    vetOut(r)=vetIn(r);
  }
}

int whichLessAbundanceCovariate(double value, NumericVector vector) {
  int res = -1;
  for (int i = 0; i < vector.size(); i++) {
    if (value < vector(i)) {
      res = i;
      break;
    }
  }
  return res;
}

NumericVector rmultinomialAbundanceCovariate(int size, NumericVector prob) {
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
    int iPos=whichLessAbundanceCovariate(random,table);
    //'Increment the matrix
    res(iPos)=res(iPos)+1;
    //'Increment the counter
    count=count+1;
  }

  return res;
}

NumericVector rdirichletAbundanceCovariate(Rcpp::NumericVector parms) {
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

NumericMatrix rdirichletAbundanceCovariate(int n, Rcpp::NumericVector parms) {
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

NumericVector invertedCumsumAbundanceCovariate(NumericVector n){
  NumericVector table(n.length());
  table(n.length()-1)=n(n.length()-1);
  for(int i=(n.length()-2);i>-1;i--){
    table(i)=table(i+1)+n(i);
  }
  return(table);
}


NumericVector countElementsAbundanceCovariate(List zList,int c,int nSpecies){
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

NumericVector meltAbundanceCovariate(NumericMatrix mat){
  //'Initialize the NumericVector
  NumericVector vec(mat.nrow()*mat.ncol());
  //'Initialize the position
  int pos=0;
  for(int col=0;col<mat.ncol();col++){
    for(int row=0;row<mat.nrow();row++){
      //'meltAbundanceCovariate the matrix
      vec(pos)=mat(row,col);
      //'Increment the position
      pos=pos+1;
    }
  }
  return(vec);
}

void updateThetaAndPhiAbundanceCovariate(NumericMatrix &ThetaGibbs,NumericMatrix Theta,NumericMatrix &PhiGibbs, NumericMatrix Phi,int gibbs){
    //'meltAbundanceCovariate the Theta and Phi matrix
    ThetaGibbs(gibbs,_)=meltAbundanceCovariate(Theta);
    PhiGibbs(gibbs,_) =meltAbundanceCovariate(Phi);
}


NumericMatrix sumarizeCommunitiesAbundanceCovariate(List zList, int n_community){
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


double pdfMultivariateNormal(arma::vec x,arma::vec mu, arma::mat Sigma){
  double k = mu.n_elem;
  double part1 = std::pow(M_2_PI,-k/2.0);
  double part2 = std::pow(arma::det(Sigma),-0.5);
  arma::mat part3 = arma::exp(-0.5*((x-mu).t()*arma::inv(Sigma)*(x-mu)));
  return(part1*part2*part3(0,0));
}

//Calculate the acceptance probability
double acceptanceProbability(arma::mat betaNew, arma::mat betaOld,NumericVector xVec,double n, arma::mat Sigma){
  //Casting
  arma::vec xVecArma(xVec.size());
  convertRcpptoARMA(xVec,xVecArma);
  //Numerator
  arma::mat part1 = xVecArma.t()*betaNew;
  arma::mat part2 = arma::exp(-0.5*((betaNew).t()*arma::inv(Sigma)*(betaNew)));
  double num = std::pow(std::exp(part1(0,0)),n)*part2(0,0);

  //Denominator
  part1 = xVecArma.t()*betaOld;
  part2 = arma::exp(-0.5*((betaOld).t()*arma::inv(Sigma)*(betaOld)));
  double den = std::pow(std::exp(part1(0,0)),n)*part2(0,0);

  return(num/den);
}

/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/


List generateThetaAbundanceCovariate(List zList,NumericMatrix xMat, int nLocations,int n_community, NumericMatrix Sigma) {
  //Do the casting
  arma::mat matSigma(Sigma.nrow(),Sigma.ncol());
  convertRcpptoARMA(Sigma, matSigma);
  //'Initialize the Phi matrix
  NumericMatrix thetaMat(nLocations,n_community);
  //'Initialize the Beta matrix
  List betaList(nLocations);

  //'Create the mMat
  NumericMatrix mMat = sumarizeCommunitiesAbundanceCovariate(zList,n_community);

  //'Foreach Specie
  for(int l=0;l<nLocations;l++){
    //'Initialize the Beta matrix
    NumericMatrix betaMat(n_community,xMat.ncol());
    //Initialize the denominator
    double denom = 0.0;
    for(int c=0;c<n_community;c++){
      //Initialize the beta vector
      arma::vec betaVec(xMat.ncol());
      betaVec.fill(0.0);
      //'nLC is the number of species in plot l that come from community c
      double nLC = mMat(l,c);
      //Do the Metropolis Hatsing
      for(int it=0;it<100;it++){
        //Store the old state
        arma::vec betaVecOld(xMat.ncol());
        betaVecOld = betaVec;

        //Generate a proposal state
        arma::mat rNorm = rmvnorm(1, betaVecOld, matSigma);
        for(int k=0;k<rNorm.n_cols;k++){
          betaVec(k) = rNorm(0,k);
        }

        //Calculate the proposal correction factor
        double c = pdfMultivariateNormal(betaVecOld,betaVec,matSigma)/pdfMultivariateNormal(betaVec,betaVecOld,matSigma);

        //Calculate the acceptance probability
        double alpha = acceptanceProbability(betaVec,betaVecOld,xMat(l,_),nLC , matSigma)*c;

        //Draw a random number
        double u = R::runif(0,1);
        if(u>alpha){
          betaVec=betaVecOld;
        }
      }
      //Casting
      NumericVector betaRes(betaVec.n_elem);
      convertARMAtoRcpp(betaVec,betaRes);
      betaMat(c,_)=betaRes;
      thetaMat(l,c)=exp(sum(xMat(l,_)*betaRes));
    }
    thetaMat(l,_)=thetaMat(l,_)/sum(thetaMat(l,_));
    betaList[l]=betaMat;
  }
  //'Store the results
  List resTemp = Rcpp::List::create(Rcpp::Named("Beta") = betaList,
                                    Rcpp::Named("Theta")  = thetaMat);
  return resTemp;
}



List generateZAbundanceCovariate(NumericMatrix size, List listTheta, NumericMatrix Phi) {
  //Pass the Theta Matrix
  NumericMatrix Theta = listTheta[1];
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
      NumericVector tmp = rmultinomialAbundanceCovariate(iSize, probability);
      //'PrintObjectLine(probability);
      zMat(s,_)=tmp;
    }
    //'Store the zMat
    zList(l)=zMat;
  }
  return zList;
}

NumericMatrix generatePhiAbundanceCovariate(int n_community, List zList, NumericVector beta) {
  //'Initialize the Phi matrix
  NumericMatrix phiMat(n_community,beta.length());
  //'Number of Species
  int nSpecies=beta.length();
  //'For each community count
  for(int c=0;c<n_community;c++){
    //'How many members are from community c and specie s (s=1,...,S)?
    NumericVector nSpeciesVec = countElementsAbundanceCovariate(zList,c,nSpecies);
    //'Create the parameter vector Dirichlet
    NumericVector parms = nSpeciesVec+beta+1;
    //'Generate the c-th phi
    NumericVector res = rdirichletAbundanceCovariate(parms);
    //'Store the results
    phiMat(c,_)=res;
  }
  return phiMat;
}



/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING PROCEDURE                  ************************************/
/***************************************************************************************************************************/


//' @name GibbsSamplingAbundanceCovariate
//' @title Gibbs Sampling for LDA AbundanceCovariate with Stick-Breaking
//' @description Compute the Gibbs Sampling for LDA AbundanceCovariate with Stick-Breaking
//' @param data - dataFrame with AbundanceCovariate
//' @param int n_community - Number of communities
//' @param beta - NumericVector for beta (Sx1)
//' @param n_gibbs - Total number of Gibbs Samples
//' @param ll_prior - Likelihood compute with Priors ?
//' @param bool display_progress=true - Should I Show the progressBar ?
//' @return List - With Theta(n_gibbs,nLocations*n_community), Phi(n_gibbs,n_community*nSpecies) and logLikelihood
// [[Rcpp::export]]
List lda_covariate(DataFrame data, DataFrame design, int n_community,NumericVector beta, int n_gibbs, bool ll_prior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix matdata = internal::convert_using_rfunction(data, "as.matrix");

  //'Convert to matrix
  NumericMatrix matDesign = internal::convert_using_rfunction(design, "as.matrix");

  //'Total number of locations
  int nLocations = matdata.nrow();

  //'Total number of species
  int nSpecies = matdata.ncol();

  //Number of covariates
  int nCovariates = matDesign.ncol();

  //'Intialize Theta
  NumericVector hyperTheta(n_community);
  hyperTheta.fill(1);
  NumericMatrix Theta = rdirichletAbundanceCovariate(nSpecies,hyperTheta);

  //'Initialize the Beta matrix
  NumericMatrix betaMat(n_community,matDesign.ncol());
  betaMat.fill(0.0);

  //'Initialize the Beta matrix
  List betaList(nLocations);

  //Initialize the betas
  for(int l=0;l<nLocations;l++){
    betaList[l]=betaMat;
  }

  //Initialize Sigma
  NumericMatrix Sigma(matDesign.ncol(),matDesign.ncol());
  Sigma.fill(0.0);
  Sigma.fill_diag(1.0);

  //'Store the results
  List listTheta = Rcpp::List::create(Rcpp::Named("Beta") = betaList,
                                    Rcpp::Named("Theta")  = Theta);

  //'Initialize Phi
  NumericVector hyperPhi(nSpecies);
  hyperPhi.fill(1);
  NumericMatrix Phi=rdirichletAbundanceCovariate(n_community,hyperPhi);

  //'Initialize the logLikelihood vector
  NumericVector logLikelihoodVec(n_gibbs);

  //List with betas, Theta and Phi
  List res(n_gibbs);

  //'Intialize the progressbar
  Progress p(n_gibbs, display_progress);
  for (int g = 0; g < n_gibbs; ++g) {
    //'Verify if everything is ok
    if (Progress::check_abort() )
      Rcpp::stop("Operation cancelled by interrupt.");

    //'Generate zList
    List zList  = generateZAbundanceCovariate(matdata, listTheta, Phi);

    //'Generate Theta
    listTheta =  generateThetaAbundanceCovariate(zList, matDesign, nLocations, n_community, Sigma);

    //'Generate Phi
    Phi = generatePhiAbundanceCovariate(n_community, zList, beta);

    //'Initialize the logLikelihood
    double logLikelihood=0.0;
    //    double logLikelihood=ll_priorFunctionAbundanceCovariate(zList, g, nSpecies, n_community,
    //                                                       vMat, Theta, Phi, gamma, ll_prior);

    //'Store the results
    List resGibbs = Rcpp::List::create(Rcpp::Named("Theta") = listTheta,
                                      Rcpp::Named("Phi")  = Phi,
                                      Rcpp::Named("logLikelihood")  = logLikelihood);


    res[g] = resGibbs;
    //'Store the logLikelihood
//    logLikelihoodVec(g)=logLikelihood;

    //'Increment the progress bar
    p.increment();

  }

  return res;
}

