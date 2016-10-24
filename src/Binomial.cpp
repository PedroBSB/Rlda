#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

double matchBinomial(double value,NumericVector vec){
  for(int i=0;i<vec.length();i++){
    if(value==vec(i)){
      return(i);
    }
  }
  return(-1);
}

NumericVector rdirichletBinomial(Rcpp::NumericVector parms) {
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

NumericMatrix rdirichletBinomial(int n, Rcpp::NumericVector parms) {
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


int whichLessDVBinomial(double value, NumericVector prob) {
  int res = -1;

  //'Create the categorical table
  double probcum = 0;

  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value <= probcum) {
      res = i;
      break;
    }
  }
  return res;
}

NumericVector rmultinomialDVBinomial(int size, NumericVector prob) {
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
    int iPos=whichLessDVBinomial(random,prob);
    //'Increment the matrix
    res(iPos)=res(iPos)+1;
    //'Increment the counter
    count=count+1;
  }

  return res;
}

NumericVector invertedCumsumBinomial(NumericVector n){
  NumericVector table(n.length());
  table(n.length()-1)=n(n.length()-1);
  for(int i=(n.length()-2);i>-1;i--){
    table(i)=table(i+1)+n(i);
  }
  return(table);
}

NumericVector meltBinomial(NumericMatrix mat){
  //'Initialize the NumericVector
  NumericVector vec(mat.nrow()*mat.ncol());
  //'Initialize the position
  int pos=0;
  for(int col=0;col<mat.ncol();col++){
    for(int row=0;row<mat.nrow();row++){
      //'meltBinomial the matrix
      vec(pos)=mat(row,col);
      //'Increment the position
      pos=pos+1;
    }
  }
  return(vec);
}

void updateThetaAndPhiBinomial(NumericMatrix &ThetaGibbs,NumericMatrix Theta,NumericMatrix &PhiGibbs, NumericMatrix Phi,int gibbs){
    //'meltBinomial the Theta and Phi matrix
    ThetaGibbs(gibbs,_)=meltBinomial(Theta);
    PhiGibbs(gibbs,_) =meltBinomial(Phi);
}

NumericMatrix aggregateValuesBinomial(NumericMatrix zMat,int nLocations,int nCommunity){
  //'Initialize the matrix
  NumericMatrix sum(nLocations,nCommunity);
  //'Fill the sum matrix
  sum.fill(0.0);
  for(int i=0;i<zMat.nrow();i++){
    //'Get the community
    int iCommunity = zMat(i,3);
    //'Get the location
    int iLocation = zMat(i,1);
    //'Agregate
    sum(iLocation,iCommunity)=sum(iLocation,iCommunity)+zMat(i,4);
  }
  //'Return the results
  return(sum);
}

List aggregateValuesBinomialByReflectanceBinomial(NumericMatrix zMat,int nBands,int nCommunity){
  //'Initialize the results
  List res(2);
  //'Initialize the matrix with Refletance
  NumericMatrix sumRef(nBands,nCommunity);
  //'Fill the sum matrix
  sumRef.fill(0.0);
  //'Initialize the matrix without Refletance
  NumericMatrix sumNonRef(nBands,nCommunity);
  //'Fill the sum matrix
  sumNonRef.fill(0.0);
  for(int i=0;i<zMat.nrow();i++){
    //'Get the community
    int iCommunity = zMat(i,3);
    //'Get the band
    int iBand = zMat(i,0);
    //'Get the reflectance
    int iRefletance = zMat(i,2);
    if(iRefletance==1){
      //'Agregate
      sumRef(iBand,iCommunity)=sumRef(iBand,iCommunity)+zMat(i,4);
    }
    else{
      //'Agregate
      sumNonRef(iBand,iCommunity)=sumNonRef(iBand,iCommunity)+zMat(i,4);
    }
  }
  //'Return the results
  res[0]=sumRef;
  res[1]=sumNonRef;
  return(res);

}

double sumLargestBinomial(NumericMatrix nSum, int iCommunity){
  double sum=0.0;
  //'For each location
  for(int l=0;l<nSum.nrow();l++){
    //'For each community
    for(int c=(iCommunity+1);c<nSum.ncol();c++){
      sum=sum+nSum(l,c);
    }
  }
  return(sum);
}

NumericMatrix mmultBinomial(const NumericMatrix& m1, const NumericMatrix& m2){
  if (m1.ncol() != m2.nrow()) stop ("Incompatible matrix dimensions");
  NumericMatrix out(m1.nrow(),m2.ncol());
  NumericVector rm1, cm2;
  for (size_t i = 0; i < m1.nrow(); ++i) {
    rm1 = m1(i,_);
    for (size_t j = 0; j < m2.ncol(); ++j) {
      cm2 = m2(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);
    }
  }
  return out;
}

/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/

NumericMatrix generateZBinomial(NumericMatrix binomMat,NumericMatrix populMat, NumericMatrix Theta, NumericMatrix Phi) {
  try{
    //'Number of locations
    int nLocations=binomMat.nrow();
    //'Number of Bands
    int nBands=binomMat.ncol();
    //'Create the Z matrix
    NumericMatrix zMat(2*nLocations*nBands,5);
    //'Counter
    int iCount=0;
    //'For each band
    for(int b=0;b<nBands;b++){
      //'For each location
      for(int l=0;l<nLocations;l++){
        //'Two cases: Reflectance and Not Reflectance
        //'Calculate the probability of belonging to each community:
        NumericVector prob=Phi(b,_)*Theta(l,_);
        if(Rcpp::sum(prob)==0){
          prob = Rcpp::rep(1.0/Phi.ncol(),Phi.ncol());
        }
        else{
          prob=prob/Rcpp::sum(prob);
        }
        prob=prob/sum(prob);
        if(binomMat(l,b)==0){
          //'Presence in locatiol l and band b (always will be one draw)
          int iSize=1;
          //'Store the results
          NumericVector tmp = rmultinomialDVBinomial(iSize, prob);
          //PrintObjectLine(prob);

          //'Find the community
          int iCommunity = (int) matchBinomial(1.0,tmp);

          //'Store the results
          zMat(iCount,0)=b;             //'Store the band
          zMat(iCount,1)=l;             //'Store the location
          zMat(iCount,2)=1;             //'Store the reflectance
          zMat(iCount,3)=iCommunity;    //'Store the community
          zMat(iCount,4)=binomMat(l,b); //'Store the size presence
          iCount=iCount+1;
        }
        //'Calculate the probability of belonging to each community:
        prob=(1.0-Phi(b,_))*Theta(l,_);
        if(Rcpp::sum(prob)==0){
          prob = Rcpp::rep(1.0/Phi.ncol(),Phi.ncol());
        }
        else{
          prob=prob/Rcpp::sum(prob);
        }

        if(binomMat(l,b)!=0){
          //'Absence in locatiol l and band b (always will be one draw)
          int iSize=1;
          //'Store the results
          NumericVector tmp = rmultinomialDVBinomial(iSize, prob);

          //'Find the community
          int iCommunity = (int) matchBinomial(1.0,tmp);

          //'Store the results
          zMat(iCount,0)=b;             //'Store the band
          zMat(iCount,1)=l;             //'Store the location
          zMat(iCount,2)=0;             //'Store the not reflectance
          zMat(iCount,3)=iCommunity;    //'Store the community
          zMat(iCount,4)=populMat(l,b)-binomMat(l,b); //'Store the size absence
          iCount=iCount+1;
        }
      }
    }
    return zMat;
  }
  catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

NumericMatrix generateThetaBinomial(NumericMatrix zMat,NumericMatrix& vMat, int nLocations,int nCommunity, double gamma) {
  //'Intialize the Theta matrix
  NumericMatrix Theta(nLocations,nCommunity);
  //'Calculate the number of individuals in each community and locations
  NumericMatrix sumMat = aggregateValuesBinomial(zMat, nLocations, nCommunity);
  //'For each location
  for(int l=0;l<nLocations;l++){
    //'For each community
    for(int c=0;c<nCommunity;c++){
      if(c<(nCommunity-1)){
        //'How many elements belong to community c in location l
        double nLC=sumMat(l,c);
        //'How many elements are larger than community c
        double nLCgreter=sumLargestBinomial(sumMat, c);
        vMat(l,c)=R::beta(nLC+1.0,nLCgreter+gamma);
      }
      else{
        //'All locations for the last community are one
        vMat(l,c)=1;
      }
    }
  }
  //'Create the Theta matrix
  //'Foreach location
  for(int l=0;l<nLocations;l++){
    NumericVector thetaVec(nCommunity);
    //'Update the Theta \prod_(k=1)^(c-1)(1-V_kl )
    double prod=1;
    //'For each community
    for(int c=0;c<nCommunity;c++){
      double vNumber = vMat(l,c);
      if (c == 0) prod=1;
      if (c >  0) prod=prod*(1.0-vMat(l,c-1));
      thetaVec(c)=vNumber*prod;
    }
    //'Store each row
    Theta(l,_)=thetaVec;
  }
  return(Theta);
}


NumericMatrix generatePhiBinomial(NumericMatrix zMat,int nBands,int nCommunity,double alpha0, double alpha1) {
  //'Initialize the Phi matrix
  NumericMatrix Phi(nBands,nCommunity);
  //'Get the sum matrices
  List sumList = aggregateValuesBinomialByReflectanceBinomial(zMat,nBands,nCommunity);
  NumericMatrix sumRef= sumList(0);
  NumericMatrix sumNonRef=sumList(1);
  //'Generate the Phi
  for(int b=0;b<nBands;b++){
    //'For each community
    for(int c=0;c<nCommunity;c++){
      Phi(b,c)=R::rbeta(alpha0+sumRef(b,c),alpha1+sumNonRef(b,c));
    }
  }
  return(Phi);
}


double logLikelihoodAndPriorFunctionBinomial(NumericMatrix matDATA,NumericMatrix matPOP, NumericMatrix vMat, NumericMatrix Theta, NumericMatrix Phi, double alpha0,double alpha1, double gamma, bool logLikelihoodAndPrior=true) {
  //'Initialize the logLikelihoodVec
  double logLikelihood=0;
  //'Total number of locations
  int nLocations=matDATA.nrow();
  //'Number of bands
  int nBands=matDATA.ncol();
  //'Matrix of probabilities
  NumericMatrix tPhi=Rcpp::transpose(Phi);
  NumericMatrix probs=mmultBinomial(Theta,tPhi);
  //'Calculate the Loglikelihood and Prior
  if(logLikelihoodAndPrior){
    //'Initialize the V_{cl} and Theta_{sc} prior
    double priorV=0.0;
    double priorPhi=0.0;
    //'For each location
    for(int l=0;l<nLocations;l++){
      //'For each band
      for(int b=0;b<nBands;b++){
        if(Phi(l,b)>0 && Phi(l,b)<1) priorPhi=R::dbeta(Phi(l,b),alpha0,alpha1,true);
        if(vMat(l,b)>0 && vMat(l,b)<1)  priorV=R::dbeta(vMat(l,b),1,gamma,true);
        logLikelihood=logLikelihood+R::dbinom(matDATA(l,b),matPOP(l,b),probs(l,b),true)+priorV+priorPhi;
      }
    }
  }
  else{
    //'For each location
    for(int l=0;l<nLocations;l++){
      //'For each band
      for(int b=0;b<nBands;b++){
        logLikelihood=logLikelihood+R::dbinom(matDATA(l,b),matPOP(l,b),probs(l,b),true);
      }
    }
  }
  return(logLikelihood);
}


/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING PROCEDURE                  ************************************/
/***************************************************************************************************************************/

//' @name GibbsSamplingBinomial
//' @title Compute the Gibbs Sampling for LDA Binomial
//' @description Compute the Gibbs Sampling for LDA Binomial
//' @param DATA - DataFrame with Presence and Absecence (Binomial)
//' @param POP - DataFrame with Population Size (Binomial)
//' @param int nCommunity - Number of communities
//' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
//' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
//' @param gamma - Hyperparameter  Beta(1,gamma)
//' @param nGibbs - Total number of Gibbs Samples
//' @param logLikelihoodAndPrior - Likelihood compute with Priors ?
//' @param bool display_progress=true - Should I Show the progressBar ?
//' @return List - With Theta(nGibbs,nCommunity*nSpecies), Phi(nGibbs,nLocations*nCommunity) and logLikelihood
// [[Rcpp::export]]
List GibbsSamplingBinomial(DataFrame DATA,DataFrame POP, int nCommunity, double alpha0, double alpha1, double gamma, int nGibbs, bool logLikelihoodAndPrior=true, bool display_progress=true) {

  //'Convert to matrix
  NumericMatrix matDATA = internal::convert_using_rfunction(DATA, "as.matrix");

  //'Convert to matrix
  NumericMatrix matPOP = internal::convert_using_rfunction(POP, "as.matrix");

  //'Total number of locations
  int nLocations = matDATA.nrow();

  //'Total number of bands
  int nBands = matDATA.ncol();

  //'Intialize Theta
  NumericVector hyperTheta(nCommunity);
  hyperTheta.fill(1);
  NumericMatrix Theta=rdirichletBinomial(nLocations,hyperTheta);
  NumericMatrix vMat(nLocations,nCommunity);

  //'Intialize Phi
  NumericVector hyperPhi(nCommunity);
  hyperPhi.fill(1);
  NumericMatrix Phi=rdirichletBinomial(nBands,hyperPhi);

  //'Initialize the ThetaGibbs
  NumericMatrix ThetaGibbs(nGibbs,nLocations*nCommunity);

  //'Initialize the PhiGibbs
  NumericMatrix PhiGibbs(nGibbs,nBands*nCommunity);

  //'Initialize the logLikelihood vector
  NumericVector logLikelihoodVec(nGibbs);

  //'Intialize the progressbar
  Progress p(nGibbs, display_progress);
  for (int g = 0; g < nGibbs; ++g) {
    //'Verify if everything is ok
    if (Progress::check_abort()) return -1.0;
    //'Initialize the zMAt
    NumericMatrix zMat = generateZBinomial(matDATA, matPOP, Theta, Phi);
    //'Generate Theta
    Theta =generateThetaBinomial(zMat,vMat, nLocations, nCommunity, gamma);
    //'Generate Phi
    Phi = generatePhiBinomial(zMat, nBands, nCommunity, alpha0, alpha1);
    //'Create the final Theta (nGibbs,nLocations*nCommunity) and final Phi (nGibbs,nBands*nCommunity)
    updateThetaAndPhiBinomial(ThetaGibbs, Theta, PhiGibbs, Phi, g);
    //'Initialize the logLikelihood
    double logLikelihood=logLikelihoodAndPriorFunctionBinomial(matDATA,matPOP,
                                                       vMat, Theta, Phi,
                                                       alpha0, alpha1, gamma,
                                                       logLikelihoodAndPrior);
    //'Store the logLikelihood
    logLikelihoodVec(g)=logLikelihood;

    //'Increment the progress bar
    p.increment();

  }

  //'Store the results
  List resTemp = Rcpp::List::create(Rcpp::Named("Theta") = ThetaGibbs,
                                    Rcpp::Named("Phi")  = PhiGibbs,
                                    Rcpp::Named("logLikelihood")  =logLikelihoodVec);

  return resTemp;
}


