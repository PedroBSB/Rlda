#include <iostream>
#include <ctime>
#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "progress.hpp"
// [[Rcpp::depends(RcppProgress)]]
using namespace arma;


/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

NumericMatrix convertSBtoNormal(NumericMatrix vmat,
                                int ncol, int nrow,
                                NumericVector prod) {
  NumericMatrix res(nrow,ncol);

  for(int j=0; j<ncol;j++){
    res(_,j)=vmat(_,j)*prod;
    prod=prod*(1-vmat(_,j));
  }

  return (res);
}

double fixMHRemote(double lo, double hi,double old1,double new1,double jump){
  double jold=R::pnorm5(hi,old1,jump,1,0)-R::pnorm5(lo,old1,jump,1,0);
  double jnew=R::pnorm5(hi,new1,jump,1,0)-R::pnorm5(lo,new1,jump,1,0);
  return(std::log(jold)-std::log(jnew));
}

double tnormRemote(double lo, double hi,double mu, double sig){
  double q1 = R::pnorm5(lo,mu,sig,1,0);
  double q2 = R::pnorm5(hi,mu,sig,1,0);
  double z = R::runif(q1,q2);
  z = R::qnorm5(z,mu,sig,1,0);
  return(z);
}

arma::mat generateOmegaRemote(arma::mat thetaMat, arma::mat Omega, List jumpList, arma::mat remoteMat, int maxBand, double a, double b ){
  //Get the Jump Matrix for Omega
  arma::mat jumpOmega = jumpList[0];
  //Get the Jump Matrix for V
  arma::mat jumpV = jumpList[1];
  //Get the Jump Matrix for Omega
  arma::mat jumpX = jumpList[2];
  //Get the number of Communities
  int n_community = jumpOmega.n_rows;
  //Get the number of Bands
  int n_bands = jumpOmega.n_cols;
  //Get the number of locations
  int n_locations = thetaMat.n_rows;

  //Create the Omega Matrix
  arma::mat omegaMat(n_community,n_bands);
  arma::mat omegaAdjust(n_community,n_bands);
  //Prior Matrix with Beta Distribution
  arma::mat priorOld(n_community,n_bands);
  arma::mat priorNew(n_community,n_bands);


  //Truncated Normal
  for(int c=0;c<omegaMat.n_rows;c++){
    for(int b=0;b<omegaMat.n_cols;b++){
      //Truncated Normal
      omegaMat(c,b) = tnormRemote(0.0, 1.0,Omega(c,b), jumpOmega(c,b));
      //Adjusting the Metropolis-Hasting
      omegaAdjust(c,b) = fixMHRemote(0.0,1.0,Omega(c,b),omegaMat(c,b),jumpOmega(c,b));
      //Calculate the prior
      priorOld(c,b) = R::dbeta(Omega(c,b),a,b,true);
      priorNew(c,b) = R::dbeta(omegaMat(c,b),a,b,true);
    }
  }

  //Create the Probability Matrix
  arma::mat probMatOld = thetaMat*Omega;
  arma::mat probMatNew = thetaMat*omegaMat;

  //Calculate the logLikelihood
  arma::vec llkOld(n_bands);
  arma::vec llkNew(n_bands);


  //For each band
  for(int b=0;b<n_bands;b++){
    //For each location
    for(int l=0;l<n_locations;l++){
      llkOld(b)=llkOld(b)+R::dbinom(remoteMat(l,b),maxBand,probMatOld(l,b),true);
      llkNew(b)=llkNew(b)+R::dbinom(remoteMat(l,b),maxBand,probMatNew(l,b),true);
    }
  }

  //For each community acept ou reject
  for(int c=0;c<n_community;c++){
    //Vector P0
    arma::vec p0 = llkOld+priorOld.row(c);
    //Vector P1
    arma::vec p1 = llkNew+priorNew.row(c)+omegaAdjust.row(c);
    //For each band
    for(int b=0;b<n_bands;b++){
      //Acceptance
      double a = std::exp(p1(b)-p0(b));
      double z = R::unif_rand();
      if(z<a){
        Omega(c,b) = omegaMat(c,b);
      }
    }
  }
  return(Omega);
}











