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

arma::mat convertSBtoNormal(arma::mat vmat,
                                int ncol, int nrow,
                                arma::colvec prod) {
  arma::mat res(nrow,ncol);

  for(int j=0; j<ncol;j++){
    res.col(j)=vmat.col(j)%prod;
    prod=prod%(1.0-vmat.col(j));
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

arma::vec rowSums(const arma:: mat & X){
  int nRows = X.n_rows;
  arma::vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = sum(X.row(i));
  }
  return(out);
}

/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/
// [[Rcpp::export]]
arma::mat generateOmegaRemote(arma::mat thetaMat, arma::mat Omega, List jumpList, arma::mat remoteMat, int maxBand, double a, double b ){
  //Get the Jump Matrix for Omega
  arma::mat jumpOmega = jumpList[0];
  //Get the Jump Matrix for V
  arma::mat jumpV = jumpList[1];
  //Get the Jump Matrix for X
  arma::mat jumpX = jumpList[2];
  //Get the number of Communities
  int n_community = Omega.n_rows;
  //Get the number of Bands
  int n_bands = Omega.n_cols;
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
  arma::rowvec llkOld(n_bands);
  arma::rowvec llkNew(n_bands);


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
    arma::rowvec tmp = priorOld.row(c);
    arma::rowvec p0 = llkOld + tmp;

    //Vector P1
    tmp = priorNew.row(c)+omegaAdjust.row(c);
    arma::rowvec p1 = llkNew + tmp;

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

// [[Rcpp::export]]
arma::mat generatePhiRemote(arma::mat Theta, arma::mat matX,arma::mat forestMat, List jumpList, arma::vec bPhi, double aPhi){
  //Get the Jump Matrix for Omega
  arma::mat jumpOmega = jumpList[0];
  //Get the Jump Matrix for V
  arma::mat jumpV = jumpList[1];
  //Get the Jump Matrix for X
  arma::mat jumpX = jumpList[2];
  //Get the total number of Species
  int n_species = matX.n_cols;
  //Get the total number of communities
  int n_community = matX.n_rows;
  //Create the new matrix
  arma::mat newXMat(matX.n_rows,matX.n_cols);
  //Fill with ones
  newXMat.fill(1.0);
  //Create the new adjuested matrix
  arma::mat xAdjust(matX.n_rows,matX.n_cols);
  //Create the old prior matrix
  arma::mat priorOld(matX.n_rows,matX.n_cols);
  //Create the new prior matrix
  arma::mat priorNew(matX.n_rows,matX.n_cols);
  //Create the prod vector
  arma::vec prod(n_community);
  for(int c=0;c<n_community;c++)prod(c)=c;

  //For each number of species
  for(int s=0;s<n_species;s++){
    //For each community
    for(int c=0;c<n_community;c++){

      //Generate new X
      newXMat(c,s) = tnormRemote(0.0, 1.0,matX(c,s), jumpX(c,s));

      //Adjusting the Metropolis-Hasting
      xAdjust(c,s) = fixMHRemote(0.0,1.0,matX(c,s),newXMat(c,s),jumpX(c,s));

      //Old Phi
      arma::mat phiOld=convertSBtoNormal(matX,n_species,n_community,prod);

      //New Phi
      arma::mat phiNew=convertSBtoNormal(newXMat,n_species,n_community,prod);

      //Calculate the old probability
      arma::mat pOld = Theta*phiOld;

      //Calculate the new probability
      arma::mat pNew = Theta*phiNew;

      //Calculate the old scalar probability
      double p1Old = arma::sum(arma::sum(forestMat%arma::log(pOld)));

      //Calculate the old scalar probability
      double p1New = arma::sum(arma::sum(forestMat%log(pNew)));

      //Simulate old prior
      priorOld(c,s)= R::dbeta(matX(c,s),aPhi,bPhi(s),true);

      //Simulate new prior
      priorNew(c,s)= R::dbeta(newXMat(c,s),aPhi,bPhi(s),true);

      if(s<n_species-1){
        //Probability p0
        double p0 = p1Old+priorOld(c,s);

        //Probability p1
        double p1 = p1New+priorNew(c,s)+xAdjust(c,s);

        //Acceptance
        double a = std::exp(p1-p0);
        double z = R::unif_rand();
        if(z<a){
          matX(c,s) = newXMat(c,s);
        }
      }
    }
  }

  //Recreate the Phi Matrix
  arma::mat phiMat = convertSBtoNormal(matX,n_species,n_community,prod);

  return(phiMat);
}

// [[Rcpp::export]]
arma::mat generateThetaRemote(arma::mat &vMatrix, arma::mat Omega,arma::mat Phi, arma::mat forestMat, List jumpList, double gamma){
  //Get the Jump Matrix for Omega
  arma::mat jumpOmega = jumpList[0];
  //Get the Jump Matrix for V
  arma::mat jumpV = jumpList[1];
  //Get the Jump Matrix for Omega
  arma::mat jumpX = jumpList[2];
  //Get the number of communities
  int n_community = jumpV.n_cols;
  //Get the total number of locations
  int n_locations = jumpV.n_rows;
  //Initialize the new V matrix
  arma::mat newVMat(n_locations,n_community);
  newVMat.fill(1.0);
  arma::mat vAdjustedMat(n_locations,n_community);

    //Create the old prior matrix
  arma::mat priorOld(vMatrix.n_rows,vMatrix.n_cols);
  //Create the new prior matrix
  arma::mat priorNew(vMatrix.n_rows,vMatrix.n_cols);

  //Create the prod vector
  arma::vec prod(n_locations);
  for(int l=0;l<n_locations;l++)prod(l)=l;

  for(int l=0;l<n_locations;l++){
    for(int c=0;c<n_community;c++){
      //Generate new X
      if(c<n_community-1) newVMat(l,c) = tnormRemote(0.0, 1.0,vMatrix(l,c), jumpV(l,c));
      //Adjusting the Metropolis-Hasting
      vAdjustedMat(l,c) = fixMHRemote(0.0,1.0,vMatrix(l,c),newVMat(l,c),jumpV(l,c));

      //Simulate old prior
      priorOld(l,c)= R::dbeta(vMatrix(l,c),1.0,gamma,true);
      //Simulate new prior
      priorNew(l,c)= R::dbeta(newVMat(l,c),1.0,gamma,true);

      //Old Theta
      arma::mat thetaOld=convertSBtoNormal(vMatrix,n_community,n_locations,prod);
      //New Theta
      arma::mat thetaNew=convertSBtoNormal(newVMat,n_community,n_locations,prod);

      //Calculate the old probability
      arma::mat pOld = thetaOld*Omega;
      //Calculate the new probability
      arma::mat pNew = thetaNew*Omega;

      //Calculate the vector old probability
      arma::vec pOld2 = rowSums(forestMat%arma::log(thetaOld*Phi));
      //Calculate the vector new probability
      arma::vec pNew2 = rowSums(forestMat%arma::log(thetaNew*Phi));

      if(c<n_community-1){
        //Probability p0
        double p0 = pOld2(l)+priorOld(l,c);
        //Probability p1
        double p1 = pNew2(l)+priorNew(l,c)+vAdjustedMat(l,c);
        //Acceptance
        double a = std::exp(p1-p0);
        double z = R::unif_rand();
        if(z<a){
          vMatrix(l,c) = newVMat(l,c);
        }
      }
    }
  }
  arma::mat Theta = convertSBtoNormal(vMatrix,n_community,n_locations,prod);
  return(Theta);
}






