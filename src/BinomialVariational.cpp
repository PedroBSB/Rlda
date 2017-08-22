#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "progress.hpp"
#include <iostream>
#include <ctime>
#include <fstream>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;



/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/


arma::mat digammaMat(arma::mat x){
  arma::mat xdi(x.n_rows,x.n_cols);
  for(int r=0;r<x.n_rows;r++){
    for(int c=0;c<x.n_cols;c++){
      xdi(r,c)=R::digamma(x(r,c));
    }
  }
  return(xdi);
}

arma::colvec digammaVec(arma::colvec x){
  arma::colvec xdi(x.n_rows);
  for(int r=0;r<x.n_rows;r++){
    xdi(r)=R::digamma(x(r));
  }
  return(xdi);
}


/***************************************************************************************************************************/
/*********************************            GIBBS SAMPLING FUNCTIONS           *******************************************/
/***************************************************************************************************************************/

/*
 get.m1m0=function(nloc,ncommun,nspp,a,b,c,d){
#get stb prior piece
 stb=digamma(a)-digamma(a+b)
stb[,ncommun]=0 #because V_{lK}=1

res=rep(0,nloc)
for (i in 2:ncommun){
res=res+digamma(b[,i-1])-digamma(a[,i-1]+b[,i-1])
stb[,i]=stb[,i]+res
}

#get data part
pdat1=digamma(c)-digamma(c+d)
pdat0=digamma(d)-digamma(c+d)

m0=m1=array(NA,dim=c(nloc,nspp,ncommun),
            dimnames=list(paste('loc',1:nloc,sep=''),
                          paste('spp',1:nspp,sep=''),
paste('comm',1:ncommun,sep='')))
for (i in 1:nloc){
for (j in 1:nspp){
m1.tmp=exp(pdat1[,j]+stb[i,])
m1[i,j,]=m1.tmp/sum(m1.tmp)

m0.tmp=exp(pdat0[,j]+stb[i,])
m0[i,j,]=m0.tmp/sum(m0.tmp)
}
}
list(m1=m1,m0=m0)
}
*/


// [[Rcpp::export]]
Rcpp::List generateM1M0matrix(int nLocations,int n_community,int n_species, arma::mat matA, arma::mat matB, arma::mat matC, arma::mat matD) {
  //' Calculate digamma
  arma::mat matAdi = digammaMat(matA);
  arma::mat matABdi = digammaMat(matA+matB);
  //' Calculate the difference
  arma::mat stb=matAdi-matABdi;
  //' Normalize last column
  stb.col(n_community-1).fill(0.0);
  arma::colvec res(stb.n_rows);
  res.fill(0.0);
    for (int c=0;c<stb.n_cols-1;c++){
      arma::colvec sumAB = matA.col(c)+matB.col(c);
      res = res + digammaVec(matB.col(c)) - digammaVec(sumAB);
      stb.col(c+1)=stb.col(c+1)+res;
    }

  //' Data hyperparameters
  arma::mat pdat1=digammaMat(matC)-digammaMat(matC+matD);
  arma::mat pdat0=digammaMat(matD)-digammaMat(matC+matD);
  //' Initialize M1 and M0
  arma::cube matM0(nLocations,n_species,n_community);
  arma::cube matM1(nLocations,n_species, n_community);
  for (int i=0;i<nLocations;i++){
    arma::vec stbVec =  arma::conv_to< arma::vec >::from(stb.row(i));
    for (int j=0;j<n_species;j++){
      arma::vec pdat1Vec =  arma::conv_to< arma::vec >::from(pdat1.col(j));
      arma::vec pdat0Vec =  arma::conv_to< arma::vec >::from(pdat0.col(j));
      //'M1
      arma::vec m1vec = arma::exp(pdat1Vec+stbVec);
      double sum1 = arma::sum(m1vec);
      //'M0
      arma::vec m0vec = arma::exp(pdat0Vec+stbVec);
      double sum0 = arma::sum(arma::exp(pdat0Vec+stbVec));
      for(int k=0;k<n_community;k++){
        matM1(i,j,k)=m1vec(k)/sum1;
        matM0(i,j,k)=m0vec(k)/sum0;
      }
    }
  }

  //'Store the results
  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("M1") = matM1,
                                    Rcpp::Named("M0")  = matM0);

  return(resTemp);
}


arma::mat generateThetaBinomialVariational(int nLocations,int n_community, arma::mat matA, arma::mat matB) {
  //'Intialize the Theta matrix
  arma::mat Theta(nLocations,n_community);
  //'Create the mean
  arma::mat meanMatDenominator=matA+matB;
  arma::mat meanMat=matA/meanMatDenominator;
  //'Fix in one the last community
  meanMat.col(n_community-1).fill(1.0);
  //' Define the first two communities
  Theta.col(0)=meanMat.col(0);
  Theta.col(1)=meanMat.col(1)%(1.0-meanMat.col(0));
  for(int c=2;c<n_community;c++){
    arma::mat tmpMat = 1.0-meanMat.cols(0,c-1);
    arma::colvec prod = arma::prod(tmpMat, 1);
    Theta.col(c)=meanMat.col(c)%prod;
  }
  return(Theta);
}

