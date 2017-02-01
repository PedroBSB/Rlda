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

update.omega=function(param,jump){
  omega.orig=omega.old=param$omega
  tmp=tnorm(nbands*ncommun,lo=0,hi=1,mu=omega.old,sig=jump)
  novos=matrix(tmp,ncommun,nbands)

  tmp=fix.MH(lo=0,hi=1,omega.old,novos,jump)
  ajuste=matrix(tmp,ncommun,nbands)

  prior.old=matrix(dbeta(omega.old,a.omega,b.omega,log=T),ncommun,nbands)
  prior.new=matrix(dbeta(novos    ,a.omega,b.omega,log=T),ncommun,nbands)

  for (i in 1:ncommun){
    omega.new=omega.old
    omega.new[i,]=novos[i,]

    prob.old=param$theta%*%omega.old
    llk.old=colSums(dbinom(remote,size=ndig.values,prob=prob.old,log=T))

    prob.new=param$theta%*%omega.new
    llk.new=colSums(dbinom(remote,size=ndig.values,prob=prob.new,log=T))

    k=acceptMH(llk.old+prior.old[i,],
               llk.new+prior.new[i,]+ajuste[i,],
                                           omega.old[i,],omega.new[i,],F)
    omega.old[i,]=k$x
  }
  list(omega=omega.old,accept=omega.old!=omega.orig)
}











