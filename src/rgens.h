#ifndef RGENS_H
#define RGENS_H
arma::mat rmvnorm(unsigned int n, const arma::vec& mu, const arma::mat& S);

arma::mat rwishart(unsigned int df, const arma::mat& S);

arma::mat riwishart(unsigned int df, const arma::mat& S);
  
#endif