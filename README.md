## Rlda package

[![Travis-CI Build Status](https://travis-ci.org/PedroBSB/Rlda.svg?branch=master)](https://travis-ci.org/PedroBSB/Rlda)

This project resulted in a new R package called Rlda for which these types of data can be used to mixed-membership clustering through a Bayesian paradigm. https://CRAN.R-project.org/package=Rlda 

With the functions and methods presented here it is possible to work with mixed-membership clustering based on different types of data (i.e., Multinomial, Bernoulli, and Binomial entries). This includes most of applications that uses categorical data as the main information.

In a program language perspective, this packaged used the Rcpp library and C++ code as a support to speed up the gibbs samples, which resulted in an efficient way to obtain the posterior samples for those three types of distributions. 
