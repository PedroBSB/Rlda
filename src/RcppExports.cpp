// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// lda_multinomial
List lda_multinomial(DataFrame data, int n_community, NumericVector beta, double gamma, int n_gibbs, bool ll_prior, bool display_progress);
RcppExport SEXP _Rlda_lda_multinomial(SEXP dataSEXP, SEXP n_communitySEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP n_gibbsSEXP, SEXP ll_priorSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n_community(n_communitySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type n_gibbs(n_gibbsSEXP);
    Rcpp::traits::input_parameter< bool >::type ll_prior(ll_priorSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_multinomial(data, n_community, beta, gamma, n_gibbs, ll_prior, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// getks
Rcpp::List getks(IntegerMatrix z, int ncommun, IntegerMatrix dat);
RcppExport SEXP _Rlda_getks(SEXP zSEXP, SEXP ncommunSEXP, SEXP datSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type ncommun(ncommunSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type dat(datSEXP);
    rcpp_result_gen = Rcpp::wrap(getks(z, ncommun, dat));
    return rcpp_result_gen;
END_RCPP
}
// getlk
IntegerMatrix getlk(IntegerMatrix z, IntegerVector locid, int ncommun, int nloc);
RcppExport SEXP _Rlda_getlk(SEXP zSEXP, SEXP locidSEXP, SEXP ncommunSEXP, SEXP nlocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type locid(locidSEXP);
    Rcpp::traits::input_parameter< int >::type ncommun(ncommunSEXP);
    Rcpp::traits::input_parameter< int >::type nloc(nlocSEXP);
    rcpp_result_gen = Rcpp::wrap(getlk(z, locid, ncommun, nloc));
    return rcpp_result_gen;
END_RCPP
}
// samplez
IntegerMatrix samplez(NumericMatrix ltheta, NumericMatrix l1minustheta, NumericMatrix lphi, IntegerMatrix dat1, IntegerVector locid, NumericMatrix randu, int ncommun, int nloc);
RcppExport SEXP _Rlda_samplez(SEXP lthetaSEXP, SEXP l1minusthetaSEXP, SEXP lphiSEXP, SEXP dat1SEXP, SEXP locidSEXP, SEXP randuSEXP, SEXP ncommunSEXP, SEXP nlocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ltheta(lthetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type l1minustheta(l1minusthetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lphi(lphiSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type dat1(dat1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type locid(locidSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type randu(randuSEXP);
    Rcpp::traits::input_parameter< int >::type ncommun(ncommunSEXP);
    Rcpp::traits::input_parameter< int >::type nloc(nlocSEXP);
    rcpp_result_gen = Rcpp::wrap(samplez(ltheta, l1minustheta, lphi, dat1, locid, randu, ncommun, nloc));
    return rcpp_result_gen;
END_RCPP
}
// convertVtoTheta
NumericMatrix convertVtoTheta(NumericMatrix vmat, NumericVector prod);
RcppExport SEXP _Rlda_convertVtoTheta(SEXP vmatSEXP, SEXP prodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type vmat(vmatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prod(prodSEXP);
    rcpp_result_gen = Rcpp::wrap(convertVtoTheta(vmat, prod));
    return rcpp_result_gen;
END_RCPP
}
// lda_binomial
List lda_binomial(DataFrame data, DataFrame pop, int n_community, double alpha0, double alpha1, double gamma, int n_gibbs, bool ll_prior, bool display_progress);
RcppExport SEXP _Rlda_lda_binomial(SEXP dataSEXP, SEXP popSEXP, SEXP n_communitySEXP, SEXP alpha0SEXP, SEXP alpha1SEXP, SEXP gammaSEXP, SEXP n_gibbsSEXP, SEXP ll_priorSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type pop(popSEXP);
    Rcpp::traits::input_parameter< int >::type n_community(n_communitySEXP);
    Rcpp::traits::input_parameter< double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type n_gibbs(n_gibbsSEXP);
    Rcpp::traits::input_parameter< bool >::type ll_prior(ll_priorSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_binomial(data, pop, n_community, alpha0, alpha1, gamma, n_gibbs, ll_prior, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// lda_binomial_var
Rcpp::List lda_binomial_var(arma::mat data, int n_community, int maxit, int n_obs, double gamma, double a_phi, double b_phi, double thresh, double delta_elbo, arma::cube m1, arma::cube m0);
RcppExport SEXP _Rlda_lda_binomial_var(SEXP dataSEXP, SEXP n_communitySEXP, SEXP maxitSEXP, SEXP n_obsSEXP, SEXP gammaSEXP, SEXP a_phiSEXP, SEXP b_phiSEXP, SEXP threshSEXP, SEXP delta_elboSEXP, SEXP m1SEXP, SEXP m0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n_community(n_communitySEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type a_phi(a_phiSEXP);
    Rcpp::traits::input_parameter< double >::type b_phi(b_phiSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< double >::type delta_elbo(delta_elboSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type m0(m0SEXP);
    rcpp_result_gen = Rcpp::wrap(lda_binomial_var(data, n_community, maxit, n_obs, gamma, a_phi, b_phi, thresh, delta_elbo, m1, m0));
    return rcpp_result_gen;
END_RCPP
}
// lda_bernoulli
List lda_bernoulli(DataFrame data, int n_community, double alpha0, double alpha1, double gamma, int n_gibbs, bool ll_prior, bool display_progress);
RcppExport SEXP _Rlda_lda_bernoulli(SEXP dataSEXP, SEXP n_communitySEXP, SEXP alpha0SEXP, SEXP alpha1SEXP, SEXP gammaSEXP, SEXP n_gibbsSEXP, SEXP ll_priorSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n_community(n_communitySEXP);
    Rcpp::traits::input_parameter< double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type n_gibbs(n_gibbsSEXP);
    Rcpp::traits::input_parameter< bool >::type ll_prior(ll_priorSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(lda_bernoulli(data, n_community, alpha0, alpha1, gamma, n_gibbs, ll_prior, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// convertSBtoNormal
NumericMatrix convertSBtoNormal(NumericMatrix vmat, int ncol, int nrow, NumericVector prod);
RcppExport SEXP _Rlda_convertSBtoNormal(SEXP vmatSEXP, SEXP ncolSEXP, SEXP nrowSEXP, SEXP prodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type vmat(vmatSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prod(prodSEXP);
    rcpp_result_gen = Rcpp::wrap(convertSBtoNormal(vmat, ncol, nrow, prod));
    return rcpp_result_gen;
END_RCPP
}
// aggregatesum
NumericVector aggregatesum(NumericVector Tobesum, int nind, int nobs, IntegerVector ind);
RcppExport SEXP _Rlda_aggregatesum(SEXP TobesumSEXP, SEXP nindSEXP, SEXP nobsSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Tobesum(TobesumSEXP);
    Rcpp::traits::input_parameter< int >::type nind(nindSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregatesum(Tobesum, nind, nobs, ind));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rlda_lda_multinomial", (DL_FUNC) &_Rlda_lda_multinomial, 7},
    {"_Rlda_getks", (DL_FUNC) &_Rlda_getks, 3},
    {"_Rlda_getlk", (DL_FUNC) &_Rlda_getlk, 4},
    {"_Rlda_samplez", (DL_FUNC) &_Rlda_samplez, 8},
    {"_Rlda_convertVtoTheta", (DL_FUNC) &_Rlda_convertVtoTheta, 2},
    {"_Rlda_lda_binomial", (DL_FUNC) &_Rlda_lda_binomial, 9},
    {"_Rlda_lda_binomial_var", (DL_FUNC) &_Rlda_lda_binomial_var, 11},
    {"_Rlda_lda_bernoulli", (DL_FUNC) &_Rlda_lda_bernoulli, 8},
    {"_Rlda_convertSBtoNormal", (DL_FUNC) &_Rlda_convertSBtoNormal, 4},
    {"_Rlda_aggregatesum", (DL_FUNC) &_Rlda_aggregatesum, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rlda(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
