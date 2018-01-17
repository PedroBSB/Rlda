#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _Rlda_aggregatesum(SEXP, SEXP, SEXP, SEXP);
extern SEXP _Rlda_convertSBtoNormal(SEXP, SEXP, SEXP, SEXP);
extern SEXP _Rlda_convertVtoTheta(SEXP, SEXP);
extern SEXP _Rlda_getks(SEXP, SEXP, SEXP);
extern SEXP _Rlda_getlk(SEXP, SEXP, SEXP, SEXP);
extern SEXP _Rlda_lda_bernoulli(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Rlda_lda_binomial(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Rlda_lda_binomial_var(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Rlda_lda_multinomial(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Rlda_samplez(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_Rlda_aggregatesum",      (DL_FUNC) &_Rlda_aggregatesum,       4},
  {"_Rlda_convertSBtoNormal", (DL_FUNC) &_Rlda_convertSBtoNormal,  4},
  {"_Rlda_convertVtoTheta",   (DL_FUNC) &_Rlda_convertVtoTheta,    2},
  {"_Rlda_getks",             (DL_FUNC) &_Rlda_getks,              3},
  {"_Rlda_getlk",             (DL_FUNC) &_Rlda_getlk,              4},
  {"_Rlda_lda_bernoulli",     (DL_FUNC) &_Rlda_lda_bernoulli,      8},
  {"_Rlda_lda_binomial",      (DL_FUNC) &_Rlda_lda_binomial,       9},
  {"_Rlda_lda_binomial_var",  (DL_FUNC) &_Rlda_lda_binomial_var,  11},
  {"_Rlda_lda_multinomial",   (DL_FUNC) &_Rlda_lda_multinomial,    7},
  {"_Rlda_samplez",           (DL_FUNC) &_Rlda_samplez,            8},
  {NULL, NULL, 0}
};

void R_init_Rlda(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
