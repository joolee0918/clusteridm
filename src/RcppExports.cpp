// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// NloglikFD1
double NloglikFD1(NumericVector par, List outdata_F, NumericVector outdata_proband, NumericVector Age, NumericVector Cal, DataFrame lam03, bool full, Function fgau, Function fdpexp, Function fppexp, Function combn);
RcppExport SEXP _clusteridm_NloglikFD1(SEXP parSEXP, SEXP outdata_FSEXP, SEXP outdata_probandSEXP, SEXP AgeSEXP, SEXP CalSEXP, SEXP lam03SEXP, SEXP fullSEXP, SEXP fgauSEXP, SEXP fdpexpSEXP, SEXP fppexpSEXP, SEXP combnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< List >::type outdata_F(outdata_FSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type outdata_proband(outdata_probandSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Age(AgeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Cal(CalSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type lam03(lam03SEXP);
    Rcpp::traits::input_parameter< bool >::type full(fullSEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    Rcpp::traits::input_parameter< Function >::type fdpexp(fdpexpSEXP);
    Rcpp::traits::input_parameter< Function >::type fppexp(fppexpSEXP);
    Rcpp::traits::input_parameter< Function >::type combn(combnSEXP);
    rcpp_result_gen = Rcpp::wrap(NloglikFD1(par, outdata_F, outdata_proband, Age, Cal, lam03, full, fgau, fdpexp, fppexp, combn));
    return rcpp_result_gen;
END_RCPP
}
// NloglikFD2
double NloglikFD2(NumericVector par, List outdata_F, NumericVector outdata_proband, NumericVector Age, NumericVector Cal, DataFrame lam03, Function fgau, Function fdpexp, Function fppexp, Function combn, Function ppch);
RcppExport SEXP _clusteridm_NloglikFD2(SEXP parSEXP, SEXP outdata_FSEXP, SEXP outdata_probandSEXP, SEXP AgeSEXP, SEXP CalSEXP, SEXP lam03SEXP, SEXP fgauSEXP, SEXP fdpexpSEXP, SEXP fppexpSEXP, SEXP combnSEXP, SEXP ppchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< List >::type outdata_F(outdata_FSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type outdata_proband(outdata_probandSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Age(AgeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Cal(CalSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type lam03(lam03SEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    Rcpp::traits::input_parameter< Function >::type fdpexp(fdpexpSEXP);
    Rcpp::traits::input_parameter< Function >::type fppexp(fppexpSEXP);
    Rcpp::traits::input_parameter< Function >::type combn(combnSEXP);
    Rcpp::traits::input_parameter< Function >::type ppch(ppchSEXP);
    rcpp_result_gen = Rcpp::wrap(NloglikFD2(par, outdata_F, outdata_proband, Age, Cal, lam03, fgau, fdpexp, fppexp, combn, ppch));
    return rcpp_result_gen;
END_RCPP
}
// loglikFD1
double loglikFD1(NumericVector par, List outdata_F, NumericVector outdata_proband, NumericVector Age, NumericVector Cal, DataFrame lam03, bool full, Function fgau, Function fdpexp, Function fppexp, Function combn);
RcppExport SEXP _clusteridm_loglikFD1(SEXP parSEXP, SEXP outdata_FSEXP, SEXP outdata_probandSEXP, SEXP AgeSEXP, SEXP CalSEXP, SEXP lam03SEXP, SEXP fullSEXP, SEXP fgauSEXP, SEXP fdpexpSEXP, SEXP fppexpSEXP, SEXP combnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< List >::type outdata_F(outdata_FSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type outdata_proband(outdata_probandSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Age(AgeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Cal(CalSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type lam03(lam03SEXP);
    Rcpp::traits::input_parameter< bool >::type full(fullSEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    Rcpp::traits::input_parameter< Function >::type fdpexp(fdpexpSEXP);
    Rcpp::traits::input_parameter< Function >::type fppexp(fppexpSEXP);
    Rcpp::traits::input_parameter< Function >::type combn(combnSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikFD1(par, outdata_F, outdata_proband, Age, Cal, lam03, full, fgau, fdpexp, fppexp, combn));
    return rcpp_result_gen;
END_RCPP
}
// loglikFD2
double loglikFD2(NumericVector par, List outdata_F, NumericVector outdata_proband, NumericVector Age, NumericVector Cal, DataFrame lam03, Function fgau, Function fdpexp, Function fppexp, Function combn, Function ppch);
RcppExport SEXP _clusteridm_loglikFD2(SEXP parSEXP, SEXP outdata_FSEXP, SEXP outdata_probandSEXP, SEXP AgeSEXP, SEXP CalSEXP, SEXP lam03SEXP, SEXP fgauSEXP, SEXP fdpexpSEXP, SEXP fppexpSEXP, SEXP combnSEXP, SEXP ppchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< List >::type outdata_F(outdata_FSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type outdata_proband(outdata_probandSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Age(AgeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Cal(CalSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type lam03(lam03SEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    Rcpp::traits::input_parameter< Function >::type fdpexp(fdpexpSEXP);
    Rcpp::traits::input_parameter< Function >::type fppexp(fppexpSEXP);
    Rcpp::traits::input_parameter< Function >::type combn(combnSEXP);
    Rcpp::traits::input_parameter< Function >::type ppch(ppchSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikFD2(par, outdata_F, outdata_proband, Age, Cal, lam03, fgau, fdpexp, fppexp, combn, ppch));
    return rcpp_result_gen;
END_RCPP
}
// loglikR
double loglikR(DataFrame outdata_R, NumericVector par, List LAM03R, List cutR, Function fgau, Function fdpexp, Function fppexp);
RcppExport SEXP _clusteridm_loglikR(SEXP outdata_RSEXP, SEXP parSEXP, SEXP LAM03RSEXP, SEXP cutRSEXP, SEXP fgauSEXP, SEXP fdpexpSEXP, SEXP fppexpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type outdata_R(outdata_RSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< List >::type LAM03R(LAM03RSEXP);
    Rcpp::traits::input_parameter< List >::type cutR(cutRSEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    Rcpp::traits::input_parameter< Function >::type fdpexp(fdpexpSEXP);
    Rcpp::traits::input_parameter< Function >::type fppexp(fppexpSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikR(outdata_R, par, LAM03R, cutR, fgau, fdpexp, fppexp));
    return rcpp_result_gen;
END_RCPP
}
// loglikS
double loglikS(DataFrame outdata_S, NumericVector par, List LAM03S, List cutS, Function fgau, Function fdpexp, Function fppexp);
RcppExport SEXP _clusteridm_loglikS(SEXP outdata_SSEXP, SEXP parSEXP, SEXP LAM03SSEXP, SEXP cutSSEXP, SEXP fgauSEXP, SEXP fdpexpSEXP, SEXP fppexpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type outdata_S(outdata_SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< List >::type LAM03S(LAM03SSEXP);
    Rcpp::traits::input_parameter< List >::type cutS(cutSSEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    Rcpp::traits::input_parameter< Function >::type fdpexp(fdpexpSEXP);
    Rcpp::traits::input_parameter< Function >::type fppexp(fppexpSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikS(outdata_S, par, LAM03S, cutS, fgau, fdpexp, fppexp));
    return rcpp_result_gen;
END_RCPP
}
// NloglikS
double NloglikS(DataFrame outdata_S, NumericVector par, List LAM03S, List cutS, Function fgau, Function fdpexp, Function fppexp);
RcppExport SEXP _clusteridm_NloglikS(SEXP outdata_SSEXP, SEXP parSEXP, SEXP LAM03SSEXP, SEXP cutSSEXP, SEXP fgauSEXP, SEXP fdpexpSEXP, SEXP fppexpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type outdata_S(outdata_SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< List >::type LAM03S(LAM03SSEXP);
    Rcpp::traits::input_parameter< List >::type cutS(cutSSEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    Rcpp::traits::input_parameter< Function >::type fdpexp(fdpexpSEXP);
    Rcpp::traits::input_parameter< Function >::type fppexp(fppexpSEXP);
    rcpp_result_gen = Rcpp::wrap(NloglikS(outdata_S, par, LAM03S, cutS, fgau, fdpexp, fppexp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_clusteridm_NloglikFD1", (DL_FUNC) &_clusteridm_NloglikFD1, 11},
    {"_clusteridm_NloglikFD2", (DL_FUNC) &_clusteridm_NloglikFD2, 11},
    {"_clusteridm_loglikFD1", (DL_FUNC) &_clusteridm_loglikFD1, 11},
    {"_clusteridm_loglikFD2", (DL_FUNC) &_clusteridm_loglikFD2, 11},
    {"_clusteridm_loglikR", (DL_FUNC) &_clusteridm_loglikR, 7},
    {"_clusteridm_loglikS", (DL_FUNC) &_clusteridm_loglikS, 7},
    {"_clusteridm_NloglikS", (DL_FUNC) &_clusteridm_NloglikS, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_clusteridm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
