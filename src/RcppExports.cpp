// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// hpc
NumericVector hpc(NumericVector x, NumericVector levels, NumericVector cuts, int logf);
RcppExport SEXP _clusteridm_hpc(SEXP xSEXP, SEXP levelsSEXP, SEXP cutsSEXP, SEXP logfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type levels(levelsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cuts(cutsSEXP);
    Rcpp::traits::input_parameter< int >::type logf(logfSEXP);
    rcpp_result_gen = Rcpp::wrap(hpc(x, levels, cuts, logf));
    return rcpp_result_gen;
END_RCPP
}
// Hpc
NumericVector Hpc(NumericVector x, NumericVector levels, NumericVector cuts, int logf);
RcppExport SEXP _clusteridm_Hpc(SEXP xSEXP, SEXP levelsSEXP, SEXP cutsSEXP, SEXP logfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type levels(levelsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cuts(cutsSEXP);
    Rcpp::traits::input_parameter< int >::type logf(logfSEXP);
    rcpp_result_gen = Rcpp::wrap(Hpc(x, levels, cuts, logf));
    return rcpp_result_gen;
END_RCPP
}
// ppc
double ppc(double q, NumericVector levels, NumericVector cuts, int lower, int logf);
RcppExport SEXP _clusteridm_ppc(SEXP qSEXP, SEXP levelsSEXP, SEXP cutsSEXP, SEXP lowerSEXP, SEXP logfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type levels(levelsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cuts(cutsSEXP);
    Rcpp::traits::input_parameter< int >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< int >::type logf(logfSEXP);
    rcpp_result_gen = Rcpp::wrap(ppc(q, levels, cuts, lower, logf));
    return rcpp_result_gen;
END_RCPP
}
// vppc
NumericVector vppc(NumericVector q, NumericVector levels, NumericVector cuts, int lower, int logf);
RcppExport SEXP _clusteridm_vppc(SEXP qSEXP, SEXP levelsSEXP, SEXP cutsSEXP, SEXP lowerSEXP, SEXP logfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type levels(levelsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cuts(cutsSEXP);
    Rcpp::traits::input_parameter< int >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< int >::type logf(logfSEXP);
    rcpp_result_gen = Rcpp::wrap(vppc(q, levels, cuts, lower, logf));
    return rcpp_result_gen;
END_RCPP
}
// dpc
double dpc(double x, NumericVector levels, NumericVector cuts, int logf);
RcppExport SEXP _clusteridm_dpc(SEXP xSEXP, SEXP levelsSEXP, SEXP cutsSEXP, SEXP logfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type levels(levelsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cuts(cutsSEXP);
    Rcpp::traits::input_parameter< int >::type logf(logfSEXP);
    rcpp_result_gen = Rcpp::wrap(dpc(x, levels, cuts, logf));
    return rcpp_result_gen;
END_RCPP
}
// vdpc
NumericVector vdpc(NumericVector x, NumericVector levels, NumericVector cuts, int logf);
RcppExport SEXP _clusteridm_vdpc(SEXP xSEXP, SEXP levelsSEXP, SEXP cutsSEXP, SEXP logfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type levels(levelsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cuts(cutsSEXP);
    Rcpp::traits::input_parameter< int >::type logf(logfSEXP);
    rcpp_result_gen = Rcpp::wrap(vdpc(x, levels, cuts, logf));
    return rcpp_result_gen;
END_RCPP
}
// order_cpp
NumericVector order_cpp(NumericVector x);
RcppExport SEXP _clusteridm_order_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(order_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// pG0
double pG0(arma::vec r_id, NumericVector G, double p);
RcppExport SEXP _clusteridm_pG0(SEXP r_idSEXP, SEXP GSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type r_id(r_idSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(pG0(r_id, G, p));
    return rcpp_result_gen;
END_RCPP
}
// pG
double pG(arma::vec r_id, NumericVector G, double p);
RcppExport SEXP _clusteridm_pG(SEXP r_idSEXP, SEXP GSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type r_id(r_idSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(pG(r_id, G, p));
    return rcpp_result_gen;
END_RCPP
}
// loglikFD2_pch
double loglikFD2_pch(NumericVector par, List Y_F, List X_F, NumericMatrix Y_proband, NumericMatrix X_proband, NumericVector Age, NumericVector Cal, NumericVector cut_F, DataFrame lam03, Function fgau, Function combn);
RcppExport SEXP _clusteridm_loglikFD2_pch(SEXP parSEXP, SEXP Y_FSEXP, SEXP X_FSEXP, SEXP Y_probandSEXP, SEXP X_probandSEXP, SEXP AgeSEXP, SEXP CalSEXP, SEXP cut_FSEXP, SEXP lam03SEXP, SEXP fgauSEXP, SEXP combnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< List >::type Y_F(Y_FSEXP);
    Rcpp::traits::input_parameter< List >::type X_F(X_FSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_proband(Y_probandSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X_proband(X_probandSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Age(AgeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Cal(CalSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cut_F(cut_FSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type lam03(lam03SEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    Rcpp::traits::input_parameter< Function >::type combn(combnSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikFD2_pch(par, Y_F, X_F, Y_proband, X_proband, Age, Cal, cut_F, lam03, fgau, combn));
    return rcpp_result_gen;
END_RCPP
}
// loglikFD2_pch_gene
double loglikFD2_pch_gene(NumericVector par, List Y_F, List X_F, NumericMatrix Y_proband, NumericMatrix X_proband, NumericVector Age, NumericVector Cal, NumericVector cut_F, DataFrame lam03, Function fgau, Function combn);
RcppExport SEXP _clusteridm_loglikFD2_pch_gene(SEXP parSEXP, SEXP Y_FSEXP, SEXP X_FSEXP, SEXP Y_probandSEXP, SEXP X_probandSEXP, SEXP AgeSEXP, SEXP CalSEXP, SEXP cut_FSEXP, SEXP lam03SEXP, SEXP fgauSEXP, SEXP combnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< List >::type Y_F(Y_FSEXP);
    Rcpp::traits::input_parameter< List >::type X_F(X_FSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_proband(Y_probandSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X_proband(X_probandSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Age(AgeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Cal(CalSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cut_F(cut_FSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type lam03(lam03SEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    Rcpp::traits::input_parameter< Function >::type combn(combnSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikFD2_pch_gene(par, Y_F, X_F, Y_proband, X_proband, Age, Cal, cut_F, lam03, fgau, combn));
    return rcpp_result_gen;
END_RCPP
}
// loglikR_pch
double loglikR_pch(NumericVector par, NumericVector cut_F, NumericMatrix Y_R, NumericMatrix X_R, List LAM03R, List cutR, Function fgau);
RcppExport SEXP _clusteridm_loglikR_pch(SEXP parSEXP, SEXP cut_FSEXP, SEXP Y_RSEXP, SEXP X_RSEXP, SEXP LAM03RSEXP, SEXP cutRSEXP, SEXP fgauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cut_F(cut_FSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_R(Y_RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X_R(X_RSEXP);
    Rcpp::traits::input_parameter< List >::type LAM03R(LAM03RSEXP);
    Rcpp::traits::input_parameter< List >::type cutR(cutRSEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikR_pch(par, cut_F, Y_R, X_R, LAM03R, cutR, fgau));
    return rcpp_result_gen;
END_RCPP
}
// loglikR_pch_gene
double loglikR_pch_gene(NumericVector par, NumericVector cut_F, NumericMatrix Y_R, NumericMatrix X_R, List LAM03R, List cutR, Function fgau);
RcppExport SEXP _clusteridm_loglikR_pch_gene(SEXP parSEXP, SEXP cut_FSEXP, SEXP Y_RSEXP, SEXP X_RSEXP, SEXP LAM03RSEXP, SEXP cutRSEXP, SEXP fgauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cut_F(cut_FSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_R(Y_RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X_R(X_RSEXP);
    Rcpp::traits::input_parameter< List >::type LAM03R(LAM03RSEXP);
    Rcpp::traits::input_parameter< List >::type cutR(cutRSEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikR_pch_gene(par, cut_F, Y_R, X_R, LAM03R, cutR, fgau));
    return rcpp_result_gen;
END_RCPP
}
// loglikS_pch
double loglikS_pch(NumericVector par, NumericVector cut_F, NumericMatrix Y_S, List LAM03S, List cutS, Function fgau);
RcppExport SEXP _clusteridm_loglikS_pch(SEXP parSEXP, SEXP cut_FSEXP, SEXP Y_SSEXP, SEXP LAM03SSEXP, SEXP cutSSEXP, SEXP fgauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cut_F(cut_FSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_S(Y_SSEXP);
    Rcpp::traits::input_parameter< List >::type LAM03S(LAM03SSEXP);
    Rcpp::traits::input_parameter< List >::type cutS(cutSSEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikS_pch(par, cut_F, Y_S, LAM03S, cutS, fgau));
    return rcpp_result_gen;
END_RCPP
}
// loglikS_pch_gene
double loglikS_pch_gene(NumericVector par, NumericVector cut_F, NumericMatrix Y_S, List LAM03S, List cutS, Function fgau);
RcppExport SEXP _clusteridm_loglikS_pch_gene(SEXP parSEXP, SEXP cut_FSEXP, SEXP Y_SSEXP, SEXP LAM03SSEXP, SEXP cutSSEXP, SEXP fgauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cut_F(cut_FSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_S(Y_SSEXP);
    Rcpp::traits::input_parameter< List >::type LAM03S(LAM03SSEXP);
    Rcpp::traits::input_parameter< List >::type cutS(cutSSEXP);
    Rcpp::traits::input_parameter< Function >::type fgau(fgauSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikS_pch_gene(par, cut_F, Y_S, LAM03S, cutS, fgau));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_clusteridm_hpc", (DL_FUNC) &_clusteridm_hpc, 4},
    {"_clusteridm_Hpc", (DL_FUNC) &_clusteridm_Hpc, 4},
    {"_clusteridm_ppc", (DL_FUNC) &_clusteridm_ppc, 5},
    {"_clusteridm_vppc", (DL_FUNC) &_clusteridm_vppc, 5},
    {"_clusteridm_dpc", (DL_FUNC) &_clusteridm_dpc, 4},
    {"_clusteridm_vdpc", (DL_FUNC) &_clusteridm_vdpc, 4},
    {"_clusteridm_order_cpp", (DL_FUNC) &_clusteridm_order_cpp, 1},
    {"_clusteridm_pG0", (DL_FUNC) &_clusteridm_pG0, 3},
    {"_clusteridm_pG", (DL_FUNC) &_clusteridm_pG, 3},
    {"_clusteridm_loglikFD2_pch", (DL_FUNC) &_clusteridm_loglikFD2_pch, 11},
    {"_clusteridm_loglikFD2_pch_gene", (DL_FUNC) &_clusteridm_loglikFD2_pch_gene, 11},
    {"_clusteridm_loglikR_pch", (DL_FUNC) &_clusteridm_loglikR_pch, 7},
    {"_clusteridm_loglikR_pch_gene", (DL_FUNC) &_clusteridm_loglikR_pch_gene, 7},
    {"_clusteridm_loglikS_pch", (DL_FUNC) &_clusteridm_loglikS_pch, 6},
    {"_clusteridm_loglikS_pch_gene", (DL_FUNC) &_clusteridm_loglikS_pch_gene, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_clusteridm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
