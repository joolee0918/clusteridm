#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>

//#ifdef _OPENMP
#include <omp.h>
//#endif

#include "commonf.h"

//[[Rcpp::export()]]
double loglikR(DataFrame outdata_R, NumericVector par, List LAM03R, List cutR, Function fgau, Function fdpexp, Function fppexp){

  int i;
  double theta = exp(par[0]);
  double lam01 = exp(par[2]);


  double auxtmp1 = 0;

  for(i=0; i<LAM03R.size(); i++) {
    double C0 = as<NumericVector>(outdata_R["C0"])[i];
    double X = as<NumericVector>(outdata_R["X"])[i];
    double Y = as<NumericVector>(outdata_R["Y"])[i];
    int del3 = as<IntegerVector>(outdata_R["del3"])[i];

    NumericVector LAM03_R = LAM03R[i];
    NumericVector LAM12_R = theta*LAM03_R;
    NumericVector cut_R = cutR[i];;
    NumericMatrix gauss_quad = fgau(20, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);

    auxtmp1 += (del3*log(R::dexp(X, 1/lam01, 0.0)*as<NumericVector>(fppexp(X, LAM03_R, cut_R,  0.0, 0.0))[0]*as<NumericVector>(fdpexp(Y, LAM12_R, cut_R, 0.0))[0]/as<NumericVector>(fppexp(X, LAM12_R, cut_R,  0.0, 0.0))[0])
                  + (1-del3)*log(R::dexp(X, 1/lam01, 0.0)*as<NumericVector>(fppexp(X, LAM03_R, cut_R, 0.0,  0.0))[0]*as<NumericVector>(fppexp(Y, LAM12_R, cut_R, 0.0))[0]/as<NumericVector>(fppexp(X, LAM12_R, cut_R, 0.0, 0.0))[0])
                  -  log(sum(w*dexp(u, lam01)*as<NumericVector>(fppexp(u, LAM03_R, cut_R, 0.0, 0.0))*as<NumericVector>(fppexp(C0, LAM12_R, cut_R,  0.0, 0.0))[0]/as<NumericVector>(fppexp(u, LAM12_R, cut_R, 0.0, 0.0)))));
  }

  return(auxtmp1);

}



//[[Rcpp::export()]]
double NloglikR_pch(NumericVector par, NumericVector cut_F, NumericMatrix Y_R, NumericMatrix X_R,  Function fgau){

  int i;
  NumericVector lam01 = exp(par);


  double auxtmp1 = 0;

  for(i=0; i<Y_R.nrow(); i++) {
    double C0 = X_R(i, 0);
    double X = Y_R(i, 0);

    NumericMatrix gauss_quad = fgau(20, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);

    auxtmp1 += log(dpc(X, lam01, cut_F, 0.0)) -  log(sum(w*vdpc(u, lam01, cut_F, 0.0)));
  }

  return(auxtmp1);

}

//[[Rcpp::export()]]
double loglikR_pch(NumericVector par, NumericVector cut_F, NumericMatrix Y_R, NumericMatrix X_R,  List LAM03R, List cutR, Function fgau){

  int i;
  double theta = exp(par[0]);
  NumericVector lam01 = exp(par[seq(2, par.size()-1)]);


  double auxtmp1 = 0;

  for(i=0; i<LAM03R.size(); i++) {
    double C0 = X_R(i, 0);
    double X = Y_R(i, 0);
    double Y = Y_R(i, 1);
    int del3 = Y_R(i, 2);

    NumericVector LAM03_R = LAM03R[i];
    NumericVector LAM12_R = theta*LAM03_R;
    NumericVector cut_R = cutR[i];;
    NumericMatrix gauss_quad = fgau(20, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);

    auxtmp1 += (del3*log(dpc(X, lam01, cut_F, 0.0)*ppc(X, LAM03_R, cut_R,  0.0, 0.0)*dpc(Y, LAM12_R, cut_R, 0.0)/ppc(X, LAM12_R, cut_R,  0.0, 0.0))
                  + (1-del3)*log(dpc(X, lam01, cut_F, 0.0)*ppc(X, LAM03_R, cut_R, 0.0,  0.0)*ppc(Y, LAM12_R, cut_R, 0.0, 0.0)/ppc(X, LAM12_R, cut_R, 0.0, 0.0))
                  -  log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_R, cut_R, 0.0, 0.0)*ppc(C0, LAM12_R, cut_R,  0.0, 0.0)/vppc(u, LAM12_R, cut_R, 0.0, 0.0))));
  }

  return(auxtmp1);

}


//[[Rcpp::export()]]
double loglikS(NumericVector par, DataFrame outdata_S, List LAM03S, List cutS, Function fgau, Function fdpexp, Function fppexp){

  int i;
  double theta = exp(par[0]);
  double lam01 = exp(par[1]);


  double auxtmp2 = 0;


  for(i=0; i<LAM03S.size(); i++) {
    double C0 = as<NumericVector>(outdata_S["C"])[i];
    int del2 = as<NumericVector>(outdata_S["del2"])[i];

    NumericVector LAM03_S = LAM03S[i];
    NumericVector LAM12_S = theta*LAM03_S;
    NumericVector cut_S = cutS[i];

    NumericMatrix gauss_quad = fgau(20, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);

    auxtmp2 += ((1-del2)*log(exp(-lam01*C0)*as<NumericVector>(fppexp(C0, LAM03_S, cut_S, 0.0, 0.0))[0])
                  + del2*log(sum(w*dexp(u, lam01)*as<NumericVector>(fppexp(u, LAM03_S, cut_S, 0.0, 0.0))*as<NumericVector>(fppexp(C0, LAM12_S, cut_S,  0.0, 0.0))[0]/as<NumericVector>(fppexp(u, LAM12_S, cut_S, 0.0, 0.0))))
                  - log(exp(-lam01*C0)*as<NumericVector>(fppexp(C0, LAM03_S, cut_S, 0.0, 0.0))[0] + sum(w*dexp(u, lam01)*as<NumericVector>(fppexp(u, LAM03_S, cut_S, 0.0, 0.0))*as<NumericVector>(fppexp(C0, LAM12_S, cut_S, 0.0, 0.0))[0]/as<NumericVector>(fppexp(u, LAM12_S, cut_S, 0.0, 0.0)))));

  }

  return(auxtmp2);

}


//[[Rcpp::export()]]
double loglikS_pch( NumericVector par, NumericVector cut_F,  NumericMatrix Y_S, List LAM03S, List cutS, Function fgau){

  int i;
  double theta = exp(par[0]);
  NumericVector lam01 = exp(par[seq(2, par.size()-1)]);

  double auxtmp2 = 0;


  for(i=0; i<LAM03S.size(); i++) {
    double C0 = Y_S(i, 0);
    int del2 = Y_S(i,1);

    NumericVector LAM03_S = LAM03S[i];
    NumericVector LAM12_S = theta*LAM03_S;
    NumericVector cut_S = cutS[i];

    NumericMatrix gauss_quad = fgau(20, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);

    auxtmp2 += ((1-del2)*log(ppc(C0, lam01, cut_F, 0.0, 0.0)*ppc(C0, LAM03_S, cut_S, 0.0, 0.0))
                  + del2*log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_S, cut_S, 0.0, 0.0)*ppc(C0, LAM12_S, cut_S,  0.0, 0.0)/vppc(u, LAM12_S, cut_S, 0.0, 0.0)))
                  - log(ppc(C0, lam01, cut_F, 0.0, 0.0)*ppc(C0, LAM03_S, cut_S, 0.0, 0.0) + sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_S, cut_S, 0.0, 0.0)*ppc(C0, LAM12_S, cut_S, 0.0, 0.0)/vppc(u, LAM12_S, cut_S, 0.0, 0.0))));

  }

  return(auxtmp2);

}


//[[Rcpp::export()]]
double NloglikS_pch( NumericVector par, NumericVector cut_F,  NumericMatrix Y_S, Function fgau){

  int i;
  NumericVector lam01 = exp(par);

  double auxtmp2 = 0;


  for(i=0; i<Y_S.nrow(); i++) {
    double C0 = Y_S(i, 0);
    int del2 = Y_S(i,1);


    NumericMatrix gauss_quad = fgau(20, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);

    auxtmp2 += ((1-del2)*log(ppc(C0, lam01, cut_F, 0.0, 0.0))
                  + del2*log(sum(w*vdpc(u, lam01, cut_F, 0.0)))
                  - log(ppc(C0, lam01, cut_F, 0.0, 0.0)+ sum(w*vdpc(u, lam01, cut_F, 0.0))));

  }

  return(auxtmp2);

}



