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
  double lam01 = exp(par[0]);
  double theta = exp(par[1]);

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
double loglikS(DataFrame outdata_S, NumericVector par, List LAM03S, List cutS, Function fgau, Function fdpexp, Function fppexp){

  int i;
  double lam01 = exp(par[0]);
  double theta = exp(par[1]);

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

