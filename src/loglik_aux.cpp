// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
double loglikR_pch(NumericVector par, NumericVector cut_F, NumericMatrix Y_R, NumericMatrix X_R,  List LAM03R, List LAM12R, List cutR, Function fgau){

  int i;
  //double theta = exp(par[1]);
  NumericVector lam01 = exp(par[seq(4, par.size()-1)]);


  double auxtmp1 = 0;

  for(i=0; i<LAM03R.size(); i++) {
    double C0 = X_R(i, 1);
    double X = Y_R(i, 0);
    double Y = Y_R(i, 1);
    int del3 = Y_R(i, 2);

    NumericVector LAM03_R = LAM03R[i];
    NumericVector LAM12_R = LAM12R[i];
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
double loglikR_pch_gene(NumericVector par, NumericVector cut_F, NumericMatrix Y_R, NumericMatrix X_R,   List LAM03R, List LAM12R, List cutR, Function fgau){

  int i;
  double pg0R;
 // double theta = exp(par[1]);
  double alpha = par[4];
  double p = exp(par[5]);
  NumericVector lam01 = exp(par[seq(6, par.size()-1)]);



  double auxtmp1 = 0;

  for(i=0; i<LAM03R.size(); i++) {
    double C0 = X_R(i, 1);
    double X = Y_R(i, 0);
    double Y = Y_R(i, 1);
    int del3 = Y_R(i, 2);
    double IG = X_R(i, 3);

    if(IG==0){
      pg0R =  pow((1-p),2);
    } else{
      pg0R = 1- pow((1-p),2);
    }



    NumericVector LAM03_R = LAM03R[i];
    NumericVector LAM12_R = LAM12R[i];
    NumericVector cut_R = cutR[i];;
    NumericMatrix gauss_quad = fgau(40, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);

    if(del3 == 1){
      if(IG== -999){
        auxtmp1 += (log(dpc(X, lam01, cut_F, 0.0)*ppc(X, LAM03_R, cut_R,  0.0, 0.0)*dpc(Y, LAM12_R, cut_R, 0.0)/ppc(X, LAM12_R, cut_R,  0.0, 0.0)*pow((1-p),2)
                        + dpc(X, lam01*exp(alpha), cut_F, 0.0)*ppc(X, LAM03_R, cut_R,  0.0, 0.0)*dpc(Y, LAM12_R, cut_R, 0.0)/ppc(X, LAM12_R, cut_R,  0.0, 0.0)*(1-pow((1-p),2)))
                      -  log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_R, cut_R, 0.0, 0.0)*ppc(C0, LAM12_R, cut_R,  0.0, 0.0)/vppc(u, LAM12_R, cut_R, 0.0, 0.0))*pow((1-p),2)
                          + sum(w*vdpc(u, lam01*exp(alpha), cut_F, 0.0)*vppc(u, LAM03_R, cut_R, 0.0, 0.0)*ppc(C0, LAM12_R, cut_R,  0.0, 0.0)/vppc(u, LAM12_R, cut_R, 0.0, 0.0))*(1-pow((1-p),2))));

      } else{

        auxtmp1 += (log(dpc(X, lam01*exp(alpha*IG), cut_F, 0.0)*ppc(X, LAM03_R, cut_R,  0.0, 0.0)*dpc(Y, LAM12_R, cut_R, 0.0)/ppc(X, LAM12_R, cut_R,  0.0, 0.0))
                      + log(pg0R)  -  log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_R, cut_R, 0.0, 0.0)*ppc(C0, LAM12_R, cut_R,  0.0, 0.0)/vppc(u, LAM12_R, cut_R, 0.0, 0.0))*pow((1-p),2)
                      + sum(w*vdpc(u, lam01*exp(alpha), cut_F, 0.0)*vppc(u, LAM03_R, cut_R, 0.0, 0.0)*ppc(C0, LAM12_R, cut_R,  0.0, 0.0)/vppc(u, LAM12_R, cut_R, 0.0, 0.0))*(1-pow((1-p),2)) ));

      }
    }else{
        if(IG== -999){
          auxtmp1 += (log(dpc(X, lam01, cut_F, 0.0)*ppc(X, LAM03_R, cut_R, 0.0,  0.0)*ppc(Y, LAM12_R, cut_R, 0.0, 0.0)/ppc(X, LAM12_R, cut_R, 0.0, 0.0)*pow((1-p),2)
                            + dpc(X, lam01*exp(alpha), cut_F, 0.0)*ppc(X, LAM03_R, cut_R, 0.0,  0.0)*ppc(Y, LAM12_R, cut_R, 0.0, 0.0)/ppc(X, LAM12_R, cut_R, 0.0, 0.0)*(1-pow((1-p),2)))
                        -  log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_R, cut_R, 0.0, 0.0)*ppc(C0, LAM12_R, cut_R,  0.0, 0.0)/vppc(u, LAM12_R, cut_R, 0.0, 0.0))*pow((1-p),2)
                        + sum(w*vdpc(u, lam01*exp(alpha), cut_F, 0.0)*vppc(u, LAM03_R, cut_R, 0.0, 0.0)*ppc(C0, LAM12_R, cut_R,  0.0, 0.0)/vppc(u, LAM12_R, cut_R, 0.0, 0.0))*(1-pow((1-p),2))));

        } else{
         auxtmp1 += (log(dpc(X, lam01*exp(alpha*IG), cut_F, 0.0)*ppc(X, LAM03_R, cut_R, 0.0,  0.0)*ppc(Y, LAM12_R, cut_R, 0.0, 0.0)/ppc(X, LAM12_R, cut_R, 0.0, 0.0))
                       + log(pg0R) -  log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_R, cut_R, 0.0, 0.0)*ppc(C0, LAM12_R, cut_R,  0.0, 0.0)/vppc(u, LAM12_R, cut_R, 0.0, 0.0))*pow((1-p),2)
                        + sum(w*vdpc(u, lam01*exp(alpha), cut_F, 0.0)*vppc(u, LAM03_R, cut_R, 0.0, 0.0)*ppc(C0, LAM12_R, cut_R,  0.0, 0.0)/vppc(u, LAM12_R, cut_R, 0.0, 0.0))*(1-pow((1-p),2))));

      }
      }

  }

  return(auxtmp1);

}


//[[Rcpp::export()]]
double loglikS_pch( NumericVector par, NumericVector cut_F,  NumericMatrix Y_S, List LAM03S, List LAM12S, List cutS, Function fgau){

  int i;
  //double theta = exp(par[1]);
  NumericVector lam01 = exp(par[seq(4, par.size()-1)]);

  double auxtmp2 = 0;


  for(i=0; i<LAM03S.size(); i++) {
    double C0 = Y_S(i, 0);
    int del2 = Y_S(i,1);

    NumericVector LAM03_S = LAM03S[i];
    NumericVector LAM12_S = LAM12S[i];
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
double loglikS_pch_gene(NumericVector par, NumericVector cut_F,  NumericMatrix Y_S, List LAM03S, List LAM12S, List cutS, Function fgau){

  int i;
  //double theta = exp(par[1]);
  double alpha = par[4];
  double p = exp(par[5]);
  NumericVector lam01 = exp(par[seq(6, par.size()-1)]);

  double auxtmp2 = 0;


  for(i=0; i<LAM03S.size(); i++) {
    double C0 = Y_S(i, 0);
    int del2 = Y_S(i,1);

    NumericVector LAM03_S = LAM03S[i];
    NumericVector LAM12_S = LAM12S[i];
    NumericVector cut_S = cutS[i];

    NumericMatrix gauss_quad = fgau(20, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);

    if(del2 == 0){
      auxtmp2 += (log(ppc(C0, lam01, cut_F, 0.0, 0.0)*ppc(C0, LAM03_S, cut_S, 0.0, 0.0)*pow((1-p),2)
                    +ppc(C0, lam01*exp(alpha), cut_F, 0.0, 0.0)*ppc(C0, LAM03_S, cut_S, 0.0, 0.0)*(1-pow((1-p),2)))
                    - log((ppc(C0, lam01, cut_F, 0.0, 0.0)*ppc(C0, LAM03_S, cut_S, 0.0, 0.0) + sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_S, cut_S, 0.0, 0.0)*ppc(C0, LAM12_S, cut_S, 0.0, 0.0)/vppc(u, LAM12_S, cut_S, 0.0, 0.0)))*pow((1-p),2)
                    + (ppc(C0, lam01*exp(alpha), cut_F, 0.0, 0.0)*ppc(C0, LAM03_S, cut_S, 0.0, 0.0) + sum(w*vdpc(u, lam01*exp(alpha), cut_F, 0.0)*vppc(u, LAM03_S, cut_S, 0.0, 0.0)*ppc(C0, LAM12_S, cut_S, 0.0, 0.0)/vppc(u, LAM12_S, cut_S, 0.0, 0.0)))*(1-pow((1-p),2))));

    }else{
      auxtmp2 += (log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_S, cut_S, 0.0, 0.0)*ppc(C0, LAM12_S, cut_S,  0.0, 0.0)/vppc(u, LAM12_S, cut_S, 0.0, 0.0))*pow((1-p),2)
                    +sum(w*vdpc(u, lam01*exp(alpha), cut_F, 0.0)*vppc(u, LAM03_S, cut_S, 0.0, 0.0)*ppc(C0, LAM12_S, cut_S,  0.0, 0.0)/vppc(u, LAM12_S, cut_S, 0.0, 0.0))*(1-pow((1-p),2)))
                    - log((ppc(C0, lam01, cut_F, 0.0, 0.0)*ppc(C0, LAM03_S, cut_S, 0.0, 0.0) + sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03_S, cut_S, 0.0, 0.0)*ppc(C0, LAM12_S, cut_S, 0.0, 0.0)/vppc(u, LAM12_S, cut_S, 0.0, 0.0)))*pow((1-p),2)
                    + (ppc(C0, lam01*exp(alpha), cut_F, 0.0, 0.0)*ppc(C0, LAM03_S, cut_S, 0.0, 0.0) + sum(w*vdpc(u, lam01*exp(alpha), cut_F, 0.0)*vppc(u, LAM03_S, cut_S, 0.0, 0.0)*ppc(C0, LAM12_S, cut_S, 0.0, 0.0)/vppc(u, LAM12_S, cut_S, 0.0, 0.0)))*(1-pow((1-p),2))));

    }

  }

  return(auxtmp2);

}


