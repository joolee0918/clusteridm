
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>
#include "commonf.h"

//#ifdef _OPENMP
#include <omp.h>
//#endif

using namespace Rcpp;



//[[Rcpp::export()]]
double NloglikFD1(NumericVector par, List outdata_F, NumericVector outdata_proband,
                  NumericVector Age, NumericVector Cal,  DataFrame lam03, bool full,
                 Function fgau, Function fdpexp, Function fppexp, Function combn){

  double lam01 = exp(par[0]);
  double rho = exp(par[1]);
  double newrho = rho/(1+rho);
  double rr;



  int i, j, k, nr, pair_nr;
  int nf = outdata_F.size();
  NumericVector tmp[2];
  NumericVector res(nf);
  double tmp11 = 0;
  double tmp1 =0, tmp2=0, tmp7=0;


  NumericVector pC0 = outdata_proband;

  for(i=0; i<nf; i++){

    DataFrame data = as<DataFrame>(outdata_F[i]);
    nr = data.nrow();
    pair_nr = (nr-1)*(nr-2)/2;
    NumericMatrix comb = combn(nr-1, 2);
    double C0;

    IntegerVector proband(nr);
    NumericVector B(nr);
    NumericVector X(nr);
    NumericVector Y(nr);
    NumericVector exam_age(nr);
    IntegerVector del1(nr);
    IntegerVector del2(nr);
    IntegerVector del3(nr);

    proband = data["proband"];
    B = data["B"];
    X = data["X"];
    Y = data["Y"];
    exam_age = data["exam.age"];
    del1 = data["del1"];
    del2 = data["del2"];
    del3 = data["del3"];
    C0 = pC0[i];

    IntegerVector idx(nr-1);
    idx = seq(1, nr-1);


    if(full == TRUE){
      NumericVector S = pexp(X, lam01, 0.0, 0.0);
      if(nr == 4){
        if (sum(del2)==4) {
          tmp1 +=log(hf(S, del2, 4, 4, rho));
        }  else if(sum(del2)==3) {
          tmp1 += log(hf(S, del2, 3, 4, rho));
        }   else if(sum(del2)==2) {
          tmp1 += log(hf(S, del2, 2, 4, rho));
        }  else if(sum(del2)==1) {
          tmp1 += log(hf(S, del2, 1, 4, rho));
        }
      } else{
        if (sum(del2)==6) {
          tmp1 += log(hf(S, del2, 6, 6, rho));
        }  else if(sum(del2)==5) {
          tmp1 += log(hf(S, del2, 5, 6, rho));
        }   else if(sum(del2)==4) {
          tmp1 += log(hf(S, del2, 4, 6, rho));
        }  else if(sum(del2)==3) {
          tmp1 += log(hf(S, del2, 3, 6, rho));
        }  else if(sum(del2)==2) {
          tmp1 += log(hf(S, del2, 2, 6, rho));
        }  else if(sum(del2)==1) {
          tmp1 += log(hf(S, del2, 1, 6, rho));
        }
      }
    }else{

    NumericMatrix Sf(nr-1, 2);
    NumericVector Sc(nr-1);
    NumericMatrix S(pair_nr, 2);
    IntegerMatrix fam_del2(pair_nr, 2);


    for(j=1; j<nr; j++){
      Sf(j-1, 0)  = R::pexp(X[0], 1/lam01, 0.0, 0.0);
      Sf(j-1, 1)  = R::pexp(X[j], 1/lam01, 0.0, 0.0);
      Sc[j-1]  = h1(R::pexp(X[0], 1/lam01, 0.0, 0.0), R::pexp(X[j], 1/lam01, 0.0, 0.0), rho);
    }

    for(j=0; j<pair_nr; j++){
      for(k=0; k<2; k++){
        S(j,k) =  Sc[comb(k,j)-1];
        fam_del2(j, k) = del2[comb(k,j)];
      }
       if (sum(fam_del2(j,_))==2) {
        tmp11 = dClayton(S(j,_),  newrho, 1.0) + dClayton(Sf(comb(0,j)-1,_), rho, 1.0) + dClayton(Sf(comb(1,j)-1,_), rho, 1.0);
      }  else if(sum(fam_del2(j,_))==1) {
        tmp11 = log(hf(S(j,_), fam_del2(j,_), 1, 2, newrho)) + fam_del2(j,0)*dClayton(Sf(comb(0,j)-1,_), rho, 1.0) + fam_del2(j,1)*dClayton(Sf(comb(1,j)-1,_), rho, 1.0);
      } else {
        tmp11 =  log(pClayton(S(j,_), newrho));
      }
      tmp1 += tmp11/(nr-2);
    }
    }

    for(j=0; j<nr; j++){
      tmp2 += del2[j]*R::dexp(X[j], 1/lam01, 1.0);
    }


    NumericMatrix gauss_quad = fgau(20, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);
    tmp7 += -log(sum(w*dexp(u, lam01)));
  }
  rr = tmp1 + tmp2 + tmp7;
  return(rr);
}




