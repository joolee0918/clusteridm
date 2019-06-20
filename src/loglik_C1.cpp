
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
double loglikFD1(NumericVector par, List outdata_F, NumericVector outdata_proband,
                NumericVector Age, NumericVector Cal,  DataFrame lam03, bool full,
                Function fgau, Function fdpexp, Function fppexp, Function combn){

  double theta = exp(par[0]);
  double rho = exp(par[1]);
  double lam01 = exp(par[2]);

  double newrho = rho/(1+rho);



int i, j, k, l, nr, pair_nr;
int nf = outdata_F.size();
NumericVector tmp[2];
int ncal = Cal.size();
int nage = Age.size();
NumericVector res(nf);
double tmp11 = 0;
double tmp1 =0;


NumericVector pC0 = outdata_proband;
IntegerVector lam03_Y = lam03["Year.f"];
IntegerVector lam03_A = lam03["Age.f"];
NumericVector lam03_rate = lam03["rate"];
NumericVector lam12_rate = theta*as<NumericVector>(lam03["rate"]);

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


  IntegerVector idx(nr-1);
  idx = seq(1, nr-1);

  C0 = pC0[i];

  List LAM03(nr), LAM12(nr), cut(nr);
  NumericVector cc(ncal + nage);

  for(j=0; j<nr; j++){
    NumericVector vt(nage);

    for(l=0; l<nage; l++){
      vt[l] = B[j] + Age[l];
    }

    cc = union_(Cal, vt);
    NumericVector C = sort_unique(cc);
    NumericVector Cf = C[B[j] <= C & C <= (exam_age[j] + B[j])];

    int ncf = Cf.size();
    NumericVector Af(ncf);
    for(l=0; l<ncf; l++){
      Af[l] = Cf[l] - B[j];
    }

    IntegerVector R = findInterval(Cf, Cal);
    IntegerVector A = findInterval(Af, Age);
    NumericVector LAM03_rate(R.size()), LAM12_rate(R.size());

    for(l=0; l<R.size(); l++ ){
      NumericVector lam03rate = lam03_rate[(lam03_Y==R[l] & lam03_A==A[l])];
      NumericVector lam12rate = lam12_rate[(lam03_Y==R[l] & lam03_A==A[l])];
      LAM03_rate[l] = lam03rate[0];
      LAM12_rate[l] = lam12rate[0];
    }

    LAM03[j] = LAM03_rate;
    LAM12[j] = LAM12_rate;
    cut[j] = Af;

  }


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

  for(j=1; j<nr; j++){
    tmp1 += del2[j]*R::dexp(X[j], 1/lam01, 1.0);
    tmp1 += as<NumericVector>(fppexp(X[j], LAM03[j], cut[j], 0.0, 1.0))[0];
    tmp1 += del1[j]*(1-del2[j])*log(as<NumericVector>(fdpexp(X[j], LAM03[j], cut[j], 0.0))[0]/as<NumericVector>(fppexp(X[j], LAM03[j], cut[j], 0.0, 0.0))[0]);
    tmp1 += del2[j]*log(as<NumericVector>(fppexp(Y[j], LAM12[j], cut[j], 0.0, 0.0))[0]/as<NumericVector>(fppexp(X[j], LAM12[j], cut[j], 0.0, 0.0))[0]);
    tmp1 += del3[j]*log(as<NumericVector>(fdpexp(Y[j], LAM12[j], cut[j], 0.0))[0]/as<NumericVector>(fppexp(Y[j], LAM12[j], cut[j], 0.0, 0.0))[0]);

  }


  tmp1 +=  R::dexp(X[0], 1/lam01, 1.0);
  tmp1 += as<NumericVector>(fppexp(X[0], LAM03[0], cut[0], 0.0, 1.0))[0];
  tmp1 += log(as<NumericVector>(fppexp(Y[0], LAM12[0], cut[0], 0.0, 0.0))[0]/as<NumericVector>(fppexp(X[0], LAM12[0], cut[0], 0.0, 0.0))[0]);


 NumericMatrix gauss_quad = fgau(20, 0, C0);
  NumericVector u = gauss_quad(_,0);
  NumericVector w = gauss_quad(_,1);
  tmp1 += -log(sum(w*dexp(u, lam01)*as<NumericVector>(fppexp(u, LAM03[0], cut[0], 0.0, 0.0))*as<NumericVector>(fppexp(C0, LAM12[0], cut[0],  0.0, 0.0))[0]/as<NumericVector>(fppexp(u, LAM12[0], cut[0], 0.0, 0.0))));
  }
  return(tmp1);
}






//[[Rcpp::export()]]
double loglikFD1_pch(NumericVector par, NumericVector cut_F, List outdata_F, NumericVector outdata_proband,
                 NumericVector Age, NumericVector Cal,  DataFrame lam03, bool full,
                 Function fgau, Function combn){

  double theta = exp(par[0]);
  double rho = exp(par[1]);
  double newrho = rho/(1+rho);
  NumericVector lam01 = exp(par[seq(2, par.size()-1)]);


  int i, j, k, l, nr, pair_nr;
  int nf = outdata_F.size();
  NumericVector tmp[2];
  int ncal = Cal.size();
  int nage = Age.size();
  NumericVector res(nf);
  double tmp11 = 0;
  double tmp1 =0;


  NumericVector pC0 = outdata_proband;
  IntegerVector lam03_Y = lam03["Year.f"];
  IntegerVector lam03_A = lam03["Age.f"];
  NumericVector lam03_rate = lam03["rate"];
  NumericVector lam12_rate = theta*as<NumericVector>(lam03["rate"]);

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


    IntegerVector idx(nr-1);
    idx = seq(1, nr-1);

    C0 = pC0[i];

    List LAM03(nr), LAM12(nr), cut(nr);
    NumericVector cc(ncal + nage);

    for(j=0; j<nr; j++){
      NumericVector vt(nage);

      for(l=0; l<nage; l++){
        vt[l] = B[j] + Age[l];
      }

      cc = union_(Cal, vt);
      NumericVector C = sort_unique(cc);
      NumericVector Cf = C[B[j] <= C & C <= (exam_age[j] + B[j])];

      int ncf = Cf.size();
      NumericVector Af(ncf);
      for(l=0; l<ncf; l++){
        Af[l] = Cf[l] - B[j];
      }

      IntegerVector R = findInterval(Cf, Cal);
      IntegerVector A = findInterval(Af, Age);
      NumericVector LAM03_rate(R.size()), LAM12_rate(R.size());

      for(l=0; l<R.size(); l++ ){
        NumericVector lam03rate = lam03_rate[(lam03_Y==R[l] & lam03_A==A[l])];
        NumericVector lam12rate = lam12_rate[(lam03_Y==R[l] & lam03_A==A[l])];
        LAM03_rate[l] = lam03rate[0];
        LAM12_rate[l] = lam12rate[0];
      }

      LAM03[j] = LAM03_rate;
      LAM12[j] = LAM12_rate;
      cut[j] = Af[seq(1, Af.size()-1)];

    }


    if(full == TRUE){
      NumericVector S = vppc(X, lam01, cut_F, 0.0, 0.0);
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
        Sf(j-1, 0)  = ppc(X[0], lam01, cut_F, 0.0, 0.0);
        Sf(j-1, 1)  = ppc(X[j], lam01, cut_F, 0.0, 0.0);
        Sc[j-1]  = h1(ppc(X[0], lam01, cut_F, 0.0, 0.0), ppc(X[j], lam01, cut_F, 0.0, 0.0), rho);
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

    for(j=1; j<nr; j++){
      tmp1 += del2[j]*dpc(X[j], lam01, cut_F, 1.0);
      tmp1 += ppc(X[j], LAM03[j], cut[j],  0.0, 1.0);
      tmp1 += del1[j]*(1-del2[j])*log(dpc(X[j], LAM03[j], cut[j], 0.0)/ppc(X[j], LAM03[j], cut[j], 0.0, 0.0));
      tmp1 += del2[j]*log(ppc(Y[j], LAM12[j], cut[j], 0.0, 0.0)/ppc(X[j], LAM12[j], cut[j], 0.0, 0.0));
      tmp1 += del3[j]*log(dpc(Y[j], LAM12[j], cut[j], 0.0)/ppc(Y[j], LAM12[j], cut[j], 0.0, 0.0));

    }


    tmp1 +=  dpc(X[0], lam01, cut_F, 1.0);
    tmp1 += ppc(X[0], LAM03[0], cut[0], 0.0, 1.0);
    tmp1 += log(ppc(Y[0], LAM12[0], cut[0], 0.0, 0.0)/ppc(X[0], LAM12[0], cut[0], 0.0, 0.0));


    NumericMatrix gauss_quad = fgau(20, 0, C0);
    NumericVector u = gauss_quad(_,0);
    NumericVector w = gauss_quad(_,1);
    tmp1 += -log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03[0], cut[0], 0.0, 0.0)*ppc(C0, LAM12[0], cut[0],  0.0, 0.0)/vppc(u, LAM12[0], cut[0], 0.0, 0.0)));
  }
  return(tmp1);
}

