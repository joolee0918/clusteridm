#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>

//#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
//#endif

#include "commonf.h"

using namespace Rcpp;


//[[Rcpp::export()]]
double loglikFD2(NumericVector par, List outdata_F, NumericVector outdata_proband,
               NumericVector Age, NumericVector Cal,  DataFrame lam03,
               Function fgau, Function fdpexp, Function fppexp, Function combn, Function ppch){

  double theta = exp(par[1]);
  double rho = par[0];
  double lam01 = exp(par[2]);
  double newrho = rho/(1+rho);


  int i, j, j1, j2, k, l, nr, pair_nr;
  int nf = outdata_F.size();
  NumericVector tmp(2);
  int ncal = Cal.size();
  int nage = Age.size();
  double C0;
  double tmp11 = 0;
  double tmp1=0, tmp2=0, tmp3=0, tmp5=0, tmp7=0, tmp8=0;
  double rr=0;

  NumericVector pC0 = outdata_proband;
  IntegerVector lam03_Y = lam03["Year.f"];
  IntegerVector lam03_A = lam03["Age.f"];
  NumericVector lam03_rate = lam03["rate"];
  NumericVector lam12_rate = theta*as<NumericVector>(lam03["rate"]);

  NumericMatrix gauss_quad(20, 2);
  NumericVector u(20), u1(20), u2(20);
  NumericVector w(20), w1(20), w2(20);
  IntegerVector sq(20);
  IntegerVector rep1(20*20);
  IntegerVector rep2(20*20);
  NumericVector uu1(20*20);
  NumericVector uu2(20*20);
  NumericVector ww1(20*20);
  NumericVector ww2(20*20);

  for(i=0; i<nf; i++){
    tmp1 = tmp2 = tmp3 = tmp5 = tmp7 = tmp8 = 0;

    DataFrame data = as<DataFrame>(outdata_F[i]);
    nr = data.nrow();


    IntegerVector proband(nr);
    NumericVector B(nr);
    NumericVector X(nr);
    NumericVector Y(nr);
    NumericVector exam_age(nr);
    IntegerVector del1(nr);
    IntegerVector del2(nr);
    IntegerVector del3(nr);


    B = data["B"];
    X = data["X"];
    Y = data["Y"];
    exam_age = data["exam.age"];

    del2 = data["del2"];
    del3 = data["del3"];

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

    if(nr==2){

      NumericVector S0 = pexp(X, lam01, 0.0);
      if(sum(del2)==2) tmp1 += dClayton(S0, rho, 1.0);
      else tmp1 += log(hf(S0, del2, 1, 2, rho));

      IntegerVector tmpdel = IntegerVector::create(1,0);

      tmp[0] = X[0];
      tmp[1] = exam_age[1];

      gauss_quad = fgau(20, 0, exam_age[1]);
      u = gauss_quad(_,0);
      w = gauss_quad(_,1);

      double px = R::pexp(X[0], 1/lam01, 0.0, 0.0);
      NumericVector vpx = rep(px, 20);

      tmp8 +=  -log(hf(pexp(tmp, lam01, 0.0, 0.0), tmpdel, 1, 2, rho)*as<NumericVector>(fppexp(exam_age[1], LAM03[1], cut[1], 0.0))[0]
              + sum(w*vdClayton(vpx, pexp(u, lam01, 0.0, 0.0), rho, 0.0)*dexp(u, lam01)*as<NumericVector>(fppexp(u, LAM03[1], cut[1], 0.0))*as<NumericVector>(fppexp(exam_age[1], LAM12[1], cut[1], 0.0))[0]/as<NumericVector>(fppexp(u, LAM12[1], cut[1], 0.0))));


    }else{
    pair_nr = (nr-1)*(nr-2)/2;
    NumericMatrix comb(2, pair_nr);
    comb = as<NumericMatrix>(combn(nr-1, 2));

    NumericMatrix Sf(nr-1, 2);
    NumericMatrix SAf(nr-1, 2);
    NumericVector Sc(nr-1);
    NumericVector SAc(nr-1);
    NumericMatrix S(pair_nr, 2);
    NumericMatrix SA(pair_nr, 2);
    IntegerMatrix fam_del2(pair_nr, 2);



    for(j=1; j<nr; j++){
      Sf(j-1, 0)  = SAf(j-1, 0) = R::pexp(X[0], 1/lam01, 0.0, 0.0);
      Sf(j-1, 1)  = R::pexp(X[j], 1/lam01, 0.0, 0.0);
      SAf(j-1, 1)  = R::pexp(exam_age[j], 1/lam01, 0.0, 0.0);
      Sc[j-1]  = h1(R::pexp(X[0], 1/lam01, 0.0, 0.0), R::pexp(X[j], 1/lam01, 0.0, 0.0), rho);
      SAc[j-1]  = h1(R::pexp(X[0], 1/lam01, 0.0, 0.0), R::pexp(exam_age[j], 1/lam01, 0.0, 0.0), rho);
    }


    for(j=0; j<pair_nr; j++){
      for(k=0; k<2; k++){
        S(j,k) =  Sc[comb(k,j)-1];
        SA(j,k) =  SAc[comb(k,j)-1];
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



    double px = R::pexp(X[0], 1/lam01, 0.0, 0.0);
    NumericVector vpx = rep(px, 20);
    NumericVector vpx2 = rep(px, 20*20);


     for(j=0; j<pair_nr; j++){
       j1 = comb(0,j);
       j2 = comb(1,j);
       NumericVector tmpcut1 = cut[j1];
       NumericVector cut1(tmpcut1.size()-1);
       for (k=1; k<tmpcut1.size(); k++)
         cut1[k-1] = tmpcut1[k];

       NumericVector tmpcut2 = cut[j2];
       NumericVector cut2(tmpcut2.size()-1);
       for (k=1; k<tmpcut2.size(); k++)
         cut2[k-1] = tmpcut2[k];


      gauss_quad = fgau(20, 0, exam_age[j1]);
      u1 = gauss_quad(_,0);
      w1 = gauss_quad(_,1);

      gauss_quad = fgau(20, 0, exam_age[j2]);
      u2 = gauss_quad(_,0);
      w2 = gauss_quad(_,1);
      sq = seq_len(20)-1;
      rep1 = rep_each(sq, 20);
      rep2 = rep(sq, 20);
      uu1 = u1[rep1];
      uu2 = u2[rep2];
      ww1 = w1[rep1];
      ww2 = w2[rep2];

     tmp8 -= log(pClayton(SA(j,_), newrho)*as<NumericVector>(fppexp(exam_age[j1], LAM03[j1], cut[j1], 0.0))[0]*as<NumericVector>(fppexp(exam_age[j2], LAM03[j2], cut[j2], 0.0))[0]
      + sum(w2*vh1(vh1(vpx, pexp(u2, lam01, 0.0), rho), rep(SA(j,0), 20), newrho)
       *vdClayton(vpx, pexp(u2, lam01, 0.0), rho, 0.0)*dexp(u2, lam01)*as<NumericVector>(ppch(u2, cut2, LAM03[j2], 0.0))
       /as<NumericVector>(ppch(u2, cut2, LAM12[j2],  0.0)))*as<NumericVector>(ppch(exam_age[j2], cut2, LAM12[j2], 0.0))[0]
       *as<NumericVector>(ppch(exam_age[j1], cut1, LAM03[j1],  0.0))[0]
       + sum(w1*vh1(vh1(vpx, pexp(u1, lam01, 0.0), rho), rep(SA(j,1), 20), newrho)
        *vdClayton(vpx, pexp(u1, lam01, 0.0), rho, 0.0)*dexp(u1, lam01)*as<NumericVector>(ppch(u1, cut1, LAM03[j1], 0.0))
        /as<NumericVector>(ppch(u1, cut1, LAM12[j1],  0.0)))*as<NumericVector>(ppch(exam_age[j1], cut1, LAM12[j1], 0.0))[0]
        *as<NumericVector>(ppch(exam_age[j2], cut2, LAM03[j2],  0.0))[0]
       +sum(ww1*ww2*vdClayton(vh1(vpx2, pexp(uu1, lam01, 0.0), rho) , vh1(vpx2, pexp(uu2, lam01, 0.0), rho),  newrho, 0.0)
       *vdClayton(vpx2, pexp(uu1, lam01, 0.0), rho, 0.0)*vdClayton(vpx2, pexp(uu2, lam01, 0.0), rho, 0.0)
       *dexp(uu1, lam01)*dexp(uu2, lam01)*as<NumericVector>(ppch(uu1, cut1, LAM03[j1], 0.0))*as<NumericVector>(ppch(exam_age[j1], cut1, LAM12[j1], 0.0))[0]
       /as<NumericVector>(ppch(uu1, cut1, LAM12[j1],0.0))*as<NumericVector>(ppch(uu2, cut2, LAM03[j2], 0.0))
       *as<NumericVector>(ppch(exam_age[j2], cut2, LAM12[j2], 0.0))[0]/as<NumericVector>(ppch(uu2, cut2,  LAM12[j2], 0.0))))/(nr-2);

     }

    }

    for(j=0; j<nr; j++){
      tmp2 += del2[j]*R::dexp(X[j], 1/lam01, 1.0);
      tmp3 += as<NumericVector>(fppexp(X[j], LAM03[j], cut[j], 0.0, 1.0))[0];
      tmp5 += del2[j]*log(as<NumericVector>(fppexp(Y[j], LAM12[j], cut[j], 0.0, 0.0))[0]/as<NumericVector>(fppexp(X[j], LAM12[j], cut[j], 0.0, 0.0))[0]);

    }

    gauss_quad = fgau(20, 0, C0);
    u = gauss_quad(_,0);
    w = gauss_quad(_,1);
    tmp7 += -log(sum(w*dexp(u, lam01)*as<NumericVector>(fppexp(u, LAM03[0], cut[0], 0.0, 0.0))*as<NumericVector>(fppexp(C0, LAM12[0], cut[0],  0.0, 0.0))[0]/as<NumericVector>(fppexp(u, LAM12[0], cut[0], 0.0, 0.0))));

    rr += tmp1 + tmp2 + tmp3 + tmp5 + tmp7 + tmp8;
  }

  return(rr);
}



//[[Rcpp::export()]]
double loglikFD2_pch(NumericVector par, List Y_F, List X_F,  NumericMatrix Y_proband, NumericMatrix X_proband,
                 NumericVector Age, NumericVector Cal,   NumericVector cut_F, DataFrame lam03,
                 Function fgau, Function combn){

  double theta = exp(par[1]);
  double rho = par[0];
  double newrho = rho/(1+rho);
  NumericVector lam01 = exp(par[seq(2, par.size()-1)]);

  int i, j, j1, j2, k, l, nr, pair_nr;
  int nf = Y_F.size();
  NumericVector tmp(2);
  int ncal = Cal.size();
  int nage = Age.size();
  double C0;
  double tmp11 = 0;
  double tmp1=0, tmp2=0, tmp3=0, tmp5=0, tmp6 = 0, tmp7=0, tmp8=0;
  double rr=0;

  NumericVector pC0 = X_proband(_, 1);
  IntegerVector lam03_Y = lam03["Year.f"];
  IntegerVector lam03_A = lam03["Age.f"];
  NumericVector lam03_rate = lam03["rate"];
  NumericVector lam12_rate = theta*as<NumericVector>(lam03["rate"]);

  NumericMatrix gauss_quad(20, 2);
  NumericVector u(20), u1(20), u2(20);
  NumericVector w(20), w1(20), w2(20);
  IntegerVector sq(20);
  IntegerVector rep1(20*20);
  IntegerVector rep2(20*20);
  NumericVector uu1(20*20);
  NumericVector uu2(20*20);
  NumericVector ww1(20*20);
  NumericVector ww2(20*20);

//#pragma omp parallel for private(i, j, l, k, gauss_quad, u, u1, u2, w, w1, w2, sq, rep1, rep2, uu1, uu2, ww1, ww2, tmp1, tmp2, tmp3, tmp5, tmp7, tmp8) num_threads(4) reduction(+: rr)
  for(i=0; i<nf; i++){
    tmp1 = tmp2 = tmp3 = tmp5 = tmp6 = tmp7 = tmp8 = 0;

    NumericMatrix data_Y = as<NumericMatrix>(Y_F[i]);
    DataFrame data_X = as<DataFrame>(X_F[i]);
    nr = data_Y.nrow();


    IntegerVector proband(nr);
    NumericVector B(nr);
    NumericVector X(nr);
    NumericVector Y(nr);
    NumericVector exam_age(nr);
    IntegerVector del1(nr);
    IntegerVector del2(nr);
    int del3;

    B = data_X["B"];
    X = data_Y(_,0);
    Y = data_X["exam.age"];
    Y[0] = Y_proband(i,1);
    exam_age = data_X["exam.age"];
    del2 = data_Y(_,1);
    del3 = Y_proband(i,2);

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

    if(nr==2){

      NumericVector S0 = vppc(X, lam01, cut_F, 0.0, 0.0);
      if(sum(del2)==2) tmp1 += dClayton(S0, rho, 1.0);
      else tmp1 += log(hf(S0, del2, 1, 2, rho));

      IntegerVector tmpdel = IntegerVector::create(1,0);

      tmp[0] = X[0];
      tmp[1] = exam_age[1];

      gauss_quad = fgau(20, 0, exam_age[1]);
      u = gauss_quad(_,0);
      w = gauss_quad(_,1);

      double px = ppc(X[0], lam01, cut_F, 0.0, 0.0);
      NumericVector vpx = rep(px, 20);

      tmp8 +=  -log(hf(vppc(tmp, lam01, cut_F, 0.0, 0.0), tmpdel, 1, 2, rho)*ppc(exam_age[1], LAM03[1], cut[1], 0.0, 0.0)
                      + sum(w*vdClayton(vpx, vppc(u, lam01, cut_F, 0.0, 0.0), rho, 0.0)*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03[1], cut[1], 0.0, 0.0)*ppc(exam_age[1], LAM12[1], cut[1], 0.0, 0.0)/vppc(u, LAM12[1], cut[1], 0.0, 0.0)));


    }else{
      pair_nr = (nr-1)*(nr-2)/2;
      NumericMatrix comb(2, pair_nr);
      comb = as<NumericMatrix>(combn(nr-1, 2));

      NumericMatrix Sf(nr-1, 2);
      NumericMatrix SAf(nr-1, 2);
      NumericVector Sc(nr-1);
      NumericVector SAc(nr-1);
      NumericMatrix S(pair_nr, 2);
      NumericMatrix SA(pair_nr, 2);
      IntegerMatrix fam_del2(pair_nr, 2);



      for(j=1; j<nr; j++){
        Sf(j-1, 0)  = SAf(j-1, 0) = ppc(X[0], lam01, cut_F, 0.0, 0.0);
        Sf(j-1, 1)  = ppc(X[j], lam01, cut_F, 0.0, 0.0);
        SAf(j-1, 1)  = ppc(exam_age[j], lam01, cut_F, 0.0, 0.0);
        Sc[j-1]  = h1(ppc(X[0], lam01, cut_F, 0.0, 0.0), ppc(X[j], lam01, cut_F, 0.0, 0.0), rho);
        SAc[j-1]  = h1(ppc(X[0], lam01, cut_F, 0.0, 0.0), ppc(exam_age[j], lam01, cut_F, 0.0, 0.0), rho);
      }

      for(j=0; j<pair_nr; j++){
        for(k=0; k<2; k++){
          S(j,k) =  Sc[comb(k,j)-1];
          SA(j,k) =  SAc[comb(k,j)-1];
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



      double px = ppc(X[0], lam01, cut_F, 0.0, 0.0);
      NumericVector vpx = rep(px, 20);
      NumericVector vpx2 = rep(px, 20*20);


      for(j=0; j<pair_nr; j++){
        j1 = comb(0,j);
        j2 = comb(1,j);
        NumericVector cut1 = cut[j1];
        NumericVector cut2 = cut[j2];

        gauss_quad = fgau(20, 0, exam_age[j1]);
        u1 = gauss_quad(_,0);
        w1 = gauss_quad(_,1);

        gauss_quad = fgau(20, 0, exam_age[j2]);
        u2 = gauss_quad(_,0);
        w2 = gauss_quad(_,1);
        sq = seq_len(20)-1;
        rep1 = rep_each(sq, 20);
        rep2 = rep(sq, 20);
        uu1 = u1[rep1];
        uu2 = u2[rep2];
        ww1 = w1[rep1];
        ww2 = w2[rep2];

        tmp8 -= log(pClayton(SA(j,_), newrho)*ppc(exam_age[j1], LAM03[j1], cut[j1], 0.0, 0.0)*ppc(exam_age[j2], LAM03[j2], cut[j2], 0.0, 0.0)
                      + sum(w2*vh1(vh1(vpx, vppc(u2, lam01, cut_F, 0.0, 0.0), rho), rep(SA(j,0), 20), newrho)
                              *vdClayton(vpx, vppc(u2, lam01, cut_F, 0.0, 0.0), rho, 0.0)*vdpc(u2, lam01, cut_F, 0.0)*vppc(u2, LAM03[j2], cut2, 0.0, 0.0)
                              /vppc(u2, LAM12[j2], cut2, 0.0, 0.0))*ppc(exam_age[j2],LAM12[j2],  cut2, 0.0, 0.0)*ppc(exam_age[j1], LAM03[j1], cut1, 0.0, 0.0)
                              + sum(w1*vh1(vh1(vpx, vppc(u1, lam01, cut_F, 0.0, 0.0), rho), rep(SA(j,1), 20), newrho)
                              *vdClayton(vpx, vppc(u1, lam01, cut_F, 0.0, 0.0), rho, 0.0)*vdpc(u1, lam01, cut_F, 0.0)*vppc(u1, LAM03[j1], cut1, 0.0, 0.0)
                              /vppc(u1, LAM12[j1], cut1, 0.0, 0.0))*ppc(exam_age[j1],  LAM12[j1], cut1, 0.0, 0.0)*ppc(exam_age[j2], LAM03[j2], cut2, 0.0, 0.0)
                              +sum(ww1*ww2*vdClayton(vh1(vpx2, vppc(uu1, lam01, cut_F, 0.0, 0.0), rho) , vh1(vpx2, vppc(uu2, lam01, cut_F, 0.0, 0.0), rho),  newrho, 0.0)
                              *vdClayton(vpx2, vppc(uu1, lam01, cut_F, 0.0, 0.0), rho, 0.0)*vdClayton(vpx2, vppc(uu2, lam01, cut_F, 0.0, 0.0), rho, 0.0)
                              *vdpc(uu1, lam01, cut_F, 0.0)*vdpc(uu2, lam01, cut_F, 0.0)*vppc(uu1, LAM03[j1], cut1, 0.0, 0.0)*ppc(exam_age[j1], LAM12[j1], cut1, 0.0, 0.0)
                              /vppc(uu1, LAM12[j1], cut1, 0.0, 0.0)*vppc(uu2, LAM03[j2],  cut2, 0.0, 0.0)
                              *ppc(exam_age[j2], LAM12[j2], cut2, 0.0, 0.0)/vppc(uu2, LAM12[j2], cut2, 0.0, 0.0)))/(nr-2);

      }

    }

    for(j=0; j<nr; j++){
      tmp2 += del2[j]*dpc(X[j], lam01, cut_F, 1.0);
      tmp3 += ppc(X[j], LAM03[j], cut[j], 0.0, 1.0);
      tmp5 += del2[j]*log(ppc(Y[j], LAM12[j], cut[j], 0.0, 0.0)/ppc(X[j], LAM12[j], cut[j], 0.0, 0.0));
    }
    tmp6 += del3*log(dpc(Y[0], LAM12[0], cut[0], 0.0)/ppc(Y[0], LAM12[0], cut[0], 0.0, 0.0));
    gauss_quad = fgau(20, 0, C0);
    u = gauss_quad(_,0);
    w = gauss_quad(_,1);
    tmp7 += -log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03[0], cut[0], 0.0, 0.0)*ppc(C0, LAM12[0], cut[0],  0.0, 0.0)/vppc(u, LAM12[0], cut[0], 0.0, 0.0)));
    rr += tmp1 + tmp2 + tmp3 + tmp5 + tmp6 + tmp7 + tmp8;
  }

  return(rr);
}



//[[Rcpp::export()]]
double loglikFD2_pch_gene(NumericVector par, List Y_F, List X_F,  NumericMatrix Y_proband, NumericMatrix X_proband,
                     NumericVector Age, NumericVector Cal,   NumericVector cut_F, DataFrame lam03,
                     Function fgau, Function combn){

  double theta = exp(par[1]);
  double rho = par[0];
  double newrho = rho/(1+rho);
  double alpha = exp(par[2]);
  double p = exp(par[3]);
  NumericVector lam01 = exp(par[seq(4, par.size()-1)]);

  int i, j, j1, j2, k, l, nr, pair_nr;
  int nf = Y_F.size();
  NumericVector tmp(2);
  int ncal = Cal.size();
  int nage = Age.size();
  double C0;
  double tmp11 = 0;
  double tmp1=0, tmp2=0, tmp3=0, tmp5=0, tmp6=0, tmp7=0, tmp8=0, tmp9=0;
  double rr=0;
  double pg0;

  NumericVector pC0 = X_proband(_, 1);
  IntegerVector lam03_Y = lam03["Year.f"];
  IntegerVector lam03_A = lam03["Age.f"];
  NumericVector lam03_rate = lam03["rate"];
  NumericVector lam12_rate = theta*as<NumericVector>(lam03["rate"]);

  NumericMatrix gauss_quad(20, 2);
  NumericVector u(20), u1(20), u2(20);
  NumericVector w(20), w1(20), w2(20);
  IntegerVector sq(20);
  IntegerVector rep1(20*20);
  IntegerVector rep2(20*20);
  NumericVector uu1(20*20);
  NumericVector uu2(20*20);
  NumericVector ww1(20*20);
  NumericVector ww2(20*20);

  //#pragma omp parallel for private(i, j, l, k, gauss_quad, u, u1, u2, w, w1, w2, sq, rep1, rep2, uu1, uu2, ww1, ww2, tmp1, tmp2, tmp3, tmp5, tmp7, tmp8) num_threads(4) reduction(+: rr)
  for(i=0; i<nf; i++){
    Rcout<<i<<"\n";
    tmp1 = tmp2 = tmp3 = tmp5 = tmp6 = tmp7 = tmp8 = tmp9 = 0;

    NumericMatrix data_Y = as<NumericMatrix>(Y_F[i]);
    DataFrame data_X = as<DataFrame>(X_F[i]);
    nr = data_Y.nrow();


    IntegerVector proband(nr);
    NumericVector B(nr);
    NumericVector X(nr);
    NumericVector Y(nr);
    NumericVector IG(nr);
    NumericVector rid(nr);
    NumericVector exam_age(nr);
    IntegerVector del1(nr);
    IntegerVector del2(nr);
    int del3;

    B = data_X["B"];
    X = data_Y(_,0);
    Y = data_X["exam.age"];
    Y[0] = Y_proband(i,1);
    exam_age = data_X["exam.age"];
    del2 = data_Y(_,1);
    del3 = Y_proband(i,2);

    IG = data_X["G"];
    rid = data_X["rel.id"];
    if(IG[0]==0){
      pg0 = pow((1-p),2);
    }else{
      pg0 = 1-pow((1-p),2);
    }

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

    if(nr==2){
      Rcout<<"two"<<"\n";

      NumericVector S0(2);
      S0[0] = ppc(X[0], lam01*exp(alpha*IG[0]), cut_F, 0.0, 0.0);
      S0[1] = ppc(X[1], lam01*exp(alpha*IG[1]), cut_F, 0.0, 0.0);
      if(sum(del2)==2) tmp1 += dClayton(S0, rho, 1.0);
      else tmp1 += log(hf(S0, del2, 1, 2, rho));
      tmp1 +=  log(pG0(rid, IG, p) / pg0);

      IntegerVector tmpdel = IntegerVector::create(1,0);

      tmp[0] = X[0];
      tmp[1] = exam_age[1];

      gauss_quad = fgau(20, 0, exam_age[1]);
      u = gauss_quad(_,0);
      w = gauss_quad(_,1);

      double px = ppc(X[0], lam01*exp(alpha*IG[0]), cut_F, 0.0, 0.0);
      NumericVector vpx = rep(px, 20);
      NumericVector tmpS1(2), tmpS2(2);
      tmpS1[0] = tmpS2[0] = px;
      tmpS1[1] = ppc(tmp[1], lam01*exp(alpha*0), cut_F, 0.0, 0.0);
      tmpS2[1] = ppc(tmp[1], lam01*exp(alpha*1), cut_F, 0.0, 0.0);
      NumericVector tmpIG1 = NumericVector::create(IG[0], 0);
      NumericVector tmpIG2 = NumericVector::create(IG[0], 1);


      tmp8 +=  -log(hf(tmpS1, tmpdel, 1, 2, rho)*ppc(exam_age[1], LAM03[1], cut[1], 0.0, 0.0)*pG0(rid, tmpIG1, p)/pg0
                  + hf(tmpS2, tmpdel, 1, 2, rho)*ppc(exam_age[1], LAM03[1], cut[1], 0.0, 0.0) *pG0(rid, tmpIG2, p)/pg0
                  + sum(w*vdClayton(vpx, vppc(u, lam01*exp(alpha*0), cut_F, 0.0, 0.0), rho, 0.0)*vdpc(u, lam01*exp(alpha*0), cut_F, 0.0)*vppc(u, LAM03[1], cut[1], 0.0, 0.0)*ppc(exam_age[1], LAM12[1], cut[1], 0.0, 0.0)/vppc(u, LAM12[1], cut[1], 0.0, 0.0))*pG0(rid, tmpIG1, p)/pg0
                  + sum(w*vdClayton(vpx, vppc(u, lam01*exp(alpha*1), cut_F, 0.0, 0.0), rho, 0.0)*vdpc(u, lam01*exp(alpha*1), cut_F, 0.0)*vppc(u, LAM03[1], cut[1], 0.0, 0.0)*ppc(exam_age[1], LAM12[1], cut[1], 0.0, 0.0)/vppc(u, LAM12[1], cut[1], 0.0, 0.0))*pG0(rid, tmpIG2, p)/pg0);


    }else{
      Rcout<<"three"<<"\n";
      pair_nr = (nr-1)*(nr-2)/2;
      NumericMatrix comb(2, pair_nr);
      comb = as<NumericMatrix>(combn(nr-1, 2));

      NumericMatrix Sf(nr-1, 2);
      NumericMatrix SAf(nr-1, 2);
      NumericVector Sc(nr-1);
      NumericVector SAc(nr-1);
      NumericVector SAg0(nr-1);
      NumericVector SAg1(nr-1);
      NumericMatrix S(pair_nr, 2);
      NumericMatrix SA(pair_nr, 2);
      IntegerMatrix fam_del2(pair_nr, 2);
      NumericMatrix fam_rel(pair_nr, 3);
      NumericMatrix fam_IG(pair_nr, 3);



      for(j=1; j<nr; j++){
        Sf(j-1, 0)  = SAf(j-1, 0) = ppc(X[0], lam01*exp(alpha*IG[0]), cut_F, 0.0, 0.0);
        Sf(j-1, 1)  = ppc(X[j], lam01*exp(alpha*IG[j]), cut_F, 0.0, 0.0);
        SAf(j-1, 1)  = ppc(exam_age[j], lam01*exp(alpha*IG[j]), cut_F, 0.0, 0.0);
        Sc[j-1]  = h1(ppc(X[0], lam01*exp(alpha*IG[0]), cut_F, 0.0, 0.0), ppc(X[j], lam01*exp(alpha*IG[j]), cut_F, 0.0, 0.0), rho);
        SAc[j-1]  = h1(ppc(X[0], lam01*exp(alpha*IG[0]), cut_F, 0.0, 0.0), ppc(exam_age[j], lam01*exp(alpha*IG[j]), cut_F, 0.0, 0.0), rho);
        SAg0[j-1]  = h1(ppc(X[0], lam01*exp(alpha*IG[0]), cut_F, 0.0, 0.0), ppc(exam_age[j], lam01*exp(alpha*0), cut_F, 0.0, 0.0), rho);
        SAg1[j-1]  = h1(ppc(X[0], lam01*exp(alpha*IG[0]), cut_F, 0.0, 0.0), ppc(exam_age[j], lam01*exp(alpha*1), cut_F, 0.0, 0.0), rho);

      }

      for(j=0; j<pair_nr; j++){
        fam_rel(j, 0) = rid[0];
        fam_IG(j, 0) = IG[0];
        for(k=0; k<2; k++){
          S(j,k) =  Sc[comb(k,j)-1];
          SA(j,k) =  SAc[comb(k,j)-1];
          fam_del2(j, k) = del2[comb(k,j)];
          fam_rel(j, k+1) = rid[comb(k,j)];
          fam_IG(j, k+1) = IG[comb(k,j)];
        }

        if (sum(fam_del2(j,_))==2) {
          tmp11 = dClayton(S(j,_),  newrho, 1.0) + dClayton(Sf(comb(0,j)-1,_), rho, 1.0) + dClayton(Sf(comb(1,j)-1,_), rho, 1.0);
        }  else if(sum(fam_del2(j,_))==1) {
          tmp11 = log(hf(S(j,_), fam_del2(j,_), 1, 2, newrho)) + fam_del2(j,0)*dClayton(Sf(comb(0,j)-1,_), rho, 1.0) + fam_del2(j,1)*dClayton(Sf(comb(1,j)-1,_), rho, 1.0);
        } else {
          tmp11 =  log(pClayton(S(j,_), newrho));
        }
        tmp1 += tmp11/(nr-2) + log(pG(fam_rel(j,_), fam_IG(j,_), p) / pg0)/(nr-2) ;

      }



      double px = ppc(X[0], lam01*exp(alpha*IG[0]), cut_F, 0.0, 0.0);
      NumericVector vpx = rep(px, 20);
      NumericVector vpx2 = rep(px, 20*20);


      for(j=0; j<pair_nr; j++){
        j1 = comb(0,j);
        j2 = comb(1,j);
        NumericVector cut1 = cut[j1];
        NumericVector cut2 = cut[j2];

        gauss_quad = fgau(20, 0, exam_age[j1]);
        u1 = gauss_quad(_,0);
        w1 = gauss_quad(_,1);

        gauss_quad = fgau(20, 0, exam_age[j2]);
        u2 = gauss_quad(_,0);
        w2 = gauss_quad(_,1);
        sq = seq_len(20)-1;
        rep1 = rep_each(sq, 20);
        rep2 = rep(sq, 20);
        uu1 = u1[rep1];
        uu2 = u2[rep2];
        ww1 = w1[rep1];
        ww2 = w2[rep2];

       NumericVector IG1 = NumericVector::create(IG[0], 0, 0);
       NumericVector IG2 = NumericVector::create(IG[0], 0, 1);
       NumericVector IG3 = NumericVector::create(IG[0], 1, 1);
       NumericVector IG4 = NumericVector::create(IG[0], 1, 1);

       NumericVector SS1 = NumericVector::create(SAg0[j1-1], SAg0[j2-1]);
       NumericVector SS2 = NumericVector::create(SAg0[j1-1], SAg1[j2-1]);
       NumericVector SS3 = NumericVector::create(SAg1[j1-1], SAg0[j2-1]);
       NumericVector SS4 = NumericVector::create(SAg1[j1-1], SAg1[j2-1]);

       NumericVector rel = fam_rel(j,_);


      tmp8 -= log(ff1(j1, j2, vpx,vpx2, alpha, lam01, newrho, rho, exam_age,cut_F, LAM03,  LAM12, cut1, cut2, IG1,  SS1, rel, pg0, p, w1, w2, u1, u2, ww1, ww2, uu1, uu2)
                    + ff1(j1, j2, vpx,vpx2, alpha, lam01, newrho, rho, exam_age,cut_F, LAM03, LAM12, cut1, cut2, IG2,  SS2, rel, pg0, p, w1, w2, u1, u2, ww1, ww2, uu1, uu2)
                    + ff1(j1, j2, vpx,vpx2, alpha, lam01, newrho, rho, exam_age,cut_F, LAM03, LAM12, cut1, cut2, IG3,  SS3, rel, pg0, p, w1, w2, u1, u2, ww1, ww2, uu1, uu2)
                    + ff1(j1, j2, vpx,vpx2, alpha, lam01, newrho, rho, exam_age,cut_F, LAM03, LAM12, cut1, cut2, IG4,  SS4, rel, pg0, p, w1, w2, u1, u2, ww1, ww2, uu1, uu2))/(nr-2);

      }

    }

    for(j=0; j<nr; j++){
      tmp2 += del2[j]*dpc(X[j], lam01*exp(alpha*IG[j]), cut_F, 1.0);
      tmp3 += ppc(X[j], LAM03[j], cut[j], 0.0, 1.0);
      tmp5 += del2[j]*log(ppc(Y[j], LAM12[j], cut[j], 0.0, 0.0)/ppc(X[j], LAM12[j], cut[j], 0.0, 0.0));
    }

    tmp6 += del3*log(dpc(Y[0], LAM12[0], cut[0], 0.0)/ppc(Y[0], LAM12[0], cut[0], 0.0, 0.0));

    gauss_quad = fgau(20, 0, C0);
    u = gauss_quad(_,0);
    w = gauss_quad(_,1);
    tmp7 += -log(sum(w*vdpc(u, lam01, cut_F, 0.0)*vppc(u, LAM03[0], cut[0], 0.0, 0.0)*ppc(C0, LAM12[0], cut[0],  0.0, 0.0)/vppc(u, LAM12[0], cut[0], 0.0, 0.0))*pow(1-p,2)
                 + sum(w*vdpc(u, lam01*exp(alpha), cut_F, 0.0)*vppc(u, LAM03[0], cut[0], 0.0, 0.0)*ppc(C0, LAM12[0], cut[0],  0.0, 0.0)/vppc(u, LAM12[0], cut[0], 0.0, 0.0))*(1-pow(1-p,2)));
    tmp9 += log(pg0);

    rr += tmp1 + tmp2 + tmp3 + tmp5 + tmp6 + tmp7 + tmp8;
  }

  return(rr);
}


