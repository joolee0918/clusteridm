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

using namespace Rcpp;


//[[Rcpp::export()]]
double loglikFD2(NumericVector par, List outdata_F, NumericVector outdata_proband,
               NumericVector Age, NumericVector Cal,  DataFrame lam03,
               Function fgau, Function fdpexp, Function fppexp, Function combn, Function ppch){

  double lam01 = exp(par[0]);
  double theta = exp(par[1]);
  double rho = exp(par[2]);
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

    proband = data["proband"];
    B = data["B"];
    X = data["X"];
    Y = data["Y"];
    exam_age = data["exam.age"];
    del1 = data["del1"];
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


