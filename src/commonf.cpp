#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>


using namespace Rcpp;


IntegerVector findInterval(NumericVector x, NumericVector breaks) {
  IntegerVector out(x.size());
  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;
  for(it = x.begin(), out_it = out.begin(); it != x.end();
  ++it, ++out_it) {
    pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
    *out_it = std::distance(breaks.begin(), pos);
  }
  return out;
}



double dClayton(NumericVector u, double rho, bool logf){
  double res;
  double u1 = u[0];
  double u2 = u[1];
  if(rho==0) {
    res = 1;
  }else{
  res  = (1+rho)*pow(u1*u2, -rho-1)*pow((pow(u1, -rho) + pow(u2, -rho)-1), -1/rho-2);
  }
  if(logf == 1.0) res = log(res);
  return(res);
}


NumericVector vdClayton(NumericVector u1, NumericVector u2, double rho, bool logf){
  NumericVector res;
  if(rho == 0) {
    res.fill(1);
  } else{
  res  = (1+rho)*pow(u1*u2, -rho-1)*pow((pow(u1, -rho) + pow(u2, -rho)-1), -1/rho-2);
  }
  if(logf == 1.0) res = log(res);
  return(res);
}


double pClayton(NumericVector u, double rho){
  double res;
  double u1 = u[0];
  double u2 = u[1];
  if(rho==0) {
    res = u1*u2;
  } else{
  res = pow(pow(u1,-rho) + pow(u2,-rho)-1, -1/rho);
  }
  return(res);
}



double hf(NumericVector u, IntegerVector del, int k, int mi, double rho){

  int j;

  double term1, term2, tmp;
  double res;


  NumericVector term(mi);

  if(rho == 0){
    res = 1;
  }else{
  for(j=0; j<mi; j++){
    term[j] = pow(u[j], -rho);
  }
  if(k>1) {
    term1=1;
    for(j=1; j<k; j++) {
      tmp = (1+j*rho);
      term1 = term1*tmp;
    }
  } else {
    term1=1;
  }
  term2 = 1;
  for(j=0; j<mi; j++){
    term2 *= pow(u[j], del[j]);
  }

  res = term1*pow(term2, -rho-1)*pow(sum(term)-mi+1, -1/rho-k);
  }
  return(res);
}


double h1(double u1, double u2, double rho){
  double res;
  if(rho == 0){
    res = u2;
  }else{
  res = pow(u1, -rho-1)*pow((pow(u2, -rho) + pow(u1, -rho) -1), -1/rho-1);
  }
  return(res);
}

NumericVector vh1(NumericVector u1, NumericVector u2, double rho){
  NumericVector res;
  if(rho == 0){
    res = u2;
  }else{
  res = pow(u1, -rho-1)*pow((pow(u2, -rho) + pow(u1, -rho) -1), -1/rho-1);
  }
  return(res);
}

NumericVector cumsum3(NumericVector x) {
  return cumsum(x); // compute + return result
}


NumericVector cumprod3(NumericVector x) {
  return cumprod(x); // compute + return result
}


//[[Rcpp::export()]]
NumericVector hpc(NumericVector x, NumericVector levels, NumericVector cuts, int logf)
{

 NumericVector cut = sort_unique(cuts);
 int p = levels.size();
 NumericVector y(x.size());
 cut.push_front(0);
 cut.push_back(R_PosInf);

  y[(cut[0] <= x) & (x < cut[1])] = levels[0];
      if (p > 1.5) {
        for (int i=1; i<p; i++) {
          y[(cut[i] <= x) & (x < cut[i + 1])] = levels[i];
        }
      }
      if (logf)
        y = log(y);
 return(y);
}

//[[Rcpp::export()]]
NumericVector Hpc(NumericVector x,  NumericVector levels, NumericVector cuts, int logf)
{
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  NumericVector y(x.size());
  cut.push_front(0);
  cut.push_back(R_PosInf);
  NumericVector z;
  LogicalVector who = (cut[0] <= x) & (x < cut[0 + 1]);
        if (sum(who)) {
          y[who] = x[who];
          y = y*levels[0];
        }
        double su = levels[0] * cut[1];
        if (p > 1.5) {
          for (int i = 1; i<p; i++) {
            who = (cut[i] <= x) & (x < cut[i + 1]);
            if (sum(who)) {
              NumericVector xwho= x[who];
              NumericVector tmpx = su + levels[i] * (xwho - cut[i]);
              y[who] = tmpx;
            }
            su = su + levels[i] * (cut[i + 1] - cut[i]);
          }
        }
        if (logf)
          y = log(y);
    return(y);
}


//[[Rcpp::export()]]
double ppc(double q, NumericVector levels,  NumericVector cuts, int lower, int logf)
{
  double y;
  if (cuts[0]==0) {
     y = R::pexp(q, 1/levels[0], 0.0, 0.0);
    }else{
  NumericVector qq(1);
  qq[0] = q;
   y = Hpc(qq,  levels, cuts, 0.0)[0];
    if (logf) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}


//[[Rcpp::export()]]
NumericVector vppc(NumericVector q, NumericVector levels,  NumericVector cuts, int lower, int logf)
{
  NumericVector y(1);
  if (cuts[0]==0) {
    y = pexp(q, levels[0], lower, logf);
  }else{
    y = Hpc(q,  levels, cuts, 0.0);
    if (logf) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}


//[[Rcpp::export()]]
double dpc(double x, NumericVector levels,  NumericVector cuts, int logf)
{

  double y;
  if (cuts[0]==0) {
    y = R::dexp(x, 1/levels[0], 0.0);
  }else{
    NumericVector xx(1);
    xx[0] = x;
  y = hpc(xx, levels, cuts,  0.0)[0] * ppc(x, levels, cuts, 0.0, 0.0);
  }
    if (logf)
      y = log(y);

    return(y);
}

//[[Rcpp::export()]]
NumericVector vdpc(NumericVector x, NumericVector levels,  NumericVector cuts, int logf)
{

  NumericVector y(x.size());
  if (cuts[0]==0) {
    y = dexp(x, levels[0]);
  }else{
    y = hpc(x, levels, cuts,  0.0) * vppc(x, levels, cuts, 0.0, 0.0);
  }
  if (logf)
    y = log(y);

  return(y);
}
