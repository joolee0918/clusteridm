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
  res  = (1+rho)*pow(u1*u2, -rho-1)*pow((pow(u1, -rho) + pow(u2, -rho)-1), -1/rho-2);
  if(logf == 1.0) res = log(res);
  return(res);
}


NumericVector vdClayton(NumericVector u1, NumericVector u2, double rho, bool logf){
  NumericVector res;
  res  = (1+rho)*pow(u1*u2, -rho-1)*pow((pow(u1, -rho) + pow(u2, -rho)-1), -1/rho-2);
  if(logf == 1.0) res = log(res);
  return(res);
}


double pClayton(NumericVector u, double rho){
  double res;
  double u1 = u[0];
  double u2 = u[1];
  res = pow(pow(u1,-rho) + pow(u2,-rho)-1, -1/rho);
  return(res);
}



double hf(NumericVector u, IntegerVector del, int k, int mi, double rho){

  int j;

  double term1, term2, tmp;
  double res;

  NumericVector term(mi);

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
  return(res);
}


double h1(double u1, double u2, double rho){
  double res;
  res = pow(u1, -rho-1)*pow((pow(u2, -rho) + pow(u1, -rho) -1), -1/rho-1);
  return(res);
}

NumericVector vh1(NumericVector u1, NumericVector u2, double rho){
  NumericVector res;
  res = pow(u1, -rho-1)*pow((pow(u2, -rho) + pow(u1, -rho) -1), -1/rho-1);
  return(res);
}

NumericVector cumsum3(NumericVector x) {
  return cumsum(x); // compute + return result
}


NumericVector cumprod3(NumericVector x) {
  return cumprod(x); // compute + return result
}

