#ifndef COMMONF_H
#define COMMONF_H

#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>


using namespace Rcpp;


IntegerVector findInterval(NumericVector x, NumericVector breaks);
double dClayton(NumericVector u, double rho, bool logf);
NumericVector vdClayton(NumericVector u1, NumericVector u2, double rho, bool logf);
double pClayton(NumericVector u, double rho);
double hf(NumericVector u, IntegerVector del, int k, int mi, double rho);
double h1(double u1, double u2, double rho);
NumericVector vh1(NumericVector u1, NumericVector u2, double rho);
NumericVector cumsum3(NumericVector x);
NumericVector cumprod3(NumericVector x);
NumericVector hpc(NumericVector x, NumericVector levels, NumericVector cuts,  int logf);
NumericVector Hpc(NumericVector x, NumericVector levels, NumericVector cuts,  int logf);
double ppc(double q, NumericVector levels, NumericVector cuts, int lower, int logf);
double dpc(double x, NumericVector levels, NumericVector cuts, int logf);
NumericVector vppc(NumericVector q, NumericVector levels, NumericVector cuts, int lower, int logf);
NumericVector vdpc(NumericVector x, NumericVector levels, NumericVector cuts, int logf);
#endif

