// include guard
#ifndef __poly_h_INCLUDED__   // if poly.h hasn't been included yet...
#define __poly_h_INCLUDED__   // #define this so the compiler knows it has been included

// see https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/stat_tut/weg/error_eg.html
// and https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/pol_tutorial/changing_policy_defaults.html
// return NAN for invalid inputs
#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
#include <boost/math/distributions.hpp>       // gamma, normal, lognormal distributions

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// next two lines must be in order (RcppNumerical precedes secr.h)
#include <RcppNumerical.h>

#include "secr.h"

using namespace Rcpp;

rpoint getxycpp(
        const double l, 
        const std::vector<double> &cumd, 
        const RcppParallel::RMatrix<double> &line, 
        const int n1, 
        const int n2);
    
//--------------------------
// original using Rdqags
//--------------------------

double hintegral1Dcpp (int fn, const NumericVector par);

double integral1Dcpp (
        int fn, 
        int m, 
        int c, 
        NumericVector gsbval, 
        int cc, 
        NumericMatrix traps,
        NumericMatrix mask, 
        int n1, 
        int n2, 
        int kk, 
        int mm);

// double integral2Dcpp  (
//         const int &fn,
//         const int &m,
//         const int &c,
//         const RcppParallel::RMatrix<double> &gsbval,
//         const RcppParallel::RMatrix<double> &poly,
//         const RcppParallel::RMatrix<double> &mask,
//         const int &n1,
//         const int &n2,
//         double ex[]);

//-----------------------------------------------------
// alternative 2-D using RcppNumerical Numer::integrate
//-----------------------------------------------------

double hintegral1DNRcpp (
        const int fn, 
        const std::vector<double> &gsb);

double hintegral2DNRcpp (
        const int fn, 
        const std::vector<double> &gsb); 

double integral1DNRcpp (
        const int fn, 
        const int m, 
        const int c, 
        const RcppParallel::RMatrix<double> &gsbval, 
        const RcppParallel::RMatrix<double> &traps,
        const RcppParallel::RMatrix<double> &mask, 
        const int n1, 
        const int n2);

double integral2DNRcpp  (
        const int &fn,
        const int &m,
        const int &c,
        const RcppParallel::RMatrix<double> &gsbval,
        const RcppParallel::RMatrix<double> &poly,
        const RcppParallel::RMatrix<double> &mask,
        const int &n1,
        const int &n2,
        const bool &convex);

#endif

