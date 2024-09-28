#include "secr.h"

using namespace std;
using namespace Rcpp;
using namespace RcppParallel; 


//===============================================================================

struct Hckm : public Worker {
  
  // input data
  const int detectfn;
  const RMatrix<double> gsbval;
  const RMatrix<double> dist2;
  const RVector<double> miscparm;
  
  // output vector to write to
  RVector<double> gk;
  RVector<double> hk;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  Hckm(const int detectfn,
       const NumericMatrix gsbval, 
       const NumericMatrix dist2, 
       const NumericVector miscparm, 
       NumericVector gk,
       NumericVector hk)
    : detectfn(detectfn), gsbval(gsbval), 
      dist2(dist2), miscparm(miscparm), gk(gk), hk(hk) {
    
  }
  
  // hazard (fn 14:19) or g ()fn 0:11) (distance lookup)
  double zLcpp (
      const int c, 
      const int k, 
      const int m)
  {
      double r = 0.0;
      double r2;
      double temp;
      r2 = dist2(k,m);
      if ((detectfn == 0) || (detectfn == 14)) {        // halfnormal or hazard halfnormal
          return (gsbval(c,0) * std::exp(-r2 / 2 / gsbval(c,1)/ gsbval(c,1)));
      }
      else if (detectfn == 3) {                         // compound halfnormal
          temp = gsbval(c,0) * std::exp(- r2  / 2 / gsbval(c,1) / gsbval(c,1));
          if (round(gsbval(c,2)) > 1) temp = 1 - pow(1 - temp, gsbval(c,2));
          return (temp);
      }
      else {
          r = std::sqrt(r2);
          if ((detectfn == 1) || (detectfn == 15)) {       // hazard rate or hazard hazard rate
              return (gsbval(c,0) * ( 1 - std::exp(- pow(r /gsbval(c,1), - gsbval(c,2)))));
          }
          else if ((detectfn == 2) || (detectfn == 16)) {  // exponential or hazard exponential
              return (gsbval(c,0) * std::exp(-r / gsbval(c,1)));
          }
          else if (detectfn == 4) {                        // uniform
              if (r<gsbval(c,1)) return (gsbval(c,0));
              else return (0);
          }
          else if (detectfn == 5) {                        // w exponential
              if (r<gsbval(c,2)) return (gsbval(c,0));
              else return (gsbval(c,0) * std::exp(-(r-gsbval(c,2)) / gsbval(c,1)));
          }
          else if ((detectfn == 6) || (detectfn == 17)) {  // annular normal or hazard annular normal
              return (gsbval(c,0) * std::exp(-(r-gsbval(c,2))*(r-gsbval(c,2)) / 2 /
                      gsbval(c,1) / gsbval(c,1)));
          }
          else if (detectfn == 7) {                        // cumulative lognormal
              double CV2, meanlog, sdlog;
              CV2 = gsbval(c,2)*gsbval(c,2)/gsbval(c,1)/gsbval(c,1);
              meanlog = log(gsbval(c,1)) - log(1 + CV2)/2;
              sdlog = std::sqrt(log(1 + CV2));
              boost::math::lognormal_distribution<> ln(meanlog,sdlog);
              return (gsbval(c,0) * boost::math::cdf(complement(ln,r))); 
          }
          else if ((detectfn == 8) || (detectfn == 18)) {  // cumulative gamma or hazard cumulative gamma
              boost::math::gamma_distribution<> gam(gsbval(c,2), gsbval(c,1)/gsbval(c,2));
              return (gsbval(c,0) * boost::math::cdf(complement(gam,r))); 
          }
          else if (detectfn == 9) {                       // binary signal strength
              // b0 = gsbval(c,0)
              // b1 = gsbval(c,1)
              temp = -(gsbval(c,0) + gsbval(c,1) * r);
              boost::math::normal_distribution<> n;
              return (boost::math::cdf(complement(n,temp)));    // upper 
          }
          else if (detectfn == 10 || detectfn == 11) {   // signal strength, signal strength spherical
              double mu, gam;
              // beta0 = gsbval(c,0)
              // beta1 = gsbval(c,1)
              // sdS = gsbval(c,2)
              // cut = miscparm(0)
              if (detectfn == 10) 
                  mu = gsbval(c,0) + gsbval(c,1) * r;
              else
                  mu = gsbval(c,0) + gsbval(c,1) * (r-1) - 10 * log(r*r) / M_LN10;
              gam = (miscparm[0] - mu) / gsbval(c,2);
              // return (R::pnorm(gam,0,1,0,0));    
              boost::math::normal_distribution<> n;
              return (boost::math::cdf(complement(n,gam)));
          }
          else if (detectfn == 19) {  // hazard variable power
              return (gsbval(c,0) * std::exp(- pow(r /gsbval(c,1), gsbval(c,2))));
          }
          else (Rcpp::stop("unknown or invalid detection function"));
      }
  }
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
      int cc = gsbval.nrow();
      int kk = dist2.nrow();
      for (std::size_t m = begin; m < end; m++) {
          for (int k=0; k < kk; k++) {
              for (int c=0; c < cc; c++) {
                  int gi = cc * (kk * m + k) + c;  
                  if (detectfn<12) {  // g output from zLcpp
                      gk[gi] = zLcpp(c, k, m);
                      hk[gi] = -log(1-gk[gi]);
                  }
                  else {              // hazard output from zLcpp
                      hk[gi] = zLcpp(c, k, m);
                      gk[gi] = 1 - std::exp(-hk[gi]);
                  }
              }
          }
      }
  }
};

// [[Rcpp::export]]
List makegkPointcpp (const int detectfn, 
                     const int grain,
                     const int ncores,
                     const NumericMatrix& gsbval, 
                     const NumericMatrix& dist2,
                     const NumericVector& miscparm
) 
{
    NumericVector hk(gsbval.nrow() * dist2.size()); 
    NumericVector gk(gsbval.nrow() * dist2.size()); 
    
    Hckm hckm (detectfn, gsbval, dist2, miscparm, gk, hk);
    
    if (ncores>1) {
        parallelFor(0, dist2.ncol(), hckm, grain, ncores);
    }
    else {
        // for debugging avoid multithreading to allow R calls
        hckm.operator()(0,dist2.ncol());    
    }
    return List::create(Named("gk") = gk, Named("hk") = hk);
}
//==============================================================================

// [[Rcpp::export]]
List cappedgkhkcpp (
        const int cc, 
        const int nk, 
        const double area,
        const NumericVector &D,
        NumericVector &gk, 
        NumericVector &hk)  
{ 
    // adjust gk, hk for competition among animals for detectors
    int gi, c,i,k,m,mm;
    double H;        // detector-specific total hazard for capped detectors
    double pHH;      // probability of at least one detection / H
    mm = D.size();
    for (c=0; c<cc; c++) {
        for (k=0; k<nk; k++) {
            H = 0;
            for (m=0; m<mm; m++) {
                gi = i3(c,k,m, cc,nk);
                H += hk[gi] * D[m] * area;
            }
            pHH = (1 - std::exp(-H))/H;
            for (m=0; m<mm; m++) {
                gi = i3(c,k,m, cc,nk);
                gk[gi] = pHH * hk[gi];
            }
        }
    }
    // outside k loop to avoid messing up hk
    for (i=0; i<hk.size(); i++) hk[i] = -log(1-gk[i]);
    return List::create(Named("gk") = gk, Named("hk") = hk);
}
