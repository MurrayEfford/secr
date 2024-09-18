#include "poly.h"
// #include <algorithm>  

using namespace Rcpp;
using namespace RcppParallel;

//===============================================================================

// 2019-11-25
// FTHL.fit <- secr.fit(hornedlizardCH, buffer = 80, trace = TRUE)
//     Warning: stack imbalance in '<-', 2 then -4

// Error in .localstuff$iter <- .localstuff$iter + 1 : 
//     R_Reprotect: only 1 protected item, can't reprotect index -1

struct Hckmpoly : public Worker {
  
  // input data
  const int             detectfn;
  const int             dim;
  const int             grain;
  const bool            convex;
  const RMatrix<double> gsbval;
  const RVector<int>    cumk;
  const RMatrix<double> traps;
  const RMatrix<double> mask;
  
  // output vector to write to
  RVector<double>       H;
  RVector<double>       gk;
  RVector<double>       hk;
  
  int cc, kk, nk, mm, npar;
  // ex is used only by Rdqags in integral2Dcpp
  // double *ex;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  Hckmpoly(const int detectfn,
      const int dim,
      const int grain,
      const bool convex,
      const NumericMatrix gsbval, 
      const IntegerVector cumk, 
      const NumericMatrix traps, 
      const NumericMatrix mask, 
      NumericVector H,
      NumericVector gk,
      NumericVector hk)
      : detectfn(detectfn), dim(dim), grain(grain), convex(convex), gsbval(gsbval), 
          cumk(cumk), traps(traps), mask(mask), H(H), gk(gk), hk(hk) {
      
    cc = gsbval.nrow();
    kk = traps.nrow();
    mm = mask.nrow();
    nk = cumk.size()-1;
    npar = gsbval.ncol();
  }

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) 
  {
      for (int c=0; c<cc; c++) {
          std::vector<double> gsb(4,0);
          for (int i=0; i<npar; i++) gsb[i] = gsbval(c, i);
          
          if (dim==1)
              H[c] = hintegral1DNRcpp(detectfn, gsb);    // unbounded integrated hazard from radial function
          else
              H[c] = hintegral2DNRcpp(detectfn, gsb);    // unbounded integrated hazard from radial function
          for (int k=0; k<nk; k++) {                    // over parts 
              for (std::size_t m = begin; m < end; m++) {
                  int gi = i3(c, k, m, cc, nk);
                  // strictly use hazard form : expected detections of animals at m 
                  // gsb[0] only makes sense here if fn is HHN, HHR, HEX, HAN HCG 
                  
                  // subset traps rows cumk[k] : (cumk[k+1]-1)
                  int n1 = cumk[k];
                  int n2 = cumk[k+1]-1;
                  
                  if (dim==1) {
                      hk[gi] = gsb[0] * integral1DNRcpp (detectfn, m, c, gsbval, traps, mask, n1, n2) / H[c];
                  }
                  else {
                      hk[gi] = gsb[0] * integral2DNRcpp (detectfn, m, c, gsbval, traps, mask, n1, n2, convex) / H[c];
                  }
                  
                  gk[gi] = 1 - exp(-hk[gi]);
              }
          }
      }
  }
};

// [[Rcpp::export]]
List makegkPolygoncpp (const int detectfn, 
                       const int dim,
                       const bool convex,
                       const int grain,
                       const int ncores,
                       const NumericMatrix& gsbval, 
                       const IntegerVector& cumk,
                       const NumericMatrix& traps,
                       const NumericMatrix& mask
) 
{
  NumericVector H(gsbval.nrow()); 
  NumericVector gk(gsbval.nrow() * (cumk.size()-1) * mask.nrow()); 
  NumericVector hk(gsbval.nrow() * (cumk.size()-1) * mask.nrow()); 
  
  Hckmpoly hckm (detectfn, dim, grain, convex, gsbval, cumk, traps, mask, H, gk, hk);
  
  if (ncores>1) {
      parallelFor(0, mask.nrow(), hckm, grain, ncores);
  }
  else {
      hckm.operator()(0,mask.nrow());    // for debugging avoid multithreading to allow R calls
  }
  return List::create(Named("H") = H, Named("gk") = gk, Named("hk") = hk);
}
//==============================================================================
