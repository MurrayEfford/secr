#include <Rcpp.h>
#include <RcppParallel.h>
#include "secr.h"

//==============================================================================
// 2020-01-26
struct polygonfxi : public Worker {
  // input data
  const int             nc;
  const int             detectfn;
  const int             grain;
  const double          minp;
  const RVector<int>    binomN;   // s
  const RVector<int>    w;        // n x s x k
  const RMatrix<double> xy;
  const RVector<int>    start;    // starting position in xy of detections of animal i,s,k
  const RVector<int>    group;    // n
  const RVector<double> hk;
  const RVector<double> H;
  const RMatrix<double> gsbval;
  const RMatrix<double> pID;
  const RMatrix<double> mask;
  const RMatrix<double> density;  // n x g
  const RVector<int>    PIA;      // 1,n,s,k,1
  const RMatrix<double> Tsk; 
  const RMatrix<double> h;
  const RMatrix<int>    hindex;
  
  const RMatrix<int>    mbool;      // appears cannot use RMatrix<bool>

  // working variables
  int  mm, nk, ss, cc;
  
  // output likelihoods
  RMatrix<double> output;
  
  // Constructor to initialize an instance of Somehistories
  // The RMatrix class can be automatically converted to from the Rcpp matrix type
  polygonfxi(
    const int           nc,
    const int           detectfn,
    const int           grain,
    const double        minp,
    const IntegerVector binomN,
    const IntegerVector w,
    const NumericMatrix xy,
    const IntegerVector start,
    const IntegerVector group,
    const NumericVector hk,
    const NumericVector H,
    const NumericMatrix gsbval,
    const NumericMatrix pID,
    const NumericMatrix mask,
    const NumericMatrix density,
    const IntegerVector PIA,
    const NumericMatrix Tsk, 
    const NumericMatrix h,
    const IntegerMatrix hindex, 
    
    const LogicalMatrix mbool,
    
    NumericMatrix output)
    :
    nc(nc), detectfn(detectfn), grain(grain), minp(minp), 
    binomN(binomN), w(w), xy(xy), start(start), group(group), hk(hk), H(H), gsbval(gsbval), 
    pID(pID), mask(mask), density(density), PIA(PIA), Tsk(Tsk),  h(h), hindex(hindex), mbool(mbool),
    output(output) {
    // now can initialise these derived counts
    mm = mask.nrow();       // number of polygons (detectors)
    nk = Tsk.nrow();        // number of polygons (detectors)
    ss = Tsk.ncol();        // number of occasions
    cc = gsbval.nrow();     // number of parameter combinations
  }
  //==============================================================================
  
  double d2Rcpp (
      const int k,
      const int m,
      const RMatrix<double> &A1,
      const RMatrix<double> &A2)
    // return squared distance between two points given by row k in A1
    // and row m in A2, where A1 and A2 have respectively A1rows and A2rows
  {
    return(
      (A1(k,0) - A2(m,0)) * (A1(k,0) - A2(m,0)) +
        (A1(k,1) - A2(m,1)) * (A1(k,1) - A2(m,1))
    );
  }
  //--------------------------------------------------------------------------
  
  // hazard (fn 14:19) (distance)
  double zcpp (const int k, const int m, const int c, const RMatrix<double> &gsbval, 
               const RMatrix<double> &xy, const RMatrix<double> &mask)
  {
    double r, r2;
    r2 = d2Rcpp(k, m, xy, mask);
    if (detectfn == 14) {  // hazard halfnormal
      return (gsbval(c,0) *  exp(-r2 / 2 / gsbval(c,1) / gsbval(c,1)));    
    }
    else {
      r = std::sqrt(r2);
      if (detectfn == 15) {  // hazard hazard rate
        return (gsbval(c,0) * ( 1 - exp(- pow(r /gsbval(c,1), - gsbval(c,2)))));
      }
      else if (detectfn == 16) {  // hazard exponential
        return (gsbval(c,0) * exp(-r / gsbval(c,1)));
      }
      else if (detectfn == 17) {  // hazard annular normal
        return (gsbval(c,0) * exp(-(r-gsbval(c,2))*(r-gsbval(c,2)) / 
                2 / gsbval(c,1)/ gsbval(c,1)));
      }
      else if (detectfn == 18) {  // hazard cumulative gamma
        //return (gsbval(c,0) * R::pgamma(r,gsbval(c,2),gsbval(c,1)/gsbval(c,2),0,0)); 
        boost::math::gamma_distribution<> gam(gsbval(c,2),gsbval(c,1)/gsbval(c,2));
        return (gsbval(c,0) * boost::math::cdf(complement(gam,r))); 
        
      }
      else if (detectfn == 19) {  // hazard variable power
        return (gsbval(c,0) * exp(- pow(r /gsbval(c,1), gsbval(c,2))));
      }
      else (Rcpp::stop("unknown or invalid detection function"));
    }
  }
  
  void prwpolygonXfxi (const int n, std::vector<double> &pm) {
      // Likelihood component due to capture history n (0 <= n < nc)
      // given that animal's range centre is at m
      // EXCLUSIVE POLYGON DETECTOR
      {
          int s;   // index of occasion  0 <= s < ss  
          int k;   // index of part 0 <= k < nk  
          int c, m, w2, w3, gi;
          double hint, Tski, Htemp;
          bool dead = false;
          for (s = 0; s < ss; s++) {   // over occasions
              w2 = s * nc + n;
              k = w[w2]; 
              dead = k < 0;  
              k = abs(k)-1;         // detector number 0..nk-1; k = -1 if not caught 
              // Not found at any detector on occasion s 
              if (k < 0) {
                  for (m=0; m<mm; m++) {
                      if (mbool(n,m)) {
                          Htemp = h(m, hindex(n,s));
                          //if (Htemp > fuzz)
                              pm[m] *= exp(-Htemp);
                      }
                      else {
                          pm[m] = 0.0;
                      }
                  }
              }
              // detected at detector k on occasion s
              else {
                  w3 = i3(n, s, k, nc, ss);   
                  c = PIA[w3] - 1;
                  if (c >= 0) {    // ignore unused detectors 
                      Tski = Tsk(k,s);
                      for (m=0; m<mm; m++) {
                          if (mbool(n,m)) {
                              gi  = i3(c,k,m,cc,nk);
                              Htemp = h(m, hindex(n,s));
                              pm[m] *=  Tski * (1-exp(-Htemp)) *  hk[gi] / Htemp;
                              // if ((grain==0) && (m==1000)) {
                              //     Rprintf("pm[m]  %10.7e hk[gi] %10.7e gsbval(c,0) %10.7e H[c] %10.7e \n", pm[m], hk[gi], gsbval(c,0), H[c]);
                              // }
                              // for each detection compute pdf(xy) | detected 
                              if (pm[m] > minp) {               // avoid underflow 
                                  // retrieve hint = integral2D(zfn(x) over k)) 
                                  hint = hk[gi] / gsbval(c,0) * H[c];  
                                  pm[m] *= zcpp(start[w3], m, c, gsbval, xy, mask) / hint;
                                  // if ((grain==0) & (m==1000)) Rprintf("n %4d j %4d pm[m]  %10.7e hk[gi] %10.7e gsbval(c,0) %10.7e H[c] %10.7e hint %10.7e with zcpp\n", n, start[i3(n,s,k,nc,ss)], pm[m], hk[gi], gsbval(c,0), H[c], hint);
                              }
                          }
                          else {
                              pm[m] = 0.0;
                          }
                      }
                  }
              }
              if (dead) break;   // out of s loop
          }
      }
  }    
  //==============================================================================
  void prwpolygonfxi (const int n, std::vector<double> &pm) {
      // Likelihood component due to capture history n (0 <= n < nc)
      // given that animal's range centre is at m
      // POLYGON DETECTOR
      {
          int s;   // index of occasion  0 <= s < ss  
          int k;   // index of part 0 <= k < nk  
          int j;   // index of xy record 
          int c, m, w3, gi;
          long count;
          bool dead = false;
          double hint;
          double Tski;
          
          for (s=0; s<ss; s++) {  // over occasions
              if (binomN[s] < 0) Rcpp::stop ("negative binomN < 0 not allowed in C++ fn prwpolygonfxi");
              for (k=0; k<nk; k++) {   // over polygons
                  w3 = i3(n,s,k,nc,ss);
                  count = w[w3];
                  dead = count<0;
                  count = abs(count);
                  c = PIA[w3] - 1;
                  if (c >= 0) {                          // skip if this polygon not used 
                      Tski = Tsk(k,s);
                      for (m=0; m<mm; m++) {
                          if (mbool(n,m)) {
                              gi  = i3(c,k,m,cc,nk);
                              pm[m] *= pski(binomN[s], count, Tski, hk[gi], 1.0);
                              
                              // for each detection, pdf(xy) | detected 
                              if ((pm[m] > minp) && (count>0)) {       // avoid underflow
                                  // retrieve hint = integral2D(zfn(x) over k)) 
                                  hint = hk[gi] / gsbval(c,0) * H[c];  
                                  for (j=start[w3]; j < start[w3]+count; j++) {
                                      pm[m] *= zcpp(j, m, c, gsbval, xy, mask) / hint;
                                  }
                              }
                          }
                          else {
                              pm[m] = 0.0;
                          }
                      }
                  }
              }
              if (dead==1) break;
          }
      }
  }    
  //==============================================================================
  void prwtransectfxi (const int n, std::vector<double> &pm) {
      // Likelihood component due to capture history n (0 <= n < nc)
      // given that animal's range centre is at m
      // TRANSECT DETECTOR
      {
          int s;   // index of occasion  0 <= s < ss  
          int k;   // index of part 0 <= k < nk  
          int j;   // index of xy record 
          int c, m, w3, gi;
          long count;
          bool dead = false;
          double hint;
          double Tski;
          
          for (s=0; s<ss; s++) {  // over occasions
              if (binomN[s] < 0) Rcpp::stop ("negative binomN < 0 not allowed in C++ fn prwitransectfxi");
              for (k=0; k<nk; k++) {   // over transects
                  w3 = i3(n,s,k,nc,ss);
                  count = w[w3];
                  dead = count<0;
                  count = abs(count);
                  c = PIA[w3] - 1;
                  if (c >= 0) {                          // skip if this transect not used 
                      Tski = Tsk(k,s);
                      for (m=0; m<mm; m++) {
                          if (mbool(n,m)) {
                              gi  = i3(c,k,m,cc,nk);
                              pm[m] *= pski(binomN[s], count, Tski, hk[gi], 1.0);
                              
                              // for each detection, pdf(xy) | detected 
                              if ((pm[m] > minp) && (count>0)) {       // avoid underflow
                                  // retrieve hint = integral2D(zfn(x) over k)) 
                                  hint = hk[gi] / gsbval(c,0) * H[c];  
                                  for (j=start[w3]; j < start[w3]+count; j++) {
                                      pm[m] *= zcpp(j, m, c, gsbval, xy, mask) / hint;
                                  }
                              }
                          }
                          else {
                              pm[m] = 0.0;
                          }
                      }
                  }
              }
              if (dead==1) break;
          }
      }
  }    
  //==============================================================================
  
  std::vector<double> onehistorymm (int n) {
      std::vector<double> pm(mm, 1.0);
      if (binomN[0] < 0)
          prwpolygonXfxi(n,pm);
      else
          prwpolygonfxi(n,pm);
      for (int m=0; m<mm; m++) {
          pm[m] *= density(m,group[n]);
      }
      return pm;    // mm vector
  }
  //==============================================================================
  
  // function call operator that works for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
      std::vector<double> pm(mm, 1.0);
      for (std::size_t n = begin; n < end; n++) {
          pm = onehistorymm (n);
          for (int m=0; m<mm; m++) output(n,m) = pm[m];
      }
  }
  //==============================================================================
};

// [[Rcpp::export]]
NumericVector polygonfxicpp (
    const int           nc,
    const int           detectfn,
    const int           grain,
    const int           ncores,
    const double        minp,
    const IntegerVector binomN,
    const IntegerVector w,
    const NumericMatrix xy,
    const IntegerVector start,
    const IntegerVector group,
    const NumericVector hk,
    const NumericVector H,
    const NumericMatrix gsbval,
    const NumericMatrix pID,
    const NumericMatrix mask,
    const NumericMatrix density,
    const IntegerVector PIA,
    const NumericMatrix Tsk,
    const NumericMatrix h,
    const IntegerMatrix hindex, 
    const LogicalMatrix mbool
    
) {
    
    int mm = mask.nrow();
    NumericMatrix output(nc, mm); 
    
    // Construct and initialise
    polygonfxi somehist (nc, detectfn, grain, minp, binomN, w, xy, start, group, 
                         hk, H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, mbool,
                         output);

    if (ncores>1) {
     // Run operator() on multiple threads
     parallelFor(0, nc, somehist, grain, ncores);
   }
   else {
    // for debugging avoid multithreading and allow R calls e.g. Rprintf
    somehist.operator()(0,nc);
   }
  // Return consolidated result
  return output;
}
//==============================================================================
