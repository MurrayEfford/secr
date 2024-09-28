#include "secr.h"

//==============================================================================
// 2019-08-10

struct simplehistoriesfxi : public Worker {
  
  // input data
  const int   x;
  const int   mm;
  const int   nc;
  const int   cc; // number of parameter combinations
  const int   grain;
  const RVector<int>    binomN;     // s 
  const RVector<int>    w;          // n x s x k or n x s (multi)
  const RVector<int>    group;      // g
  const RVector<double> gk; 
  const RVector<double> hk; 
  const RMatrix<double> density;    // n x g
  const RVector<int>    PIA;
  const RMatrix<double> Tsk;        // k x s
  const RMatrix<double> h;
  const RMatrix<int>    hindex;
  
  // working variables
  int  kk, ss;

  // output likelihoods
  RMatrix<double> output;
  
  // Constructor to initialize an instance of simplehistoriesfxi 
  // The RMatrix class can be automatically converted to from the Rcpp matrix type
  simplehistoriesfxi(
    int x, 
    int mm, 
    int nc, 
    int cc,
    int grain,                    
    const IntegerVector binomN,  
    const IntegerVector w,
    const IntegerVector group,
    const NumericVector gk, 
    const NumericVector hk, 
    const NumericMatrix density,
    const IntegerVector PIA,
    const NumericMatrix Tsk,
    const NumericMatrix h,
    const IntegerMatrix hindex, 
    NumericMatrix output)    
    : 
    x(x), mm(mm), nc(nc), cc(cc), grain(grain), 
    binomN(binomN), w(w), group(group), gk(gk), hk(hk), density(density), PIA(PIA), 
    Tsk(Tsk), h(h), hindex(hindex), 
    output(output) {
    
    // now can initialise these derived counts
    kk = Tsk.nrow();             // number of detectors
    ss = Tsk.ncol();             // number of occasions
    
  }
  //==============================================================================
  
  void prw (const int n, std::vector<double> &pm) {
    int c, gi, k, m, s, wi, wxi, count;
    double Tski,H;
    bool dead = false;
    for (s = 0; s < ss; s++) {   // over occasions
      if (binomN[s] == -2) {     // multi-catch traps
        wi = w[s * nc + n]; 
        dead = wi < 0;  
        k = abs(wi)-1;         // trap number 0..kk-1; k = -1 if not caught 
        // Not captured in any trap on occasion s 
        if (k < 0) {
          for (m=0; m<mm; m++) {
              H = h(m, hindex(n,s));
              if (H > fuzz)
		  pm[m] *= exp(-H);
          }
        }
        // Captured in trap k on occasion s
        else {
          wxi = i4(n, s, k, x, nc, ss, kk);   
          c = PIA[wxi] - 1;
          if (c >= 0) {    // ignore unset traps 
            Tski = Tsk(k,s);
            for (m=0; m<mm; m++) {
                H = h(m, hindex(n,s));
                gi  = i3(c, k, m, cc, kk);
                pm[m] *=  Tski * (1-exp(-H)) *  hk[gi] / H;
            }
          }
        }
      }
      else {  // all other point detectors
        for (k=0; k<kk; k++) {
          wxi =  i4(n, s, k, x, nc, ss, kk);
          c = PIA[wxi] - 1;
          if (c >= 0) {    // ignore unset traps 
            Tski = Tsk(k,s);
            wi = i3(n, s, k, nc, ss);
            count = w[wi];
            if (count<0) {count = -count; dead = true; }
            for (m=0; m<mm; m++) {
                gi  = i3(c, k, m, cc, kk);
                pm[m] *= pski(binomN[s], count, Tski, gk[gi], 1.0);
            }
          }
        }
      }
      if (dead) break;   // out of s loop
    }
  }
  //==============================================================================
  
  std::vector<double> onehistorymm (int n) {
    std::vector<double> pm(mm, 1.0);
    prw (n, pm);           
    for (int m=0; m<mm; m++) {
      pm[m] *= density(m,group[n]); 
    }
    return pm;   // mm vector
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
NumericMatrix simplehistoriesfxicpp (
    const int x, 
    const int mm, 
    const int nc, 
    const int cc, 
    const int grain, 
    const int ncores, 
    const IntegerVector binomN, 
    const IntegerVector w,
    const IntegerVector group,
    const NumericVector gk, 
    const NumericVector hk, 
    const NumericMatrix density,
    const IntegerVector PIA, 
    const NumericMatrix Tsk, 
    const NumericMatrix h,
    const IntegerMatrix hindex
    ) {
  
  NumericMatrix output(nc, mm); 
  
  // Construct and initialise
  simplehistoriesfxi somehist (x, mm, nc, cc, grain, 
                               binomN, w, group, gk, hk, 
                              density, PIA, Tsk, h, hindex, output);
  
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
