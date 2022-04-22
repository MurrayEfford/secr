#include "secr.h"     
using namespace std;
using namespace Rcpp;

//==============================================================================


// 'naive' functions are used to estimate auto initial values

// 2019-12-08 naivecap3cpp is new, simpler function for expected number of 
//            detections per detected animal. It is based on the hazard formulae 
//            Efford & Boulanger 2019.

// [[Rcpp::export]]
double naivedcpp (
    const double sigma,          // Parameter : detection scale 
    const IntegerVector &wt,      // integer trap weights 
    const NumericMatrix &traps,   // x,y locations of traps (first x, then y)  
    const NumericMatrix &animals, // x,y locations of animals (first x, then y) 
    const int    fn              // code 0 = halfnormal ONLY
)
{
  int    kk;      // number of traps
  int    nc;       // number of animals
  
  double truncate2 = (2.45 * sigma) * (2.45 * sigma);
  double sump  = 0;
  double sumdp = 0;
  double x,y;
  double dij, d21, d22, p1p2;
  int i,j,n;
  
  if (fn != 0)
    Rcpp::stop ("invalid detection function in external function naivedcpp");
  
  kk = traps.nrow();
  nc = animals.nrow();
  
  for (n=0; n<nc; n++)
  {
    x = animals[n];
    y = animals[n + nc];
    
    for (i=0; i<kk; i++) {
      if (wt[i] > 0) {
        for (j=0; j<(i-1); j++) {
          dij = (traps[i] - traps[j]) * (traps[i] - traps[j]) +
            (traps[i+kk] - traps[j+kk]) * (traps[i+kk] - traps[j+kk]);
          d21 = (traps[i] - x) * (traps[i] - x) + (traps[i+kk] - y) * (traps[i+kk] - y);
          d22 = (traps[j] - x) * (traps[j] - x) + (traps[j+kk] - y) * (traps[j+kk] - y);
          
          if ((d21<=truncate2) && (d22<=truncate2)) {
            p1p2 = exp(-(d21+d22) / 2 / sigma / sigma);
            if (wt[i] > 1)
              p1p2 = 1 - pow(1-p1p2, wt[i]);   
          }
          else
            p1p2 = 0;
          sump  += p1p2;
          sumdp += p1p2 * std::sqrt(dij);
        }
      }
    }
    for (i=0; i<kk; i++) {  // diagonal 
        d21 = (traps[i] - x) * (traps[i] - x) + (traps[i+kk] - y) * (traps[i+kk] - y);
        if (d21<=truncate2)                                     // d21=d22
	    sump += exp(-2*d21 /2 / sigma / sigma)/2;
    }
    
  }
  return(sumdp/sump);
}

//==============================================================================

// [[Rcpp::export]]
double naivecap3cpp (
        const int    detect,        // scalar code 0 = multicatch, 1 = proximity, 2 = count
        const double lambda0,       // Parameter : detection magnitude 
        const double sigma,         // Parameter : detection scale 
        const NumericMatrix &Tsk,    // trap usage k x s (applied directly to cumulative hazard)
        const NumericMatrix &traps, // x,y locations of traps (first x, then y)  
        const NumericMatrix &mask,  // x,y locations of mask points (first x, then y) 
        const int    fn             // only hazard halfnormal ()14) is supported at present      
)
{
    int kk = traps.nrow();  // number of traps
    int mm = mask.nrow();   // number of mask points
    int ss = Tsk.ncol();    // number of occasions
    
    double ls;
    double Lambdas = 0;
    double Lambda = 0;
    double captures = 0;
    double animals = 0;
    std::vector<double> h(kk,0);
    
    if (fn != 14)
        Rcpp::stop ("invalid detection function in naivecap3cpp");
    if (detect>2) 
        Rcpp::stop ("unrecognised detector in naivecap3cpp");
    
    for (int m=0; m<mm; m++) {
        for (int k=0; k<kk; k++) {
            h[k] = lambda0 * exp(-d2cpp(k, m, traps, mask) / 2 / sigma / sigma);
        }
        Lambda = 0;
        for (int s=0; s<ss; s++) {
            Lambdas = 0;
            for (int k=0; k<kk; k++) {
                ls = h[k] * Tsk(k,s);
                Lambdas += ls;
                if (detect==1) captures += 1 - exp(-ls);
            }
            if (detect==0) captures += 1 - exp(-Lambdas);
            else if (detect==2) captures += Lambdas;    // assuming Poisson counts
            Lambda += Lambdas; 
        }
        animals += 1 - exp(-Lambda);
    }
    if (animals<=0)
        return(0);    
    else
        return(captures/ animals);
}

/*==============================================================================*/
