#include "poly.h"
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

//=============================================================
struct hdotpoly : public Worker {
    
    // input data
    const int detectfn;
    const bool convex;
    const int dim;
    const RVector<double> gsbR;
    const RMatrix<double> gsbvalR;
    const RVector<int> cumk;
    const RVector<int> markocc;
    const RMatrix<double> trapsR;
    const RMatrix<double> xyR;
    const RMatrix<double> TskR;
    RVector<double> hdot;
    
    double H;
    
    // output vector to write to
    
    int kk, nk, mm, npar, ss;
    bool allsighting = true;
    
    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    hdotpoly(const int detectfn,
             const bool convex,
             const int dim,
             const NumericVector gsb, 
             const NumericMatrix gsbval, 
             const IntegerVector cumk,
             const IntegerVector markocc,
             const NumericMatrix traps, 
             const NumericMatrix xy, 
             const NumericMatrix Tsk, 
             NumericVector hdot)
        : detectfn(detectfn), convex(convex), dim(dim), gsbR(gsb), gsbvalR(gsbval),
          cumk(cumk), markocc(markocc), trapsR(traps), xyR(xy), TskR(Tsk), 
          hdot(hdot) {
        
        nk = cumk.size()-1;
        npar = gsb.size();
        ss = Tsk.ncol();
        
        for (int s=0; s<ss; s++) {
            if (markocc[s]>0) allsighting = false;                    // capture occasions 
        }
        
        if (dim==1)
            H = hintegral1DNRcpp(detectfn, as<std::vector<double>>(gsb));
        else
            H = hintegral2DNRcpp(detectfn, as<std::vector<double>>(gsb));     
        
    }
    
    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end)
    {
        double hk = 0.0;
        double sumhk;
        double Tski;
        
        for (std::size_t i = begin; i<end; i++) {
            sumhk = 0.0;
            for (int s=0; s<ss; s++) {
                if ((markocc[s]>0) || allsighting) {                     
                    for (int k=0; k<nk; k++) {
                        Tski = TskR(k,s);
                        // subset traps rows cumk[k] : (cumk[k+1]-1)
                        int n1 = cumk[k];
                        int n2 = cumk[k+1]-1;
                        
                        if (Tski > 1e-10) {  
                            if (dim == 1)
                                hk = gsbR[0] * integral1DNRcpp (detectfn, i, 0, gsbvalR, trapsR, xyR, n1, n2) / H;
                            else // dim == 2
                                hk = gsbR[0] * integral2DNRcpp (detectfn, i, 0, gsbvalR, trapsR, xyR, n1, n2, convex) / H;
                            sumhk += hk * Tski;
                        }
                    }
                }
            }
            hdot[i] = sumhk;
        }
    }
};

// [[Rcpp::export]]
NumericVector hdotpolycpp (
        const NumericMatrix &xy, 
        const NumericMatrix &traps, 
        const NumericMatrix &Tsk, 
        const IntegerVector &markocc, 
        const IntegerVector &cumk, 
        const int &detectfn, 
        const NumericVector &gsb,
        const bool &convex,
        const int &dim,
        const int &grain,
        const int &ncores)
    {
    int nxy = xy.nrow();
    NumericMatrix gsbval(1,gsb.size());
    NumericVector hdot(nxy);
    for (int i=0; i<gsb.size(); i++) gsbval(0,i) = gsb(i);
    
    hdotpoly hpoly (detectfn, convex, dim, gsb, gsbval, cumk, markocc, traps, xy, Tsk, hdot);
    
    if (ncores>1) {
        parallelFor(0, nxy, hpoly, grain, ncores);
    }
    else {
        hpoly.operator()(0,nxy);    // for debugging avoid multithreading to allow R calls
    }
    
    return(wrap(hdot));
}

//=============================================================

