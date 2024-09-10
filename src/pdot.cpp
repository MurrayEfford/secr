#include "secr.h"

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

struct pdotpoint : public Worker {
    
    const RMatrix<double> xyR;
    const RMatrix<double> trapsR;
    const RMatrix<double> dist2R;                   
    const RVector<int>    detectR; 
    const RMatrix<double> TskR;
    const RVector<int>    markoccR;
    const int             fn;
    const RMatrix<double> gl0R;
    const RMatrix<double> sigR;
    const RVector<double> miscparmR;
    const double          w2;
    const RVector<int>    binomNR;
    RVector<double>       pdot;
    
    int  ss,kk;
    bool allsighting = true;
    std::vector<double> gsbS {0,0,0,0};
    std::vector<double> miscparmS;
    
    pdotpoint(
        const NumericMatrix &xy,
        const NumericMatrix &traps, 
        const NumericMatrix &dist2,
        const IntegerVector &detect,
        const NumericMatrix &Tsk, 
        const IntegerVector &markocc,
        const int           &fn,
        const NumericMatrix &gl0,            // new
        const NumericMatrix &sig,            // new
        const NumericVector &otherdetpar,    // new    
        const NumericVector &miscparm,
        const double        &w2,
        const IntegerVector &binomN,
        NumericVector pdot)
        
        :  xyR(xy), trapsR(traps), dist2R(dist2), detectR(detect), TskR(Tsk), 
           markoccR(markocc), fn(fn), gl0R(gl0), sigR(sig), miscparmR(miscparm), 
           w2(w2), binomNR(binomN), pdot(pdot) {
        
        ss = Tsk.ncol();
        kk = traps.nrow();
        for (int s=0; s<ss; s++) {
            if (markocc[s]>0) allsighting = false;    /* capture occasions */
        }
        
        gsbS[2] = otherdetpar[2];    // assumed constant
        gsbS[3] = otherdetpar[3];    // assumed constant
        miscparmS = as<std::vector<double>>(miscparm);
    }
    
    double onepoint (int i) {
        double tempval = 1.0;
        double p;
        double Tski = 1.0;
        for (int s=0; s<ss; s++) {
            if (((markoccR[s]>0) || allsighting) && (detectR[s]!=13)) {                     
                for (int k=0; k<kk; k++) {
                    Tski = TskR(k,s);
                    if (Tski > 1e-10) {
                        gsbS[0] = gl0R(k,s);
                        gsbS[1] = sigR(k,s);
                        p = pfnS(fn, dist2R(k,i), gsbS, miscparmS, w2);
                        /* counts */
                        if (detectR[s] == 2) {    
                            if (binomNR[s] == 0)
                                p = 1 - countp(0, 0, Tski * hazard(p));
                            else if (binomNR[s] == 1)
                                p = 1 - countp(0, round(Tski), p);
                            else {
                                if (fabs(Tski-1) > 1e-10)
                                    p = 1 - pow(1-p, Tski);
                                p = 1 - countp(0, binomNR[s], p);
                            }
                        }
                        else {                
                            if (fabs(Tski-1) > 1e-10)
                                p = 1 - pow(1-p, Tski);
                        }
                        tempval *= 1 - p;
                    }
                }
            }
        }
        return (1 - tempval);   
    }
    
    // function call operator that works for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end)
    {
        for (std::size_t i = begin; i<end; i++) {
            pdot[i] = onepoint(i);
        }
    }
};

// [[Rcpp::export]]
NumericVector pdotpointcpp (
        
        const NumericMatrix &xy, 
        const NumericMatrix &traps, 
        const NumericMatrix &dist2,                   
        const IntegerVector &detect, 
        const NumericMatrix &Tsk, 
        const IntegerVector &markocc, 
        const int &fn, 
        const NumericMatrix &gl0,            // new
        const NumericMatrix &sig,            // new
        const NumericVector &otherdetpar,     // new 2-vector    
        const NumericVector &miscparm,
        const double &w2, 
        const IntegerVector &binomN,
        const int &grain,
        const int &ncores)
    
{
    
    if (anypolygon(detect) || anytransect(detect)) {
        Rcpp::stop("pdotpoint not for polygon or transect detectors");
    }
    if (fn>20) {
        Rcpp::stop("pdotpointcpp requires detectfn < 21");
    }
    
    int nxy = xy.nrow();
    NumericVector pdot(nxy);
    
    pdotpoint ppoint (xy, traps, dist2, detect, Tsk, markocc, fn, 
                      gl0, sig, otherdetpar, miscparm, w2, binomN, pdot);
    
    if (ncores>1) {
        parallelFor(0, nxy, ppoint, grain, ncores);
    }
    else {
        ppoint.operator()(0,nxy);    // for debugging avoid multithreading to allow R calls
    }
    
    return(wrap(pdot));
}

