#include "secr.h"
// using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

struct nkpoint : public Worker {
    
    const RVector<double> DR;
    const RMatrix<double> dist2R;                
    const RVector<int>    detectR; 
    const RMatrix<double> TskR;
    const RVector<int>    markoccR;
    const int             fn;
    const RVector<double> gsbR;
    const RVector<double> miscparmR;
    const double          w2;
    const RVector<int>    binomNR;
    RVector<double>       nk;
    RMatrix<double>       H;
    
    int  kk, ss, mm;
    bool allsighting = true;
    bool multicatch = false;
    std::vector<double> gsbS;
    std::vector<double> miscparmS;
    
    nkpoint(
        const NumericVector &D,
        const NumericMatrix &dist2,
        const IntegerVector &detect,
        const NumericMatrix &Tsk, 
        const IntegerVector &markocc,
        const int           &fn,
        const NumericVector &gsb, 
        const NumericVector &miscparm,
        const double        &w2,
        const IntegerVector &binomN,
        NumericVector &nk,
        NumericMatrix &H)
        
        :  DR(D), dist2R(dist2), detectR(detect), TskR(Tsk), 
            markoccR(markocc), fn(fn), gsbR(gsb), miscparmR(miscparm), w2(w2), binomNR(binomN), 
            nk(nk), H(H) {
        
        ss = Tsk.ncol();
        kk = dist2.nrow();
        mm = dist2.ncol();
        int s,k,m;
        
        for (s=0; s<ss; s++) {
            if (markocc[s]>0) allsighting = false;                    /* capture occasions */
            if (detect[s]==0) multicatch = true;
        }
        gsbS = as<std::vector<double>>(gsb);
        miscparmS = as<std::vector<double>>(miscparm);

        if (multicatch) {
            /* compute total hazard of detection for each mask point m and occasion s */
            for (m=0; m<mm; m++) {
                for (s=0; s<ss; s++) {
                    for (k=0; k<kk; k++) {
                        H(m,s) += Tsk(k,s) * hazard(pfnS(fn, dist2(k,m), gsbS, miscparmS, w2));
                    }
                }
            }
        }
    }
    
    double onetrap (int k) {
        double tempval = 1.0;
        double tempsum = 0.0;
        double p;
        double Tski = 1.0;
        for (int m=0; m<mm; m++) {
            tempval = 1.0;
            for (int s=0; s<ss; s++) {
                if (((markoccR[s]>0) || allsighting) && (detectR[s]!=13)) {                     
                    Tski = TskR(k,s);
                    if (Tski > 1e-10) {
                        p = pfnS(fn, dist2R(k,m), gsbS, miscparmS, w2);
                        
                        if (detectR[s] == 0) {         
                            /* multi-catch trap */
                            if (H(m,s) > 0) {
                                p = (1 - std::exp(-H(m,s))) * Tski * hazard(p) / H(m,s);
                            }
                            else {
                                p = 0;
                            }
                        }
                    else if (detectR[s] == 2) {    
                        /* counts */
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
                        /* binary proximity detector */
                        if (fabs(Tski-1) > 1e-10)
                            p = 1 - pow(1-p, Tski);
                    }
                    tempval *= 1 - p;
                    }
                }
            }
            tempsum += (1-tempval) * DR[m];
        }
        return (tempsum);   
    }
    
    // function call operator that works for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end)
    {
        for (std::size_t k = begin; k<end; k++) {
            nk[k] = onetrap(k);
        }
    }
};

// [[Rcpp::export]]
NumericVector nkpointcpp (
        
        const NumericVector &D, 
        const NumericMatrix &dist2,                   
        const IntegerVector &detect, 
        const NumericMatrix &Tsk, 
        const IntegerVector &markocc, 
        const int           &fn, 
        const NumericVector &gsb, 
        const NumericVector &miscparm,
        const double        &w2, 
        const IntegerVector &binomN,
        const int           &grain,
        const int           &ncores)
    
{
    
    if (anypolygon(detect) || anytransect(detect)) {
        Rcpp::stop("nkpoint not for polygon or transect detectors");
    }
    if (fn>20) {
        Rcpp::stop("nkpointcpp requires detectfn < 21");
    }
    
    int ss = Tsk.ncol();
    int kk = dist2.nrow();
    int mm = dist2.ncol();
    NumericVector nk(kk);
    NumericMatrix H(mm,ss);
    
    nkpoint npoint (D, dist2, detect, Tsk, markocc, fn, gsb, miscparm, w2, binomN, nk, H);
    
    if (ncores>1) {
        parallelFor(0, kk, npoint, grain, ncores);
    }
    else {
        npoint.operator()(0, kk);    // for debugging avoid multithreading to allow R calls
    }
    
    return(wrap(nk));
}
