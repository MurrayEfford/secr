#include <Rcpp.h>
#include "poly.h"
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
    const RVector<double> gsbR;
    const RVector<double> miscparmR;
    const double          w2;
    const RVector<int>    binomNR;
    RVector<double>       pdot;
    
    int  ss,kk;
    bool allsighting = true;
    std::vector<double> gsbS;
    std::vector<double> miscparmS;
    
    pdotpoint(
        const NumericMatrix &xy,
        const NumericMatrix &traps, 
        const NumericMatrix &dist2,
        const IntegerVector &detect,
        const NumericMatrix &Tsk, 
        const IntegerVector &markocc,
        const int           &fn,
        const NumericVector &gsb, 
        const NumericVector &miscparm,
        const double        &w2,
        const IntegerVector &binomN,
        NumericVector pdot)
        
        :  xyR(xy), trapsR(traps), dist2R(dist2), detectR(detect), TskR(Tsk), 
           markoccR(markocc), fn(fn), gsbR(gsb), miscparmR(miscparm), w2(w2), binomNR(binomN), 
           pdot(pdot) {
        
        ss = Tsk.ncol();
        kk = traps.nrow();
        for (int s=0; s<ss; s++) {
            if (markocc[s]>0) allsighting = false;                    /* capture occasions */
        }
        gsbS = as<std::vector<double>>(gsb);
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
        const NumericVector &gsb, 
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
    
    pdotpoint ppoint (xy, traps, dist2, detect, Tsk, markocc, fn, gsb, miscparm, w2, binomN, pdot);
    
    if (ncores>1) {
        parallelFor(0, nxy, ppoint, grain, ncores);
    }
    else {
        ppoint.operator()(0,nxy);    // for debugging avoid multithreading to allow R calls
    }
    
    return(wrap(pdot));
}

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
            H = hintegral1Ncpp(detectfn, as<std::vector<double>>(gsb));
        else
            H = hintegral2Ncpp(detectfn, as<std::vector<double>>(gsb));     
        
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
