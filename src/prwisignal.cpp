#include "secr.h"
using namespace RcppParallel;

//==============================================================================
// 2019-08-18
//
struct signalhistories : public Worker {
    // input data
    int   mm;
    int   nc;
    int   detectfn;
    int   grain;
    const RVector<int>    binomN;     // s
    const RVector<int>    w;          // n x s x k
    const RMatrix<double> signal;
    const RVector<int>    group;      // n
    const RVector<double> gk;
    const RMatrix<double> gsbval;
    const RMatrix<double> dist2;
    const RMatrix<double> density;    // n x g
    const RVector<int>    PIA;        // 1,n,s,k,1   for given x
    const RVector<double> miscparm;
    const RMatrix<int>    mbool;      // appears cannot use RMatrix<bool>
    
    // working variables
    int  kk, ss, cc;
    
    // output likelihoods
    RVector<double> output;
    
    // Constructor to initialize an instance of Somehistories
    // The RMatrix class can be automatically converted to from the Rcpp matrix type
    signalhistories(
        int mm,
        int nc,
        int detectfn,
        int grain,
        const IntegerVector binomN,
        const IntegerVector w,
        const NumericMatrix signal,
        const IntegerVector group,
        const NumericVector gk,
        const NumericMatrix gsbval,
        const NumericMatrix dist2,
        const NumericMatrix density,
        const IntegerVector PIA,
        const NumericVector miscparm,
        const LogicalMatrix mbool,
        NumericVector output)
        :
        mm(mm), nc(nc), detectfn(detectfn), grain(grain), 
        binomN(binomN), w(w), signal(signal), group(group), gk(gk), gsbval(gsbval), 
        dist2(dist2), density(density), PIA(PIA), miscparm(miscparm),
        mbool(mbool), output(output) {
        // now can initialise these derived counts
        kk = dist2.nrow();          // number of detectors
        ss = 1;                     // number of occasions
        cc = gsbval.nrow();         // number of parameter combinations
    }
    //==============================================================================
    // local mufnL uses RMatrix input
    double mufnL (
            const int k,
            const int m,
            const double b0,
            const double b1,
            const RMatrix<double> &dist2,
            const bool spherical)
        // Return predicted signal strength at k for source at point m,
        // given strength at source of b0 dB and attenuation of b1 dB/m.
        // Spherical spreading is included if spherical > 0
        // Uses distance lookup in dist2
        
    {
        double d2val;
        d2val = dist2(k,m);
        if (spherical <= 0)
            return (b0 + b1 * std::sqrt(d2val));
        else {
            if (d2val>1) {
                return (b0 - 10 * log ( d2val ) / 2.302585 + b1 * (std::sqrt(d2val)-1)); 
            }
            else
                return (b0);
        }
    }
    
    // hazard (fn 14:19) or g ()fn 0:11) (distance lookup)
    double zLcpp (
            const int c,
            const int k,
            const int m)
    {
        double r2 = dist2(k,m);
        double r = std::sqrt(r2);
        double mu, gam;
        double temp;
        if (detectfn == 9) {    // binary signal strength
            // b0 = gsbval(c,0)
            // b1 = gsbval(c,1)
            temp = -(gsbval(c,0) + gsbval(c,1) * r);
            // return (R::pnorm(temp,0,1,0,0));    // upper
            boost::math::normal_distribution<> n;
            return (boost::math::cdf(complement(n,temp)));    // upper
        }
        else if (detectfn == 10 || detectfn == 11) {   // signal strength, signal strength spherical
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
            return (boost::math::cdf(complement(n,gam)));    // upper
        }
        else (Rcpp::stop("unknown or invalid detection function"));
    }
    void prwsignal (const int n, std::vector<double> &pm) {
        int c, gi, k, m, s, w3, count;
        double sig, mu, sdS;
        for (s=0; s < ss; s++) {   // over occasions
            for (k=0; k<kk; k++) {
                w3 = i3(n, s, k, nc, ss);
                c = PIA[w3] - 1;
                if (c >= 0) {    // ignore unset traps
                    count = abs(w[w3]);
                    if (count == 0) {   // not detected at this mic
                        for (m=0; m<mm; m++) {
                            if (mbool(n,m)) {
                                gi  = i3(c, k, m, cc, kk);
                                pm[m] *= pski(binomN[s], 0, 1, gk[gi], 1.0);
                            }
                            else {
                                pm[m] = 0.0;
                            }
                        }
                    }
                    else {   // detected at this mic
                        sig = signal(n,k);
                        for (m=0; m<mm; m++) {
                            if (mbool(n,m)) {
                                if (sig >= 0) {
                                    // valid measurement of signal
                                    mu  = mufnL (k, m, gsbval(c,0), gsbval(c,1), dist2, detectfn==11);
                                    sdS = gsbval(c,2);
                                    // pm[m] *= R::dnorm((sig - mu), 0, sdS, 0);
                                    boost::math::normal_distribution<> n(0,sdS);
                                    pm[m] *= boost::math::pdf(n,(sig - mu));
                                }
                                else  {
                                    // signal value missing; detection only
                                    gi = i3(c,k,m,cc,kk);
                                    pm[m] *= countp (1, binomN[s], gk[gi]);
                                }
                            }
                            else {
                                pm[m] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }
    
    //==============================================================================
    double onehistorycpp (int n) {
        double prwi;
        std::vector<double> pm(mm, 1.0);
        prwsignal(n,pm);
        for (int m=0; m<mm; m++) {
            pm[m] *= density(m,group[n]);
        }
        prwi = std::accumulate(pm.begin(), pm.end(), 0.0);
        return prwi;    // may be zero
    }
    // function call operator that works for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t n = begin; n < end; n++) {
            output[n] = onehistorycpp (n);
        }
    }
};

// [[Rcpp::export]]
NumericVector signalhistoriescpp (
        const int mm,
        const int nc,
        const int detectfn,
        const int grain,
        const int ncores,
        const IntegerVector binomN,
        const IntegerVector w,
        const NumericMatrix signal,
        const IntegerVector group,
        const NumericVector gk,
        const NumericMatrix gsbval,
        const NumericMatrix dist2,
        const NumericMatrix density,
        const IntegerVector PIA,
        const NumericVector miscparm,
        const LogicalMatrix mbool) {
    
    NumericVector output(nc);
    
    // Construct and initialise
    signalhistories somehist (mm, nc, detectfn, grain, binomN, w, signal, 
                              group, gk, gsbval, dist2, density, PIA, miscparm, mbool, 
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
