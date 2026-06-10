#include "secr.h"

//==============================================================================
// 2019-08-10, 2019-10-18
// detector types multi, proximity, count, telemetry

// 2026-06-06 experimental use of LSE: convert pm to log(pm)
// 2026-06-07 rename gethr to gethrcpp
// 2026-06-08 telemscale now redundant as new LSE algrithm here and in 
//            secr_generalsecrloglikfn takes care of undeflow
// 2026-06-09 gethrcpp now fills only the required (maskused) elements of hr 

// [[Rcpp::export]]
NumericVector gethrcpp(
        const int nc,
        const int fn, 
        const IntegerVector &start,
        const NumericMatrix &xy, 
        const NumericMatrix &mask, 
        const IntegerVector &maskused,
        const NumericMatrix &gsbval) {
    // precompute all telemetry-to-mask detectfn values
    // 2026-06-06 convert to logarithm
    int c,j,m,t,hri;
    double r;
    int nt = xy.nrow();
    int mm = mask.nrow();
    int cc = gsbval.nrow();
    NumericVector gsb(3);  
    NumericVector hr(nt * mm * cc);
    fnptr lfnr;
    lfnr = getlfnr(fn);                  // log!                        
    for (c=0; c<cc; c++) {
        gsb(1) = gsbval(c,1);
        // normalising coefficient 1/c, c = 2.pi.sigma^2
        if ((fn==14) || (fn==16))
            gsb(0) = 1 / (2 * M_PI * gsb[1] * gsb[1]);
        else
            Rcpp::stop ("telemetry only coded for HHN and HEX");
        for (j=0; j<maskused.size(); j++) {
            m = maskused[j];
            for (t=0; t<nt; t++) {
                r = std::sqrt(d2cpp (m, t, mask, xy)); 
                hri = i3(c, m, t, cc, mm); 
                if (hri>1e8) Rcpp::stop ("c,m,t combinations exceed 1e8 in gethrcpp");
                hr[hri] = lfnr(gsb, r);
            }	    
        }
    }
    return (hr);
}
//=============================================================

struct simplehistories : public Worker {
    
    // input data
    const int   mm;
    const int   nc;
    const int   cc; // number of parameter combinations
    const int   grain;   
    const RVector<int>    binomN;     // s 
    const RVector<int>    markocc;    // s 
    const RVector<int>    firstocc;   // s 
    const RVector<double> pID;        // s 
    const RVector<int>    w;          // n,s,k, n,s (multi)
    const RVector<int>    group;      // g
    const RVector<double> gk; 
    const RVector<double> hk; 
    const RMatrix<double> density;    // n,g
    const RVector<int>    PIA;        // 1,n,s,k,1  for given x
    const RMatrix<double> Tsk;        // k,s
    const RMatrix<double> h;
    const RMatrix<int>    hindex;
    const RVector<int>    mask_indices;
    const RVector<int>    mask_offsets;
    const RVector<int>    mask_id;       // Maps individual to mask row
    const RVector<double> telemhr; 
    const RVector<int>    telemstart; 

    // working variables
    int  kk, ss;
    bool allX = true;
    
    // output 
    RVector<double> output;

    // Constructor to initialize an instance of Somehistories 
    // The RMatrix class can be automatically converted to from the Rcpp matrix type
    simplehistories(
        const int mm, 
        const int nc, 
        const int cc,
        const int grain,    
        const IntegerVector binomN,  
        const IntegerVector markocc,  
        const IntegerVector firstocc,  
        const NumericVector pID,  
        const IntegerVector w,
        const IntegerVector group,
        const NumericVector gk, 
        const NumericVector hk, 
        const NumericMatrix density,
        const IntegerVector PIA,
        const NumericMatrix Tsk,
        const NumericMatrix h,
        const IntegerMatrix hindex, 
        const IntegerVector mask_indices, 
        const IntegerVector mask_offsets,
        const IntegerVector mask_id,
        const NumericVector telemhr,
        const IntegerVector telemstart,
        NumericVector output
        )
        : 
        mm(mm), 
        nc(nc), 
        cc(cc), 
        grain(grain),
        binomN(binomN), 
        markocc(markocc), 
        firstocc(firstocc), 
        pID(pID), 
        w(w),
        group(group), 
        gk(gk), 
        hk(hk), 
        density(density), 
        PIA(PIA), 
        Tsk(Tsk), 
        h(h), 
        hindex(hindex), 
        mask_indices(mask_indices), 
        mask_offsets(mask_offsets), 
        mask_id(mask_id),
        telemhr(telemhr), 
        telemstart(telemstart),
        output(output) {
        
        // now can initialise these derived counts
        kk = Tsk.nrow();             // number of detectors
        ss = Tsk.ncol();             // number of occasions
        
        for (int s=0; s<ss; s++) if (binomN[s] != -2) allX = false;
    }
    //==============================================================================

    void fnucpp (const int n, const int s, int &cumcount, std::vector<double> &pm) {
        
        // Probability density of telemetry locations for animal n if it belongs to
        // latent class x and is centred at m. It is assumed that telemetry data correspond
        // to the last detector (nk)
        // 
        // detect is an array with the occasion-specific detector type in positions 0:(ss-1)
        // start is an array with the starting location of fixes for each animal in positions ss:(ss+nc-1).
        // 
        // The array hr contains precomputed detection function evaluations 
        // (dimension cc x mm x nt where nt is total number of telemetry locations).
        // 
        // [The size of hr can grow very large, and maybe there needs to be an option to
        //  compute hr dynamically.]
        // 
        // Assume detections are sorted by occasion.
        
        int c = 0;
        int hri = 0; 
        int i, j, m, t, w3;
        int count;
        int m_row = mask_id[n];
        if (telemstart[n+1] > telemstart[n]) {
            w3 = i3(n, s, kk-1, nc, ss);
            count = w[w3];  // number of telemetry fixes 
            if (count>0) {
                c = PIA[w3] - 1;                
                if (c<0) {
                    Rcpp::stop ("telemetry usage zero on telemetry occasion");
                }
                for (i=cumcount; i<(cumcount+count); i++) {
                    t = telemstart[n] + i;
                    for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
                        m = mask_indices[j];
                        hri  = i3(c, m, t, cc, mm); 
                        pm[m] += telemhr[hri];         // log from 2026-06-06
                        
                    }
                }
                cumcount += count;
            }
        }
    }
    
    //----------------------------------------------------------------------------
    
    void prwX (const int n, const int s, bool &dead, std::vector<double> &pm) {
        // multi-catch traps
        int c, j, k, m, w3;
        double H;
        int m_row = mask_id[n];
        
        if (allX) {
            k = w[s * nc + n];    // all detectors are traps, CH has been compressed to 2-D
            dead = k < 0;  
            k = abs(k)-1;         // trap number 0..kk-1; k = -1 if not caught 
        }
        else {
            // find site at which trapped
            k = -1;
            for (int ik=0; ik<kk; ik++) {
                if (w[i3(n, s, ik, nc, ss)] != 0) {
                    dead = ik < 0;
                    k = abs(ik);
                    break;
                }
            }             
        }
        // Not captured in any trap on occasion s 
        if (k < 0) {
            for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
                m = mask_indices[j];
                H = h(m, hindex(n,s));
                if (H > fuzz)
                    pm[m] -= H;    // log 2026-06-06
            }
        }
        // Captured in trap k on occasion s
        else {
            w3 = i3(n, s, k, nc, ss);   
            c = PIA[w3] - 1;
            if (c >= 0) {    // ignore unset traps 
                for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
                    m = mask_indices[j];
                    H = h(m, hindex(n,s));
                    if (H>fuzz) {
                        pm[m] += log(Tsk(k,s) * (1-exp(-H)) *  hk[i3(c, k, m, cc, kk)] / H);
                    }
                    else {
                        pm[m] = -huge;    // log 2026-06-06
                    }
                }
            }
        }
    }
    //----------------------------------------------------------------------------
    
    void prw (const int n, const int s, bool &dead, std::vector<double> &pm) {
        int c, j, k, m, w3, count;
        int m_row = mask_id[n];
        
        for (k=0; k<kk; k++) {   // over detectors
            w3 =  i3(n, s, k, nc, ss);
            c = PIA[w3] - 1;
            if (c >= 0) {    // ignore unused detectors 
                count = w[w3];
                if (count<0) {count = -count; dead = true; }
                for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
                    m = mask_indices[j];
                    // bug fix 2023-03-09 for Poisson counts
                    // log 2026-06-06
                    if (binomN[s]==0)
                        pm[m] += log(pski(binomN[s], count, Tsk(k,s), hk[i3(c, k, m, cc, kk)], pID[s]));  
                    else 
                        pm[m] += log(pski(binomN[s], count, Tsk(k,s), gk[i3(c, k, m, cc, kk)], pID[s]));  
                }
            }
        }
    }
    //----------------------------------------------------------------------------
    
    // cumulative probability animal n _not yet_ marked on successive occasions 1:ss
    // void getpdots (const int n, 
    //                std::vector<double> &pds) {
    //     int c, k, m, s;
    //     double pp;
    //     for (m=0; m<mm; m++) {
    //         pp = 1.0;
    //         for (s=0; s<(ss-1); s++) {
    //             if (markocc[s] > 0) {  // marking occasions only 
    //                 if (binomN[s] == -2) {
    //                     pp *= exp(-h(m, hindex(n,s)));
    //                 }   
    //                 else {
    //                     for (k=0; k< kk; k++) {
    //                         c = PIA[i3(n,s,k,nc,ss)] - 1;
    //                         if (c >= 0) {    // drops unset traps 
    //                             // pID always 1.0 on marking occasions
    //                             pp *= pski(binomN[s], 0, Tsk(k,s), gk[i3(c, k, m, cc, kk)], 1.0);  
    //                         }
    //                     }
    //                 }
    //             }
    //             pds[m*ss + s + 1] = pp;    // Pr(not marked) at next occasion
    //         }
    //     }
    // }
    //==============================================================================
    
    // site- and occasion-specific hazard of detection for animal n, x D 
    // void Hsk (const int n, 
    //           const std::vector<double> &pm, 
    //           std::vector<double> &hsk) {
    //     int c, k, m, s;
    //     for (int i=0; i>(ss*kk); i++) hsk[i]=0.0;
    //     for (s=0; s<ss; s++) {
    //         for (k=0; k<kk; k++) {   // over detectors
    //             c = PIA[i3(n, s, k, nc, ss)] - 1;
    //             if (c >= 0) {    // ignore unused detectors 
    //                 for (m=0; m<mm; m++) {
    //                     if (mbool(n,m)) {
    //                         if (binomN[s]<0)
    //                             hsk[s*kk+k] += pm[m] * Tsk(k,s) * gk[i3(c, k, m, cc, kk)];    
    //                         else
    //                             hsk[s*kk+k] += pm[m] * Tsk(k,s) * hk[i3(c, k, m, cc, kk)];    
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    //==============================================================================
    
    // site- and occasion-specific hazard of detection for animal n (most likely arbitrary)
    double hskm (const int n, 
                 const int s,
                 const int k,
                 const int m) {
        int c;
        double value = 0.0;
        c = PIA[i3(n, s, k, nc, ss)] - 1;
        if (c >= 0) {    // ignore unused detectors 
            value = Tsk(k,s) * hk[i3(c, k, m, cc, kk)];
        }
        return value;
    }
    //----------------------------------------------------------------------------
    
    double onehistory (std::size_t n) {
        bool dead = false;
        double sumpm = 0.0;
        double maxpm = -huge;
        int j, m;
        int cumcount = 0;
        int m_row = mask_id[n];
        
        std::vector<double> pm(mm, 0.0);    // log 2026-06-06  
        for (int s = 0; s < ss; s++) {      // over occasions
            // firstocc[n] used to exclude pre-marking sightings
            if (markocc[s]>0 || (s > firstocc[n])) {         
                if (binomN[s] == -2)        // multi-catch traps
                    prwX (n, s, dead, pm);           
                else if (binomN[s] >= -1) 
                    prw (n, s, dead, pm);  
                else if (binomN[s] == -3) 
                    fnucpp(n, s, cumcount, pm);
            }
            if (dead) break;               // out of s loop
        }
        
        for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
            m = mask_indices[j];
            pm[m] += log(density(m,group[n])); 
            if (pm[m]>maxpm) maxpm = pm[m];
        }
        
        // LSE trick 2026-06-06
        for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
            m = mask_indices[j];
            sumpm += exp(pm[m] - maxpm); 
        }
        sumpm = maxpm + log(sumpm);
        
        if (grain==0) {
            Rprintf("Debug n %zu sumpm %8.6e \n", n, sumpm);
        }
        
        // if (grain==0)
        //     Rprintf("n %4d sumpm %8.6g\n", n,sumpm);
        
            // expected number of unidentified marked sightings Tm
            // sum over marked animals
            //if (Tm && n==0) {
                // 
                // if (false) {
                //     std::vector<double> hsk(kk*ss, 0.0);   // for Tmmusk 
                //     Hsk(n, pm, hsk);
                //     
                //     if (grain==0) Rprintf("n %4d firstocc[n] %4d\n", n, firstocc[n]);
                //     
                //     for (int s=0; s<ss; s++) { 
                //         if ((markocc[s]<1) && (firstocc[n]<s)) {   // sighting occasions after first marking
                //             for (int k=0; k<kk; k++) {  
                //                 // Efford and Hunter 2018 Eqn 8
                //                 // summing over animals
                //                 Tmmusk(k,s) += (1-pID[s]) * hsk[s*kk+k] / sumpm;
                //                 
                //                 // if (grain==0) {
                //                 //     Rprintf("n %4d s %4d k %4d pID %8.6f hsk %8.6e sumpm %8.6e tmcomp %8.6e \n", 
                //                 //          n, s, k, pID[s], hsk[s*kk+k], sumpm, (1-pID[s]) * hsk[s*kk+k] / sumpm);
                //                 // }
                //             }
                //         }
                //     }
                // }
                // else 
        
        return sumpm; // may be -huge
    }
    //----------------------------------------------------------------------------
    
    // function call operator that works for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {        
        for (std::size_t n = begin; n < end; n++) {
            output[n] = onehistory (n);   // 2026-06-06 log
        }
    }
    //----------------------------------------------------------------------------
};

// [[Rcpp::export]]
NumericVector simplehistoriescpp (
        const int mm, 
        const int nc, 
        const int cc, 
        const int grain,   
        const int ncores,   
        const IntegerVector binomN, 
        const IntegerVector markocc, 
        const IntegerVector firstocc, 
        const NumericVector pID, 
        const IntegerVector w,
        const IntegerVector group,
        const NumericVector gk, 
        const NumericVector hk, 
        const NumericMatrix density,        // relative density - sums to 1.0
        const IntegerVector PIA, 
        const NumericMatrix Tsk, 
        const NumericMatrix h,
        const IntegerMatrix hindex, 
        const IntegerVector mask_indices,
        const IntegerVector mask_offsets,
        const IntegerVector mask_id,       // Maps individual to mask row
        const NumericVector telemhr,
        const IntegerVector telemstart)
    {
    
    NumericVector output(nc); 

    // Construct and initialise
    simplehistories somehist (mm, nc, cc, grain, 
                              binomN, markocc, firstocc, pID, w, 
                              group, gk, hk, density, PIA, Tsk, h, hindex, 
                              mask_indices, mask_offsets, mask_id, 
                              telemhr, telemstart, 
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
    // return output (log scale);
    return output;
    
}
//==============================================================================
