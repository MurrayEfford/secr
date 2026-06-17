#include "secr.h"
#include <algorithm>

//==============================================================================
// 2019-09-05 detector type count
//            one occasion, binomial size from integer Tsk
//            no deaths, no groups
// 2020-04-05 fixed bug in indiv option
// 2021-05-11 replaced pm0base.reserve(mm) by pm0base.resize(mm) etc.
//            this is hoped to avoid test error
//            ERROR: AddressSanitizer: container-overflow on address 
//            in pr0()
// 2026-06-08 dump mbool, new individual masks
// 2026-06-15 ThreadRegistry registry;
// 2026-06-17 fillpr0 correction to pr0

struct fasthistories : public Worker {
    
    // input data
    const int   mm;
    const int   nc;
    const int   cc; // number of parameter combinations
    const int   grain;
    const bool  safeLL;
    const int   binomN;
    const bool  indiv;
    const RMatrix<int>    w;            // n x k
    const RMatrix<int>    ki;           // n x k
    const RVector<double> gk; 
    const RVector<double> hk; 
    const RVector<double> density;      // n
    const RVector<int>    PIA;
    const RVector<int>    Tsk;          // k
    const RVector<int>    mask_indices;
    const RVector<int>    mask_offsets;
    const RVector<int>    mask_id;      // Maps individual to mask row
    
    std::vector<double> pm0base;
    std::vector<double> pm0kbase;
    
    int kk;
    
    // output likelihoods
    RVector<double> output;
    ThreadRegistry& registry;
    
    // Workspace to hold thread-specific buffers
    struct ThreadWorkspace {
        std::vector<double> pm0, pm0k, pm;
        ThreadWorkspace(int mm, int kk) : pm0(mm), pm0k(kk * mm), pm(mm) {}
    };
    std::vector<ThreadWorkspace> workspaces;

    // Constructor to initialize an instance of Somehistories 
    // The RMatrix class can be automatically converted to from the Rcpp matrix type
    fasthistories(
        const int mm, 
        const int nc, 
        const int cc,
        const int grain,  
        const int ncores,
        
        const bool safeLL,
        const int binomN,
        const bool indiv,
        const IntegerMatrix w,
        const IntegerMatrix ki,
        
        const NumericVector gk, 
        const NumericVector hk, 
        const NumericVector density,
        const IntegerVector PIA,
        const IntegerVector Tsk,
        
        const IntegerVector mask_indices, 
        const IntegerVector mask_offsets,
        const IntegerVector mask_id,
        ThreadRegistry&     reg_in,
        NumericVector output)
        : 
        mm(mm), 
        nc(nc), 
        cc(cc), 
        grain(grain), 
        safeLL(safeLL),
        binomN(binomN), 
        indiv(indiv),
        w(w), 
        ki(ki), 
        gk(gk), 
        hk(hk), 
        density(density), 
        PIA(PIA), 
        Tsk(Tsk), 
        mask_indices(mask_indices), 
        mask_offsets(mask_offsets), 
        mask_id(mask_id),
        registry(reg_in),
        output(output) {
        
        kk = Tsk.size();        // assuming single occasion

        // Initialize workspaces based on available concurrency
        int n_threads = (ncores > 0) ? ncores : 1;
        for(int i = 0; i < n_threads; ++i) {
            workspaces.emplace_back(mm, kk);
        }
        
        pm0base.resize(mm);
        pm0kbase.resize(kk * mm);
        
        // pr0(0, pm0base, pm0kbase);
        ThreadWorkspace& ws = workspaces[0];
        fillpr0(0, pm0base, pm0kbase);
        
    }
    
    //==============================================================================
    
    void fillpr0 (int n, std::vector<double> &pm0, std::vector<double> &pm0k) { 
        
        int c, k, m, w3;
        std::fill(pm0.begin(), pm0.end(), 1.0);

        for (k=0; k<kk; k++) {
            w3 =  i3(n, 0, k, nc, 1);    // allow covariate of individual n or trap k 2019-12-06
            c = PIA[w3] - 1;
            if (c >= 0) {    // ignore unset traps
                for (m=0; m<mm; m++) {
                    if (binomN==0)
                        pm0k[kk * m + k] = gpois (0, Tsk[k] * hk[i3(c, k, m, cc, kk)]);
                    else
                        pm0k[kk * m + k] = gbinom (0, Tsk[k], gk[i3(c, k, m, cc, kk)]);
                    pm0[m] *= pm0k[kk * m + k];
                }
            }
        }
    }
    //==============================================================================
    void pr0 (const int n, ThreadWorkspace& ws) {
        std::vector<double>& pm0n = ws.pm0;
        std::vector<double>& pm0kn = ws.pm0k;
        fillpr0(n, pm0n, pm0kn);
    }
    //==============================================================================
    
    void prwL (const int n, ThreadWorkspace& ws) {
        int c, i, j, k, m, w3;
        
        // Use the buffers in the workspace
        std::vector<double>& pm0 = ws.pm0;
        std::vector<double>& pm0k = ws.pm0k;
        std::vector<double>& pm = ws.pm; 
        
        double pm0ktmp;
        bool base = (n==0) || !indiv;
        
        if (!base) { // (n>0) && indiv) {
            pr0(n, ws); // update for this animal - expt 2019-12-06, 2020-04-04
        }

        int m_row = mask_id[n];
        for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
            m = mask_indices[j];
            if (base) 
                pm[m] = pm0base[m];
            else 
                pm[m] = pm0[m];  
        }
        
        for (i=0; i<kk; i++) {
            k = ki(n,i);
            if (k<0) break;    // no more sites
            w3 =  i3(n, 0, k, nc, 1);
            c = PIA[w3] - 1;
            if (c >= 0) {    // ignore unset traps
                for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
                    m = mask_indices[j];
                    if (base) 
                        pm0ktmp = pm0kbase[kk * m + k]; 
                    else 
                        pm0ktmp = pm0k[kk * m + k];
                    if (binomN==0) {
                        pm[m] *= gpois (w(n,i), Tsk[k] * hk[i3(c, k, m, cc, kk)]) / pm0ktmp;
                    }
                    else {
                        pm[m] *= gbinom (w(n,i), Tsk[k], gk[i3(c, k, m, cc, kk)]) / pm0ktmp;
                    }
                }
            }
        }
    
    }
    //==============================================================================
    
    double onefasthistory (int n, ThreadWorkspace& ws) {
        std::fill(ws.pm.begin(), ws.pm.end(), 1.0);
        std::vector<double>& pm = ws.pm;
        prwL (n, ws);           
        int m_row = mask_id[n];
        double sumpm = 0.0;
        double maxpm = -huge;
        int j,m;
        
        if (safeLL) {
            for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
                m = mask_indices[j];
                pm[m] *= density[m];
                pm[m] = log(pm[m]);
                if (pm[m]>maxpm) maxpm = pm[m];
            }
            for (int j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
                int m = mask_indices[j];
                sumpm += exp(pm[m] - maxpm); 
            }
            sumpm = maxpm + log(sumpm);
        }
        else {
            for (j = mask_offsets[m_row]; j < mask_offsets[m_row+1]; ++j) {
                m = mask_indices[j];
                sumpm += pm[m] * density[m];
            }
            sumpm = log(sumpm);
        }
        
        return sumpm;
    }
    //==============================================================================
    
    // function call operator that works for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) { 
        int idx = registry.get_index();
        ThreadWorkspace& ws = workspaces[idx];
        for (std::size_t n = begin; n < end; n++) {
            output[n] = onefasthistory(n, ws);
        }
        //==============================================================================
    }
};

// [[Rcpp::export]]
NumericVector fasthistoriescpp (
        const int  mm, 
        const int  nc, 
        const int  cc, 
        const int  grain, 
        const int  ncores, 
        const bool safeLL,
        const int  binomN,
        const bool indiv,
        const IntegerMatrix w,
        const IntegerMatrix ki,
        const NumericVector gk, 
        const NumericVector hk, 
        const NumericVector density,
        const IntegerVector PIA, 
        const IntegerVector Tsk,
        const IntegerVector mask_indices,
        const IntegerVector mask_offsets,
        const IntegerVector mask_id       // Maps individual to mask row
        
        ) {
    
    NumericVector output(nc); 
    ThreadRegistry registry;
    
    // Construct and initialise
    fasthistories fasthist (
            mm, nc, cc, grain, ncores, 
            safeLL, binomN, indiv, w, ki, 
            gk, hk, density, PIA, Tsk, 
            mask_indices, mask_offsets, mask_id, registry, output); 
    
    if (ncores>1) {
        // Run operator() on multiple threads
        parallelFor(0, nc, fasthist, grain, ncores);
    }
    else {
        // for debugging avoid multithreading and allow R calls e.g. Rprintf
        fasthist.operator()(0,nc);    
    }
    
    // Return consolidated result
    // return output (log scale);
    return output;
}
//==============================================================================
