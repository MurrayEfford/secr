// Functions related to sighting likelihood October-December 2015, September 2019

// EH eq refers Efford & Hunter ms 2016-11-18 

// like 1 conditional likelihood
// like 3 unused here
// like 4 unused here
// like 5 all sighting, n0 known
// like 6 all sighting, n0 unknown

#include "secr.h"

//==============================================================================

// [[Rcpp::export]]
List Tsightinglikcpp (
        const Rcpp::IntegerMatrix &T,           // sighting count(s)
        const Rcpp::IntegerVector &markocc,     // distinguish sighting and marking occasions
        const Rcpp::IntegerVector &binomN, 
        const Rcpp::NumericMatrix &Tsk,         // usage
        const Rcpp::NumericMatrix &musk,        // expected count
        const int debug) {
    
    int s,k;
    double tempmu;
    
    // codes for count aggregation 
    // default is no aggregation (separate T_sk each occasion and detector) 
    bool TPooled;
    bool TBydetector; 
    
    int TCsk= 0;
    int nused = 0;   // number of detectors with non zero effort 
    double summu = 0;
    
    int ss = Tsk.ncol(); 
    int nk = Tsk.nrow();
    
    std::vector<int> nusedk;
    std::vector<double> summuk;
    int nsight = 0;     // number of sighting occasions 
    int firstsightocc;
    
    double Tlik = 0; 
    int T0 = T.size();
    
    //------------------------------------------------------------------------------
    // set working variables 
    firstsightocc = ss+1;
    TPooled = (T0 == 1);      // TPooled == 1 if there is a single summed count    
    TBydetector = (T0 == nk); // TBydetector == 1 if counts are summed by detector   
    if (TBydetector) {
        summuk.resize(nk);
        nusedk.resize(nk);
    }
    if (debug>1) {
        Rprintf("TPooled %4d \n", TPooled);
        Rprintf("TBydetector %4d \n", TBydetector);
    }
    //-----------------------------------------------------------------------------
    // loop over occasions 
    for (s=0; s < ss; s++) {
        if (markocc[s] < 1) {     // sighting occasions only 
            nsight += 1;
            if (s < firstsightocc) firstsightocc = s;
            for (k=0; k < nk; k++) {
                tempmu = musk[s * nk + k];
                // Tsk is effort; no relation to TCsk! 
                nused += Tsk[s * nk + k]>0;
                
                //-----------------------------------------------------------------
                // compute likelihood for this cell 
                if (!TPooled && !TBydetector) {
                    TCsk = T[s * nk + k];
                    if ((TCsk>0) && (tempmu<=0)) {
                        // Rprintf("TCsk = %d tempmu = %8.4f\n", TCsk, tempmu);
                        // Rprintf ("zero sighting probability when T number >0\n");
                        return List::create(Named("resultcode") = 53);
                    }
                    
                    // 2022-03-05 new condition to trap zero-usage case
                    // otherwise no data, none expected, no change (zero usage)
                    if (tempmu > 0) {
                        // binary (multi, proximity) 
                        if (binomN[s]<0) {
                            if (TCsk>1) TCsk = 1;
                            if (tempmu>0) { 
                                // Tlik += R::dbinom(TCsk, 1, 1-exp(-tempmu), 1);
                                boost::math::bernoulli_distribution<> bern(1-exp(-tempmu));
                                Tlik += log(pdf(bern,TCsk));
                            }
                        }
                        // count 
                        else {
                            // Tlik += R::dpois(TCsk,  tempmu, 1);
                            boost::math::poisson_distribution<> pois(tempmu);
                            Tlik += log(pdf(pois,TCsk));   // tempmu)); bugfix 2021-10-17
                        }
                    }      
                    
                    if (std::isnan(Tlik) || (Tlik < -1e6)) {
                        // Rprintf("very negative or NaN Tlik in Tsightinglik\n");
                        return List::create(Named("resultcode") = 54);
                    }
                }
                
                if (TPooled) {
                    nused += Tsk[s * nk + k]>0;
                    if (markocc[s] == 0) {
                        if (binomN[s] < 0)
                            summu += 1-exp(-tempmu); // summing probabilities? 
                        else
                            summu += tempmu;
                    }
                    // else markocc < 0, marked animals not distinguished so no increment 
                }
                if (TBydetector) {
                    nusedk[k] += Tsk[s * nk + k]>0;
                    if (markocc[s] == 0) {
                        if (binomN[s]<0) 
                            summuk[k] += 1-exp(-tempmu);
                        else
                            summuk[k] += tempmu;
                    }
                    // ??? else markocc < 0, marked animals not distinguished so no increment 
                }
                //-----------------------------------------------------------------
            }  // end k loop     
        }    
    } // end s loop 
    //-----------------------------------------------------------------------------
    
    // Likelihood for pooled counts only 
    if (TBydetector) {
        for (k=0; k<nk; k++) {
            // Increment likelihood for counts by detector 
            // sum over s for detector k 
            if (debug>0) 
                Rprintf("k %4d sumT %4d summu %8.3f nused %4d \n", k, T[k], summu, nused);
            if (binomN[firstsightocc] < 0) {  
                // assume sighting detector same all occasions 
                // and p constant over occasions (using arithmetic mean here) 

                // 2022-03-05 new condition to trap zero-usage case
                if (summuk[k]>0) {
                    // Tlik += R::dbinom(T[k],  nusedk[k], summuk[k] / nsight, 1);  // 2017-03-17 
                    boost::math::binomial_distribution<> bin(nusedk[k], summuk[k] / nsight);
                    Tlik += log(pdf(bin,T[k]));
                }
            }
            else {
                // Tlik += R::dpois(T[k],  summuk[k], 1);
                // 
                if (summuk[k]>0) {
                    boost::math::poisson_distribution<> pois(summuk[k]);
                    Tlik += log(pdf(pois,T[k]));
                }
            }
        }
    }
    else if (TPooled) {
        // first input is the sum over s,k 
        if (debug>0) Rprintf("sumT %4d summu %8.3f nused %4d \n", T[1], summu, nused);
        if (binomN[firstsightocc] < 0) {  // assume sighting detectors same on all occasions 
            // Tlik = R::dbinom(T[1],  nused, summu, 1);  // weird use of summu - to be fixed 
            boost::math::binomial_distribution<> bin(nused, summu);
            Tlik += log(pdf(bin,T[1]));
        }
        else {
            // Tlik = R::dpois(T[1],  summu, 1);
            boost::math::poisson_distribution<> pois(summu);
            Tlik += log(pdf(pois,T[1]));
        }
    } 
    
    return List::create(Named("resultcode") = 0, Named("Tlik") = Tlik);
    
}
//==============================================================================

