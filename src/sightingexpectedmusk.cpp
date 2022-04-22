#include <Rcpp.h>
#include "secr.h"

//==============================================================================
// 2019-09-14
// detector types multi, proximity, count

//struct expectedmusk : public Worker {
  class expectedmusk {
      public:
        
    // input data
    const int   nc;
    const int   cc; // number of parameter combinations
    const bool  Tu;  // true if data include sightings of unmarked animals
    const bool  Tm;  // true if data include unidentified sightings of marked animals
    const int sightmodel;
    const IntegerVector binomN;     // s 
    const IntegerVector markocc;    // s 
    const NumericVector pID;        // s 
    const IntegerVector group;      // g
    const NumericVector gk; 
    const NumericVector hk; 
    const NumericMatrix pi_density; // n x g
    const NumericMatrix Nm;         // n x g
    const IntegerVector PIA;        // 1, n, s, k, 1 (x given)
    const NumericMatrix Tsk;        // k x s
    const NumericMatrix h;
    const IntegerMatrix hindex;
    const NumericVector a0;

    // working variables
    int  kk, mm, ss;

    // output 
    NumericMatrix Tumusk;      // if Tu
    NumericMatrix Tmmusk;      // if Tm
    
    // Constructor to initialize an instance of expectedmusk 
    // The RMatrix class can be automatically converted to/from the Rcpp matrix type NumericMatrix
    expectedmusk(
        const int nc, 
        const int cc,
        const bool  Tu,              
        const bool  Tm,  
        const int sightmodel,
        const IntegerVector binomN,  
        const IntegerVector markocc,  
        const NumericVector pID,  
        const IntegerVector group,
        const NumericVector gk, 
        const NumericVector hk, 
        const NumericMatrix pi_density,
        const NumericMatrix Nm,
        const IntegerVector PIA,
        const NumericMatrix Tsk,
        const NumericMatrix h,
        const IntegerMatrix hindex, 
        const NumericVector a0,
        NumericMatrix Tumusk,
        NumericMatrix Tmmusk
    )  : 
        nc(nc), cc(cc), Tu(Tu), Tm(Tm), sightmodel(sightmodel),
        binomN(binomN), markocc(markocc), pID(pID), group(group), gk(gk), hk(hk), 
        pi_density(pi_density), Nm(Nm), PIA(PIA), Tsk(Tsk), h(h), hindex(hindex),
        a0(a0), Tumusk(Tumusk), Tmmusk(Tmmusk) {
        
        // now can initialise these derived counts
        kk = Tsk.nrow();             // number of detectors
        ss = Tsk.ncol();             // number of occasions
        mm = pi_density.nrow();      // number of mask points
        
    }
    //==============================================================================
    
    // cumulative probability animal n not yet marked on successive occasions 1:ss
    // given located at m
    void getpdots (const int n, 
                   std::vector<double> &pds) {
        int c, k, m, s;
        double pp;
        for (m=0; m<mm; m++) {
            pp = 1.0;
            for (s=0; s<(ss-1); s++) {
                if (markocc[s] > 0) {  // marking occasions only 
                    if (binomN[s] == -2) {
                        pp *= exp(-h(m, hindex(n,s)));
                    }   
                    else {
                        for (k=0; k< kk; k++) {
                            c = PIA[i3(n,s,k,nc,ss)] - 1;
                            if (c >= 0) {    // drops unset traps 
                                // pID always 1.0 on marking occasions
                                pp *= pski(binomN[s], 0, Tsk(k,s), gk[i3(c, k, m, cc, kk)], 1.0);  
                            }
                        }
                    }
                }
                pds[m*ss + s + 1] = pp;    // Pr(not marked) at next occasion
            }
        }
    }
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
    //==============================================================================
    
    void compute (NumericMatrix &Tumusk, NumericMatrix &Tmmusk) {
        
        // expected number of sightings
        // expected number of unmarked sightings Tu
        std::vector<double> pds(ss*mm, 0.0); 
        getpdots(0, pds);                        // representative animal n=0
        if (Tu) {
            for (int s=0; s<ss; s++) { 
                if (markocc[s]<1) {              // sighting occasions 
                    for (int k=0; k<kk; k++) {  
                        // Efford and Hunter 2018 Eqn 4
                        // summing over mask points for representative parameter values (arbitrary animal)
                        for (int m=0; m<mm; m++) {
                            if (sightmodel==0) 
                                Tumusk(k,s) += Nm(m,group[0]) *  pds[m*ss+s] * hskm(0,s,k,m);
                            else if (sightmodel==5)    // all pre-marked, number known
                                Tumusk(k,s) += (Nm(m,group[0]) - nc * pi_density(m, group[0])) * hskm(0,s,k,m);
                            else if (sightmodel==6)    // all pre-marked, number unknown
                                Tumusk(k,s) += (Nm(m,group[0]) - nc / a0[0] * pi_density(m, group[0])) * hskm(0,s,k,m);
                        }
                        // pID not relevant for unmarked sightings
                    }
                    
                }
            }
        }
        // expected number of unidentified marked sightings Tm
        // sum over marked animals
        if (Tm) {
            //     std::vector<double> hsk(kk*ss, 0.0);   // for Tmmusk 
            //     Hsk(n, pm, hsk);
            //     for (int s=0; s<ss; s++) { 
            //         if ((markocc[s]<1) && (firstocc[n]<s)) {   // sighting occasions after first marking
            //             for (int k=0; k<kk; k++) {  
            //                 // Efford and Hunter 2018 Eqn 8
            //                 // summing over animals
            //                 Tmmusk(k,s) += (1-pID[s]) * hsk[s*kk+k] / sumpm;
            //             }
            //         }
            //     }
            // }
            for (int s=0; s<ss; s++) { 
                if (markocc[s]<1) {              // sighting occasions 
                    for (int k=0; k<kk; k++) {  
                        // Efford and Hunter 2018 Section 3.3
                        // summing over mask points for representative parameter values (arbitrary animal)
                        for (int m=0; m<mm; m++) {
                            if (sightmodel==0) {         // not all sighting
                                Tmmusk(k,s) +=  Nm(m,group[0]) * (1-pds[m*ss+s]) * hskm(0,s,k,m);
                            }
                            else if (sightmodel==5) {   // all pre-marked, number known
                                Tmmusk(k,s) +=  nc * pi_density(m, group[0]) * hskm(0,s,k,m);
                            }
                            else if (sightmodel==6) {   // all pre-marked, number unknown
                                Tmmusk(k,s) +=  nc / a0[0] * pi_density(m, group[0]) * hskm(0,s,k,m);
                            }
                        }
                        Tmmusk(k,s) *= 1-pID[s];   // 2019-12-16
                    }
                }
            }
        }
    }
    //==============================================================================
    
};

// [[Rcpp::export]]
List expectedmucpp (
        const int nc,                 // number rows in CH (number marked when CH includes zeros - allsighting, knownmarks)
        const int cc, 
        const bool Tu,                // true if data include sightings of unmarked animals
        const bool Tm,                // true if data include unidentified sightings of marked animals
        const int sightmodel,
        const IntegerVector binomN, 
        const IntegerVector markocc, 
        const NumericVector pID, 
        const IntegerVector group,
        const NumericVector gk, 
        const NumericVector hk, 
        const NumericMatrix pi_density,    // relative density - sums to 1.0
        const NumericMatrix Nm,            // expected number per mask cell 
        const IntegerVector PIA, 
        const NumericMatrix Tsk, 
        const NumericMatrix h,
        const IntegerMatrix hindex,
        const NumericVector a0             // \int pi(x) p.(x) dx
        ) {
    
    NumericMatrix Tumusk(Tsk.nrow(), Tsk.ncol()); 
    NumericMatrix Tmmusk(Tsk.nrow(), Tsk.ncol()); 
    
    // Construct and initialise
    expectedmusk expectedmu (nc, cc, Tu, Tm, sightmodel, binomN, markocc, pID,  
                             group, gk, hk, pi_density, Nm, PIA, Tsk, h, hindex, 
                             a0, Tumusk, Tmmusk);
    
    //expectedmu.operator(); 
    expectedmu.compute (Tumusk, Tmmusk);
    
    // Return consolidated result
    // return output;
    return List::create(Named("Tumusk") = Tumusk, Named("Tmmusk") = Tmmusk);
    
}
//==============================================================================
