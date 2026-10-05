#include "secr.h"
#include <RcppParallel.h>
#include <random>

//==============================================================================
// detector types multi, proximity, count

// for sightmodel = 7 we do not distinguish marked and unmarked animals 
// c-hat is for the total count of detections C
// all counts etc. are in 'unmarked' arrays

// 2026-09-17 Switch to C++11 <random> with thread-Local RNG
// 2026-09-17 'seed' argument (integer or double converted by Rcpp to uint64_t)

struct chat : public Worker {
    
    // input data
    const int    mm;                  // number of mask cells
    const int    nmark;
    const int    cc0;                 // number of parameter combinations
    const int    grain;
    const int    nsim;                // number of replicate simulations for chat 
    const int    sightmodel;          // 5 allsighting known n0, 6 allsighting unknown n0
    const double sumD;
    const double area;                // area of mask cell in hectares
    const int    distrib;             // 0 Poisson, 1 fixed
    const int    anytelem;
    const RVector<int>    binomN;     // s 
    const RVector<int>    markocc;    // s 
    const RMatrix<double> pID;        // s 
    const RVector<int>    group;      // g
    const RVector<double> gk0; 
    const RVector<double> hk0; 
    const RMatrix<double> density;    // n x g
    const RVector<int>    PIA0;
    const RMatrix<double> Tsk;        // k x s
    const RVector<double> pmix;
    
    // output
    RMatrix<double> chatmat;
    RVector<double> chatout;

    uint64_t base_seed;     // declared after the outputs to match initialiser order (-Wreorder)

    // working variables
    int  kk, k1, ss, resultcode;
    double sumNm;
    bool allsighting = false;
    const double tol = 1e-6;    
    int nmix = 1;
    
    // std::vector<double> pds; 
    std::vector<double> Nm; 
    std::vector<double> cumprob; 
    std::vector<double> a0; 
    
    // Constructor to initialize an instance of Somehistories 
    chat(
        const int mm, 
        const int nmark, 
        const int cc0,
        const int grain,                    
        const int nsim,        // number of replicate simulations for chat 
        const int sightmodel,  // 5 allsighting known n0, 6 allsighting unknown n0, 7 telemetry 'marking'
        const double sumD,
        const double area,
        const int distrib,
        const int anytelem,
        const IntegerVector binomN,  
        const IntegerVector markocc,  
        const NumericMatrix pID,  
        const IntegerVector group,
        const NumericVector gk0, 
        const NumericVector hk0, 
        const NumericMatrix density,
        const IntegerVector PIA0,
        const NumericMatrix Tsk,
        const NumericVector pmix, 
        uint64_t base_seed,
        
        NumericMatrix chatmat,
        NumericVector chatout)
        : 
        mm(mm), nmark(nmark), cc0(cc0), grain(grain), nsim(nsim), sightmodel(sightmodel),
        sumD(sumD), area(area), distrib(distrib), anytelem(anytelem),
        binomN(binomN), markocc(markocc), pID(pID), group(group), gk0(gk0), hk0(hk0), 
        density(density),PIA0(PIA0), Tsk(Tsk), pmix(pmix),
        chatmat(chatmat), chatout(chatout), base_seed(base_seed) {
        
        // now can initialise these derived counts
        kk = Tsk.nrow();             // number of detectors
        ss = Tsk.ncol();             // number of occasions
        nmix = pmix.size();
        k1 = kk - anytelem;          // number of detectors excl. telemetry

        allsighting = (sightmodel >= 5);
        if (nsim < 2)
            Rcpp::stop ("nsim for chat must be at least 2, and preferably much more!");
        if (sightmodel == 1) { // conditional likelihood incompatible with unresolved sightings 
            for (int s=0; s<ss; s++) {
                if (markocc[s] < 0) resultcode = 2;
            }
        }
        
        Nm = getNm();
        a0 = geta0();  // INTERIM
        sumNm = std::accumulate(Nm.begin(), Nm.end(), 0.0);
        if (distrib) cumprob = getcumprob();  // only for fixed N multinomial 

        RNGScope scope;        
    }
    //==============================================================================
    
    std::vector<double> getNm() {
        std::vector<double> Nm(mm);
        for (int m=0; m<mm; m++) {
            Nm[m] = density[m] * sumD * area;
        }
        return Nm;
    }
    //==============================================================================
    
    std::vector<double> geta0() {
        std::vector<double> a0(nmix);
        for (int x=0; x<nmix; x++) {
            a0[x] = 0;
        }
        return a0;
    }
    //==============================================================================
    
    // cell membership fixed N
    std::vector<double> getcumprob() {
        std::vector<double> cumprob(mm);
        cumprob[0] = Nm[0] / sumNm;
        for (int m=1; m<mm; m++) {
            cumprob[m] = cumprob[m-1] + Nm[m] / sumNm;
        }
        return cumprob;
    }
    //==============================================================================
    
    int locate(double x, int mm, std::vector<double> cumprob){
        // more efficient search Press et al. 1989 p 101 
        // assume cumprob in ascending order 
        int ju, jm, jl;
        jl = 0;
        ju = mm;
        while (ju-jl > 1) {
            jm = (ju+jl) / 2;
            if(x > cumprob[jm]) 
                jl = jm;
            else 
                ju = jm;
            // Rprintf("x %8.6f cumprob[jm] %8.6f jl %4d jm %4d ju %4d \n", x, cumprob[jm], jl, jm, ju); 
        }
        return(jl);
    }
    //------------------------------------------------------------------------------
    
    //probability animal unmarked at s,m becomes marked at s
    double getpmark(int x, int s, int m) {
        int c,k;
        int n = 0;
        double pp = 1.0;
        double H = 0;
        if (markocc[s] > 0) {  // marking occasions only 
            for (k=0; k< k1; k++) {
                c = PIA0[i4(n,s,k,x,nmark,ss,kk)] - 1;
                if (c >= 0) {    // drops unset traps 
                    // pID always 1.0 on marking occasions
                    if (binomN[s] == -2) {
                        // accumulate hazard over traps
                        H += Tsk(k,s) * hk0[i3(c, k, m, cc0, kk)];
                    }
                    else {   
                        pp *= pski(binomN[s], 0, Tsk(k,s), gk0[i3(c, k, m, cc0, kk)], 1.0);  
                    }
                }
            }
            if (binomN[s] == -2) {
                pp = exp(-H);
            }
        }
        return (1-pp);    // Pr(detected and marked) at this occasion 
    }
    //------------------------------------------------------------------------------
    
    std::vector<double> onesim (int r, std::mt19937_64& gen) {
        
        int c, i, k, m, s, x;
        int np = 0;
        int N, Nmarked, temppop;
        double mu1, mu2, pmark, Hskx;
        
        std::vector<double> cumprobmkd(mm);    // marked cell membership
        std::vector<double> cumprobunmkd(mm);  // unmarked cell membership
        if (grain<1) {
            Rprintf("Starting replicate r %4d \n", r);
        }
        std::vector<int> initialpopunmarked(mm,0);  // initial population in each cell (unmarked) 
        std::vector<int> initialpopmarked(mm,0);    // initial population in each cell (allsighting)
        std::vector<int> popunmarked(mm*nmix);  // unmarked population in each cell 
        std::vector<int> popmarked(mm*nmix);    // marked population in each cell 
        std::vector<int> xi(3,0);
        std::vector<double> sump(3,0);
        std::vector<double> musk(3,0);
        std::vector<double> p(3,0);
        std::vector<double> out(8);
        double A = mm * area;              // total mask area 
        
        // instantiate uniform PRNG
        std::uniform_real_distribution<double> u_dist(0.0, 1.0);
        
        //-----------------------------------------------------
        
        // sighting-only, conditioning on Nm 
        if (allsighting && sightmodel != 7) {
            cumprobmkd[0] = density[0];
            cumprobunmkd[0] = (Nm[0] - nmark * density[0]) / (sumNm - nmark);
            for (m=1; m<mm; m++) {
                cumprobmkd[m] = cumprobmkd[m-1] + density[m];
                cumprobunmkd[m] = cumprobunmkd[m-1] + (Nm[m] - nmark * density[m]) / (sumNm - nmark) ;
            }
            // pre-set marked population if conditioning on nmark 
            for (i=0; i<nmark; i++) {
                m = locate(u_dist(gen), mm, cumprobmkd);  // 2026-09-17 u_dist(gen) replaces unif_rand()
                initialpopmarked[m] += 1;
            }
        }
        
        //------------------------------------------------------------------
        // Random initial population 
        // for telemetry marking sightmodel=7 we want the entire population
        
        if (distrib) {   // N 'fixed' (= binomial n);  multinomial cells 
            // simulated density is matched by discreteN only in the long-run average 
            // N = discreteN(sumNm);  
            N = sumNm;
            int tn = (int) N;
            if (N != tn) N = (tn + (u_dist(gen) < (N-tn)));
            if (allsighting && sightmodel != 7) {
                if (N>nmark) {
                    for (i=0; i<(N-nmark); i++) {
                        m = locate(u_dist(gen), mm, cumprobunmkd);
                        initialpopunmarked[m] += 1;    // unmarked animals in cell m 
                    }
                }
            }
            else {
                for (i=0; i<N; i++) {
                    m = locate(u_dist(gen), mm, cumprob);
                    if ((m<0) || (m>=mm)) {
                        Rprintf("erroneous location of simulated animal in sightingchat\n");
                        resultcode = 1;
                    }
                    initialpopunmarked[m] += 1;
                }
            }
        }
        else {   // N Poisson, cells multinomial allowing for marked 
            if (allsighting && sightmodel != 7) {
                std::poisson_distribution<int> p_dist(sumNm);
                N = p_dist(gen); // R::rpois(sumNm);
                if (N>nmark) {
                    for (i=0; i<(N-nmark); i++) {
                        m = locate(u_dist(gen), mm, cumprobunmkd);
                        initialpopunmarked[m] += 1;
                    }
                }
            }
            else {     // Poisson total, Poisson each cell 
                for (m=0; m < mm; m++) {
                    std::poisson_distribution<int> p_dist(Nm[m]);
                    initialpopunmarked[m] = p_dist(gen);  // R::rpois(Nm[m]);
                }
            }
        }
        //------------------------------------------------------------------
        
        // initial pop same for all mixture classes
        for (x=0; x<nmix; x++) {
            for (m=0;m<mm;m++) {
                if (allsighting) {
                    popunmarked[x*mm+m] = initialpopunmarked[m];
                    popmarked[x*mm+m]   = initialpopmarked[m];
                }
                else {
                    popunmarked[x*mm+m] = initialpopunmarked[m];
                    popmarked[x*mm+m]   = 0;
                }
            }
        }
        out[7] = std::accumulate(initialpopunmarked.begin(), initialpopunmarked.end(), 0.0);
        if (grain<1) {
            Rprintf("popn %8.6e \n", out[7]);
        }
        
        mu1 = 0; mu2 = 0;
        for (s=0; s < ss; s++) {
            if (grain<1) {
                // Rprintf("Starting occasion s %4d \n", s);
            }
            //-----------------------------------------------------------------------
            // marking occasions 
            // update marked and unmarked populations 
            // do not update when telemetrytype = "marking" (sightmodel 7)
            // when all animals remain nominally unmarked
            if (markocc[s]==1) { 
                if (sightmodel != 7) {
                    // update marked and unmarked populations 
                    for (x=0; x<nmix; x++) {
                        Nmarked=0;                   // for check
                        for (m=0; m < mm; m++) {
                            temppop = popunmarked[x*mm + m];
                            if (temppop > 0) {
                                pmark = getpmark(x, s, m); 
                                if (grain<1) Rprintf ("s %4d m %4d pmark %8.6e\n", s, m, pmark);
                                for (i=0; i<temppop; i++) {
                                    if (u_dist(gen) < pmark) {
                                        popunmarked[x * mm + m] --;
                                        popmarked[x * mm + m] ++;
                                        Nmarked++;   // for check
                                    }
                                }
                            }
                        }
                        if (grain<1) Rprintf ("s %4d x %3d Nmarked %8d\n", s, x, Nmarked);
                    }
                }
            }    
            //-----------------------------------------------------------------------
            
            // sighting occasions 
            // accumulate sightings 
            else {                         
                for (k=0; k < k1; k++) {
                    if (grain<1) {
                        // popn = std::accumulate(popunmarked.begin(), popunmarked.end(), 0.0);
                        // Rprintf("Starting occasion s %4d detector k %4d nmark %4d popunmarked %8d\n", 
                        //         s, k, nmark, popn);
                    }
                    for (i=0; i<3; i++) musk[i] = 0;
                    for (x=0; x<nmix; x++) {
                        mu1 = 0; mu2 = 0; 
                        // n = 0, generic individual for this mixture class x
                        c = PIA0[i4(0,s,k,x,nmark,ss,kk)] - 1;
                        if (c >= 0) {                    // drops unset traps 
                            for (m=0; m < mm; m++) {
                                Hskx = Tsk(k,s) * hk0[i3(c,k,m,cc0,kk)];	
                                
                                if (popunmarked[x*mm+m]>0) {      // any unmarked animals at m 
                                    if (sightmodel == 1) {        // CL 
                                    }
                                    else if (sightmodel >= 5) {   // sightmodel 5,6,7 all sighting
                                        mu1 += popunmarked[x*mm+m] * Hskx;
                                        if (markocc[s] == -1)
                                            mu1 += popmarked[x*mm+m] * Hskx;
                                    }
                                    else {                       // sightmodel  0,2,3,4 
                                        mu1 += popunmarked[x*mm+m] * Hskx ;
                                        if (markocc[s] == -1)
                                            mu1 += popmarked[x*mm+m] * Hskx;
                                    }
                                }    // end popunmarked
                                
                                if (popmarked[x*mm+m]>0) {       // any marked animals at m 
                                    if (sightmodel == 1) {        // CL 
                                        if (markocc[s] == 0) 
                                            mu2 +=  popmarked[x*mm+m] * Hskx;
                                    }
                                    else if (sightmodel == 5) {   // sightmodel 5 all sighting, known 
                                        if (markocc[s] == 0) 
                                            mu2 +=  density[m] * nmark * Hskx; 
                                    }
                                    else if (sightmodel == 6) {   // sightmodel 6 all sighting, unknown 
                                        if (markocc[s] == 0) 
                                            if (a0[x] > 0)   // 2015-12-31
                                                mu2 +=  density[m] * nmark * (A/a0[x]) * Hskx; 
                                    }
                                    else {                  // sightmodel  0,2,3,4 
                                        if (markocc[s] == 0) {
                                            mu2 += popmarked[x*mm+m] * Hskx;
                                        }
                                    }
                                } // end popmarked
                                
                            }   // over mask points
                        }  // if detector used
                        
                        // mu1, mu2 expected Tu, Tm for given m,x,s,k
                        musk[0] += mu1 * pmix[x];
                        musk[1] += mu2 * (1 - pID[s + ss*x]) * pmix[x];
                        musk[2] = musk[0] + musk[1];
                        
                    }  // end loop over latent classes 
                    
                    // simulate actual counts and increment xi
                    if (binomN[s] < 0) {   // multi -2, proximity -1
                        np += 1;
                        for (i=0; i<3; i++) {
                            p[i] = 1-exp(-musk[i]);
                            sump[i] += p[i];
                            if (musk[i]>tol) xi[i] += (u_dist(gen) < p[i]);
                        }
                    }
                    else if (binomN[s] == 0) { // Poisson cell counts
                        for (i=0; i<3; i++) {
                            if (musk[i]>tol) {
                                std::poisson_distribution<int> p_dist(musk[i]);
                                xi[i] += p_dist(gen); 
                            }
                        }
                    }  
                    else {
                        Rcpp::stop ("binomial sightings are not currently allowed");
                    }
                    if (grain<1) {
                        Rprintf("r %4d s %4d k %4d musk[0] %8.6e musk[1] %8.6e \n",
                                r, s, k, musk[0], musk[1]);
                    }
                }  // end loop over detectors  
            }  // sighting occasions
        }  // end loop over occasions 
        
        // copy results xi, sump to output vector
        for (i=0; i<3; i++) {
            if (i==0 || sightmodel<7) {
                out[i] = xi[i];
                out[i+3] = sump[i];
            }
            else {
                out[i] = 0;
                out[i+3] = 0;
            }
        }
        out[6] = np;
        return(out);
    }
    //==============================================================================
    
    void chatvar () {
        std::vector<double> varx(3, 0.0);
        std::vector<double> expectedvar(3);
        std::vector<double> meanp(3, 0.0);
        std::vector<double> meanx(3, 0.0);
        
        std::vector<double> xi(3);
        std::vector<double> sump(3);
        
        int np;
        double delta;
        
        np = round(chatmat(0,6));    // same for all replicates
        for (int r=0; r<nsim; r++) {
            for (int i=0; i<3; i++) {
                xi[i] = chatmat(r,i);
                sump[i] = chatmat(r,i+3);
                if (grain<1) {
                    if (i==0)
                        Rprintf("r %4d xi[i] %8.6e \n", r, xi[i]);
                }
                // Welford's algorithm for incremental mean and variance
                delta = xi[i] - meanx[i];
                meanx[i] += delta/(r+1);
                varx[i] += delta*(xi[i] - meanx[i]);
                // assuming uniform probability across cells ? 
                if (np>0) {
                    delta = sump[i]/np - meanp[i];
                    meanp[i] += delta/(r+1);
                }
            }
        }
        
        // relate simulated variance to expected variance
        for (int i=0; i<3; i++) {
            varx[i] /= nsim-1;
            if (np>0)   // multi, proximity; implicitly assume all occasions same 
                expectedvar[i] = np * meanp[i] * (1-meanp[i]);    
            else        // Poisson 
                expectedvar[i] = meanx[i];
            
            if (grain<1) {
                Rprintf("i %4d np %4d meanp[i] %8.6e varx[i] %8.6e expectedvar[i] %8.6e \n", 
                        i, np, meanp[i], varx[i], expectedvar[i]);
            }
            if (expectedvar[i] > 0) 
                chatout[i] = varx[i]/expectedvar[i]; 
            else 
                chatout[i] = 1;
        }
    }
    
    
    //==============================================================================
    
    // function call operator that works for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {    
        std::vector<double> chatvec;
        // Instantiate thread-local PRNG with thread-unique offset seed
        std::mt19937_64 gen(base_seed + begin);
        
        for (std::size_t r = begin; r < end; r++) {
            chatvec = onesim (r, gen);
            for (int i=0; i<8; i++) chatmat(r,i) = chatvec[i];
        }
    }
    //==============================================================================
};

// [[Rcpp::export]]
List sightingchatcpp (
        const int           mm, 
        const int           nc,          // number marked (number of rows in capthist incl zeroes)
        const int           cc0, 
        const int           grain, 
        const int           ncores, 
        const int           nsim,        // number of replicate simulations for chat 
        const int           sightmodel,  // 5 allsighting known n0, 6 allsighting unknown n0
        const double        sumD,
        const double        area,
        const int           distrib,
        const int           anytelem,
        const IntegerVector binomN,      // detector -2 multi, -1 proximity 0 Poisson count 1 Binomial from usage, 2...etc. 
        const IntegerVector markocc, 
        const NumericMatrix pID, 
        const IntegerVector group,       // group number for 0<=n<*nc   [full likelihood only] 
        const NumericVector gk0, 
        const NumericVector hk0, 
        const NumericMatrix density,     // relative density - sums to 1.0
        const IntegerVector PIA0, 
        const NumericMatrix Tsk,         // nk x s usage matrix 
        const NumericVector pmix,
        const uint64_t      seed,
        const bool          verbose) {
    
    NumericMatrix chatmat(nsim,8); 
    NumericVector chatout(3); 
    
    // Construct and initialise
    chat somechat (mm, nc, cc0, grain, nsim, sightmodel, sumD, area, distrib, 
                   anytelem, binomN, markocc, pID, 
                   group, gk0, hk0, density, PIA0, Tsk, pmix, seed, chatmat, chatout);
    
    if (ncores>1) {
        // Run operator() on multiple threads
        parallelFor(0, nsim, somechat, grain, ncores);
    }
    else {
        // for debugging avoid multithreading and allow R calls e.g. Rprintf
        somechat.operator()(0,nsim);    
    }
    
    //-----------------------------------------------------
    
    // Rprintf("sims completed\n");
        
    somechat.chatvar();
    
    // Return consolidated result
    if (verbose) {
        return (List::create(Named("resultcode") = 0, 
                             Named("chat") = chatout, 
                             Named("chatmat") = chatmat));
    }
    else {
        return (List::create(Named("resultcode") = 0, Named("chat") = chatout));
    }
    
}
//==============================================================================
