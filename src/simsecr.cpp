// Simulate capture histories from fitted model
// ******************* indicates weakness

#include <Rcpp.h>
#include "poly.h"
using namespace Rcpp;
//==============================================================================

// 2024-08-02 changed output order of 'value' to i,s,k
//            not clear whether this requires new sortkey for poly, signal

NumericVector getpar (
        const int i, 
        const int s, 
        const int k, 
        const int xi, 
        const int N, 
        const int ss, 
        const int nk, 
        const int fn, 
        const bool bswitch, 
        const IntegerVector &PIA0, 
        const NumericMatrix &gsb0val,
        const IntegerVector &PIA1, 
        const NumericMatrix &gsb1val)
{
    // bswitch determines whether to use naive (0) or caught before (1) 
    int wxi;
    int c;
    wxi = i4(i,s,k,xi,N,ss,nk);
    NumericVector gsb(3);      // defaults to zero
    
    if (bswitch) {
        c = PIA0[wxi]-1;
        gsb[0] = gsb0val(c,0);
        gsb[1] = gsb0val(c,1);
        if (par3(fn)) gsb[2] = gsb0val(c,2);
    }
    else {
        c = PIA1[wxi]-1;
        gsb[0] = gsb1val(c,0);
        gsb[1] = gsb1val(c,1);
        if (par3(fn)) gsb[2] = gsb1val(c,2);
    }
    return(gsb);
}
//==============================================================================

int rdiscrete (
        const int n, 
        const NumericVector &pmix)
    // return random discrete observation from distribution in pmix 
{
    std::vector<double> cumpmix(n+1);
    int x;
    double r;
    if (n<1) Rcpp::stop ("invalid n in rdiscrete");
    if (n==1) return (0);
    else {
        cumpmix[0] = 0;
        for (x=0; x<n; x++) {
            cumpmix[x+1] = cumpmix[x] + pmix[x];
        }
        r = unif_rand();
        for (x=1; x<=n; x++) if (r<cumpmix[x]) break;
        return(x);
    }
}
//==============================================================================

// btype is code for behavioural response 
// 0 none
// 1 individual
// 2 individual, trap-specific
// 3 trap-specific

int bswitch (
        const int btype, 
        const int N, 
        const int i, 
        const int k, 
        const std::vector<int> &caughtbefore)
{
    if (btype == 0)
        return(0);
    else if (btype == 1) 
        return(caughtbefore[i]);
    else if (btype == 2) 
        return(caughtbefore[k * (N-1) + i]);
    else if (btype == 3) 
        return(caughtbefore[k]);
    else 
        Rcpp::stop("unrecognised btype in simsecr");
    return(0);
}
//==============================================================================


// [[Rcpp::export]]
List simdetectpointcpp (
        const int           &detect,      // detector -1 single, 0 multi, 1 proximity, 2 count,... 
        const int           &N, 
        const int           &cc0,
        const int           &cc,
        const NumericVector &gk0, 
        const NumericVector &gk, 
        const NumericVector &hk0, 
        const NumericVector &hk, 
        const IntegerVector &PIA0,       // lookup which g0/sigma/b combination to use for given g, S, K [naive animal] 
        const IntegerVector &PIA1,        // lookup which g0/sigma/b combination to use for given n, S, K  [caught before] 
        const int           &nmix,        // number of classes 
        const IntegerVector &knownclass, // known membership of 'latent' classes 
        const NumericVector &pmix,       // membership probabilities
        const NumericMatrix &Tsk,        // ss x kk array of 0/1 usage codes or effort 
        const int           &btype,       // code for behavioural response  0 none etc. 
        const int           &Markov,      // learned vs transient behavioural response 0 learned 1 Markov 
        const IntegerVector &binomN     // number of trials for 'count' detector modelled with binomial 
)
{
    //  detect may take values -
    // -1  single-catch traps
    //  0  multi-catch traps
    //  1  binary proximity detectors
    //  2  count  proximity detectors
    
    int    kk = Tsk.nrow();            // number of detectors 
    int    ss = Tsk.ncol();            // number of occasions
    
    double p;
    int    i,k,s;
    int    ik;
    int    nc = 0;
    int    count = 0;
    double runif;
    int    wxi = 0;
    int    c = 0;
    double Tski = 1.0;  
    bool   before;
    
    std::vector<int> caughtbefore(N * kk, 0);
    std::vector<int> x(N, 0);          // mixture class of animal i 
    
    // return values
    IntegerVector caught(N);           // caught in session 
    IntegerVector value (N*ss*kk);     // return value array
    
    //========================================================
    // 'single-catch only' declarations 
    int    nanimals;
    int    ntraps;
    int    tr_an_indx = 0;
    int    anum = 0;
    int    tnum = 0;
    int    nextcombo;
    int    finished;
    int    OK;
    double lambda;
    double event_time;
    std::vector<int> occupied(kk);
    std::vector<double> intrap(N);
    std::vector<trap_animal> tran(N * kk);
    
    //========================================================
    // 'multi-catch only' declarations 
    std::vector<double> h(N * kk);        // multi-catch only 
    std::vector<double> hsum(N);          // multi-catch only 
    std::vector<double> cump(kk+1,0);     // multi-catch only 
    
    //========================================================
    // MAIN LINE 
    
    List nullresult = List::create(Named("n") = 0,
                                   Named("caught") = caught,
                                   Named("value") = value,
                                   Named("resultcode") = 2);
    
    RNGScope scope;             // Rcpp initialise and finalise random seed 
    
    if ((detect < -1) || (detect > 2)) {
        return(nullresult);
    }
    
    //----------------------------------------------------------------------------
    // mixture models 
    if (nmix>1) {
        if (nmix>2)
            Rcpp::stop("simsecr nmix>2 not implemented");
        for (i=0; i<N; i++) {
            if (knownclass[i] > 1) 
                x[i] = knownclass[i] - 2;      // knownclass=2 maps to x=0 etc. 
            else
                x[i] = rdiscrete(nmix, pmix) - 1;
        }
    }
    
    // ------------------------------------------------------------------------- 
    // MAIN LOOP 
    
    for (s=0; s<ss; s++) {
        
        // --------------------------------------------------------------------- 
        // single-catch traps 
        if (detect == -1) {
            // initialise day 
            tr_an_indx = 0;
            nanimals = N;
            ntraps   = kk;
            for (i=0; i<N; i++) intrap[i] = 0;
            for (k=0; k<kk; k++) occupied[k] = 0;
            nextcombo = 0;
            
            // make tran 
            for (i=0; i<N; i++) {  // animals randomtime
                for (k=0; k<kk; k++) { // traps 
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        before = bswitch (btype, N, i, k, caughtbefore);
                        wxi =  i4(i, s, k, x[i], N, ss, kk);
                        if (before)
                            c = PIA1[wxi] - 1;
                        else 
                            c = PIA0[wxi] - 1;
                        if (c >= 0) {    // ignore unused detectors 
                            if (before)
                                lambda = hk[i3(c, k, i, cc, kk)];
                            else
                                lambda = hk0[i3(c, k, i, cc0, kk)];
                            lambda = lambda * Tski;
                            event_time = randomtimel(lambda);
                            if (event_time <= 1) {
                                tran[tr_an_indx].time   = event_time;
                                tran[tr_an_indx].animal = i;    // 0..N-1 
                                tran[tr_an_indx].trap   = k;    // 0..kk-1 
                                tr_an_indx++;
                            }
                        }
                    }
                }
            }
            // end of make tran 
            
            if (tr_an_indx > 1) probsort (tr_an_indx, tran);
            
            while ((nextcombo < tr_an_indx) && (nanimals>0) && (ntraps>0)) {
                finished = 0;
                OK       = 0;
                while ((1-finished)*(1-OK) > 0) {      // until finished or OK 
                    if (nextcombo >= (tr_an_indx))
                        finished = 1;                  // no more to process 
                    else {
                        anum = tran[nextcombo].animal;
                        tnum = tran[nextcombo].trap;
                        OK = (1-occupied[tnum]) * (1-intrap[anum]); // not occupied and not intrap 
                        nextcombo++;
                    }
                }
                if (finished==0) {
                    // Record this capture 
                    occupied[tnum] = 1;
                    intrap[anum]   = tnum+1;         // trap = k+1 
                    nanimals--;
                    ntraps--;
                }
            }
            
            for (i=0; i<N; i++) {
                if (intrap[i]>0) {
                    if (caught[i]==0) {                    // first capture of this animal 
                        nc++;
                        caught[i] = nc;                    // nc-th animal to be captured 
                    }
                    value[i3(i, s, intrap[i]-1, N, ss)] = 1;  
                }
            }
        }
        
        // -------------------------------------------------------------------------- 
        // multi-catch trap; only one site per occasion 
        else if (detect == 0) {
            for (i=0; i<N; i++) {
                hsum[i] = 0;
                for (k=0; k<kk; k++) {
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        before = bswitch (btype, N, i, k, caughtbefore);
                        wxi =  i4(i, s, k, x[i], N, ss, kk);
                        if (before)
                            c = PIA1[wxi] - 1;
                        else 
                            c = PIA0[wxi] - 1;
                        if (c >= 0) {    // ignore unused detectors 
                            if (before)
                                h[k * N + i] = Tski * hk[i3(c, k, i, cc, kk)];
                            else
                                h[k * N + i] = Tski * hk0[i3(c, k, i, cc0, kk)];
                            hsum[i] += h[k * N + i];
                        }
                    }
                }
                
                for (k=0; k<kk; k++) {
                    cump[k+1] = cump[k] + h[k * N + i]/hsum[i];
                }
                if (unif_rand() < (1-exp(-hsum[i]))) {
                    if (caught[i]==0)  {        // first capture of this animal 
                        nc++;
                        caught[i] = nc;
                    }
                    // find trap with probability proportional to p
                    // searches cumulative distribution of p  
                    runif = unif_rand();
                    k = 0;
                    // while ((runif > cump[k]) && (k<kk)) k++;   // bug fix 2019-10-04
                    while ((runif > cump[k+1]) && (k<kk)) k++;
                    value[i3(i, s, k, N, ss)] = 1;  
                }
            }
        }
        
        // -------------------------------------------------------------------------------- 
        // the 'proximity' group of detectors 1:2 - proximity, count 
        else if ((detect >= 1) && (detect <= 2)) {
            for (i=0; i<N; i++) {
                for (k=0; k<kk; k++) {
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        before = bswitch (btype, N, i, k, caughtbefore);
                        wxi =  i4(i, s, k, x[i], N, ss, kk);
                        if (before)
                            c = PIA1[wxi] - 1;
                        else 
                            c = PIA0[wxi] - 1;
                        if (c >= 0) {    // ignore unused detectors 
                            if (before)
                                p = gk[i3(c, k, i, cc, kk)];
                            else
                                p = gk0[i3(c, k, i, cc0, kk)];
                            if (p < -0.1) { 
                                return(nullresult);
                            }  
                            if (p>0) {
                                if (detect == 1) {
                                    if (fabs(Tski-1) > 1e-10)
                                        p = 1 - pow(1-p, Tski);
                                    count = unif_rand() < p;           // binary proximity 
                                }
                                else if (detect == 2) {             // count proximity 
                                    if (binomN[s] == 1)
                                        count = rcount(round(Tski), p, 1);
                                    else
                                        count = rcount(binomN[s], p, Tski);
                                }
                                if (count>0) {
                                    if (caught[i]==0) {              // first capture of this animal 
                                        nc++;
                                        caught[i] = nc;
                                    }
                                    value[i3(i, s, k, N, ss)] = count;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if ((btype > 0) && (s < (ss-1))) {
            // update record of 'previous-capture' status 
            if (btype == 1) {
                for (i=0; i<N; i++) {
                    if (Markov) 
                        caughtbefore[i] = 0;
                    for (k=0; k<kk; k++)
                        caughtbefore[i] = R::imax2 (value[i3(i, s, k, N, ss)], caughtbefore[i]);
                }
            }
            else if (btype == 2) {
                for (i=0; i<N; i++) {
                    for (k=0; k<kk; k++) {
                        ik = k * (N-1) + i;
                        if (Markov) 
                            caughtbefore[ik] = value[i3(i, s, k, N, ss)];
                        else 
                            caughtbefore[ik] = R::imax2 (value[i3(i, s, k, N, ss)], 
                                                         caughtbefore[ik]);
                    }
                }
            }
            else {
                for (k=0;k<kk;k++) {
                    if (Markov) 
                        caughtbefore[k] = 0;
                    for (i=0; i<N; i++) 
                        caughtbefore[k] = R::imax2 (value[i3(i, s, k, N, ss)], caughtbefore[k]);
                }
            }
        }
        
    }   // loop over s 
    
    return (List::create(Named("n") = nc, 
                         Named("caught") = caught,
                         Named("value") = value,
                         Named("resultcode") = 0));
    
}
//==============================================================================

// [[Rcpp::export]]
List simdetectpolycpp (
        const int           detect,      // detector -1 single, 0 multi, 1 proximity, 2 count,... 
        const int           fn,          // code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 = uniform 
        const int           nmix,        // number of classes 
        const int           btype,       // code for behavioural response  0 none etc. 
        const int           Markov,      // learned vs transient behavioural response 0 learned 1 Markov 
        const IntegerVector &kk,         // number of vertices per polygon (zero-terminated vector)
        const NumericMatrix &animals,    // x,y points of animal range centres (first x, then y) 
        const NumericMatrix &traps,      // x,y locations of traps (first x, then y) 
        const NumericMatrix &gsb0val,    // Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal]
        const NumericMatrix &gsb1val,    // Parameter values (matrix nr= comb of g0,sigma,b nc=3) [caught before]
        const IntegerVector &PIA0,       // lookup which g0/sigma/b combination to use for given g, S, K [naive animal] 
        const IntegerVector &PIA1,        // lookup which g0/sigma/b combination to use for given n, S, K  [caught before] 
        const IntegerVector &knownclass, // known membership of 'latent' classes 
        const NumericVector &pmix,       // membership probabilities
        const NumericMatrix &Tsk,        // ss x kk array of 0/1 usage codes or effort 
        const IntegerVector &binomN,     // number of trials for 'count' detector modelled with binomial 
        const int           maxperpoly  
)
{
    //  detect may take values -
    //  3  exclusive polygon detectors
    //  4  exclusive transect detectors
    //  6  polygon detectors
    //  7  transect detectors
    
    int    ss = Tsk.ncol();            // number of animals 
    int    nk = Tsk.nrow();            // number of polygons/transects
    int    N = animals.nrow();         // number of occasions 
    if (nk!=(kk.size()-1)) Rcpp::stop("Tsk doesn't match kk");
    int    i,j,k,l,s;
    int    ik;
    int    nc = 0;
    int    nd = 0;
    int    count = 0;
    double Tski = 1.0;
    bool   before;
    
    // work arrays
    double *work = NULL;
    int    *sortorder = NULL;
    double *sortkey = NULL;
    
    std::vector<int> caughtbefore(N * nk, 0);
    std::vector<int> x(N, 0);          // mixture class of animal i 
    
    //========================================================
    std::vector<int> cumk(nk+1); 
    int    sumk = 0;           // total number of vertices 
    bool   gotcha = false;
    int    n1,n2;
    NumericVector xy(2);
    NumericVector gsb(3);
    NumericMatrix gsbval(1,3);
    double w, ws;
    double dx,dy,d;
    int    maxdet;
    double *cumd = NULL;
    rpoint *line = NULL;
    rpoint xyp;
    rpoint animal;
    double lx;
    double maxg = 0;
    double lambdak;  // temp value for Poisson rate 
    double grx;      // temp value for integral gr 
    double H;
    int    J;
    int    maybecaught;
    double pks;
    double sumhaz;
    // double *ex;
    
    //========================================================
    // MAIN LINE 

    RNGScope scope;             // Rcpp initialise and finalise random seed 
    
    cumk[0] = 0;
    for (i=0; i<=nk; i++) {                               
        if (kk[i]<=0) break;
        cumk[i+1] = cumk[i] + kk[i];
    }
    sumk = cumk[nk];
    if ((detect == 6) || (detect == 7))
        maxdet = N * ss * nk * maxperpoly;
    else
        maxdet = N * ss * nk;
    
    // return values
    IntegerVector caught(N);           // caught in session 
    NumericMatrix detectedXY(maxdet, 2); // x,y locations of detections  
    IntegerVector value (N*ss*nk);     // return value array of detections
    
    List nullresult = List::create(Named("n") = 0,
                                   Named("caught") = caught,
                                   Named("detectedXY") = detectedXY,
                                   Named("value") = value,
                                   Named("resultcode") = 2);
    
    if ((detect < 3) || (detect > 7) || (detect==5)) {
        return(nullresult);
    }
    
    //----------------------------------------------------------------------------
    
    if ((detect == 4) || (detect == 7)) {                         // transect only 
        line = (rpoint *) R_alloc(sumk, sizeof(struct rpoint));
        cumd = (double *) R_alloc(sumk+1, sizeof(double));
        // coordinates of vertices 
        for (i=0; i<sumk; i++) {
            line[i].x = traps(i,0);
            line[i].y = traps(i,1);
        }
        // cumulative distance along line; all transects end on end 
        for (k=0; k<nk; k++) {
            cumd[cumk[k]] = 0;
            for (i=cumk[k]; i<(cumk[k+1]-1); i++) {
                cumd[i+1] = cumd[i] + distance1(line[i], line[i+1]);
            }
        }
    }
    //----------------------------------------------------------------------------
    
    work = (double*) R_alloc(maxdet, sizeof(double));  
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));
    
    //----------------------------------------------------------------------------
    // mixture models 
    if (nmix>1) {
        if (nmix>2)
            Rcpp::stop("simsecr nmix>2 not implemented");
        for (i=0; i<N; i++) {
            if (knownclass[i] > 1) 
                x[i] = knownclass[i] - 2;      // knownclass=2 maps to x=0 etc. 
            else
                x[i] = rdiscrete(nmix, pmix) - 1;
        }
    }
    
    // find maximum distance between animal and detector vertex 
    if (detect==3 || detect == 6) { 
        w = 0;
        J = cumk[nk];
        for (i = 0; i< N; i++) {
            for (j = 0; j < J; j++) {
                dx = animals(i,0) - traps(j,0);
                dy = animals(i,1) - traps(j,1);
                d = std::sqrt(dx*dx + dy*dy);
                if (d > w) w = d;
            }
        } 
    }
    else w = 1;
    
    // ------------------------------------------------------------------------- 
    // MAIN LOOP 
    
    for (s=0; s<ss; s++) {
        
        // -------------------------------------------------------------------------------- 
        // exclusive polygon detectors  
        if (detect == 3) {
            
            for (i=0; i<N; i++) {
                // this implementation assumes NO VARIATION AMONG DETECTORS 
                before = bswitch (btype, N, i, 0, caughtbefore);
                gsb = getpar (i, s, 0, x[i], N, ss, nk, fn, 
                              before, PIA0, gsb0val, PIA1, gsb1val);
                maybecaught = unif_rand() < gsb(0);
                gsb(0) = 1;
                if (w > (10 * gsb(1))) 
                    ws = 10 * gsb(1);
                else 
                    ws = w;
                
                if (maybecaught) {
                    xy = gxy (fn, gsb, ws);                 // simulate one location
                    xy[0] = xy[0] + animals(i,0);
                    xy[1] = xy[1] + animals(i,1);
                    for (k=0; k<nk; k++) {                      // each polygon 
                        Tski = Tsk(k,s);
                        if (fabs(Tski) > 1e-10) {
                            n1 = cumk[k];
                            n2 = cumk[k+1]-1;
                            gotcha = insidecpp(xy, n1, n2, traps);  // assume closed 
                            if (gotcha) {
                                if (caught[i]==0) {             // first capture of this animal 
                                    nc++;
                                    caught[i] = nc;
                                }
                                nd++;
                                value[i3(i, s, k, N, ss)] = 1;  
                                work[(nd-1)*2] = xy[0];
                                work[(nd-1)*2+1] = xy[1];
                                sortkey[nd-1] = (double) (s * N + caught[i]);
                                break;   // no need to look at more poly 
                            }
                        }
                    }
                }
            }
        }
        
        // -------------------------------------------------------------------------------- 
        // exclusive transect detectors  
        else if (detect == 4) {
            // ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
            // conversions for RMatrix input to integralxDNRcpp
            const RcppParallel::RMatrix<double> trapsR(traps);
            const RcppParallel::RMatrix<double> animalsR(animals);
            
            for (i=0; i<N; i++) {                            // each animal 
                animal.x = animals[i];
                animal.y = animals[i + N];
                sumhaz = 0;
                // ------------------------------------ 
                // sum hazard 
                for (k=0; k<nk; k++) {            
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        before = bswitch (btype, N, i, k, caughtbefore);
                        gsb = getpar (i, s, k, x[i], N, ss, nk, fn, 
                                      before, PIA0, gsb0val, PIA1, gsb1val);
                        n1 = cumk[k];
                        n2 = cumk[k+1]-1;
                        H = hintegral1Ncpp(fn, as<std::vector<double>>(gsb));
                        for (i=0;i<3;i++) gsbval(0,i) = gsb(i);
                        const RcppParallel::RMatrix<double> gsbvalR(gsbval);
                        sumhaz += -log(1 - gsb(0) * 
                            integral1DNRcpp (fn, i, 0, gsbvalR, trapsR, animalsR, n1, n2) / H);
                    }
                }
                // ------------------------------------ 
                
                for (k=0; k<nk; k++) {                        // each transect 
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        before = bswitch (btype, N, i, k, caughtbefore);
                        gsb = getpar (i, s, k, x[i], N, ss, nk, fn, 
                                      before, PIA0, gsb0val, PIA1, gsb1val);
                        n1 = cumk[k];
                        n2 = cumk[k+1]-1;
                        H = hintegral1Ncpp(fn, as<std::vector<double>>(gsb));
                        for (i=0;i<3;i++) gsbval(0,i) = gsb(i);
                        const RcppParallel::RMatrix<double> gsbvalR(gsbval);
                        lambdak = gsb(0) * integral1DNRcpp (fn, i, 0, gsbvalR, trapsR, animalsR,
                                      n1, n2) / H;
                        pks = (1 - exp(-sumhaz)) * (-log(1-lambdak)) / sumhaz;
                        count = unif_rand() < pks;
                        maxg = 0;
                        if (count>0) {                       // find maximum - approximate 
                            for (l=0; l<=100; l++) {
                                lx = (cumd[n2] - cumd[n1]) * l/100;
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, gsb, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = R::fmax2(maxg, grx);
                            }
                            for (l=n1; l<=n2; l++) {
                                xyp = line[l];
                                grx = gr (fn, gsb, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = R::fmax2(maxg, grx);
                            }
                            maxg= 1.2 * maxg;                 // safety margin 
                            if (maxg<=0)
                                Rprintf("maxg stop in simsecr\n"); // not found 
                            gotcha = false;
                            l = 0;
                            while (!gotcha) {
                                lx = unif_rand() * (cumd[n2] - cumd[n1]);     // simulate location 
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, gsb, xyp, animal);
                                if (unif_rand() < (grx/maxg))    // rejection sampling 
                                    gotcha = true;
                                l++;
                                if (l % 10000 == 0)
                                    Rcpp::checkUserInterrupt();
                                if (l>1e8) gotcha = true;       // give up and accept anything!!!! 
                            }
                            if (caught[i]==0) {               // first capture of this animal 
                                nc++;
                                caught[i] = nc;
                            }
                            nd++;
                            if (nd >= maxdet) {
                                return(nullresult);  // error
                            }
                            value[i3(i, s, k, N, ss)] = 1;  
                            work[(nd-1)*2] = xyp.x;
                            work[(nd-1)*2+1] = xyp.y;
                            sortkey[nd-1] = (double) (s * N + caught[i]);
                        }
                        if (count>0) break;   // no need to look further 
                    }
                }                                             // end loop over transects 
            }                                                 // end loop over animals 
        }
        
        // -------------------------------------------------------------------------------- 
        // polygon detectors  
        else if (detect == 6) {
            for (i=0; i<N; i++) {
                // this implementation assumes NO VARIATION AMONG DETECTORS 
                before = bswitch (btype, N, i, 0, caughtbefore);
                gsb = getpar (i, s, 0, x[i], N, ss, nk, fn, 
                              before, PIA0, gsb0val, PIA1, gsb1val);
                count = rcount(binomN[s], gsb(0), 1.0); // drop 2019-10-04 Tski);
                gsb(0) = 1;
                if (w > (10 * gsb(1))) 
                    ws = 10 * gsb(1);
                else 
                    ws = w;
                for (j=0; j<count; j++) {
                    xy = gxy (fn, gsb, ws);                  // simulate one location 
                    xy[0] = xy[0] + animals(i,0);
                    xy[1] = xy[1] + animals(i,1);
                    for (k=0; k<nk; k++) {                      // each polygon 
                        Tski = Tsk(k,s);
                        if (fabs(Tski) > 1e-10) {
                            n1 = cumk[k];
                            n2 = cumk[k+1]-1;
                            gotcha = insidecpp(xy, n1, n2, traps);  // assume closed 
                            if (gotcha) {
                                if (caught[i]==0) {             // first capture of this animal 
                                    nc++;
                                    caught[i] = nc;
                                }
                                nd++;
                                if (nd > maxdet) {
                                    return(nullresult);  // error
                                }
                                value[i3(i, s, k, N, ss)]++;
                                work[(nd-1)*2] = xy[0];
                                work[(nd-1)*2+1] = xy[1];
                                sortkey[nd-1] = (double) (k * N * ss + s * N + caught[i]);
                            }
                        }
                    }
                }
            }
        }
        // -------------------------------------------------------------------------------- 
        // transect detectors  
        else if (detect == 7) {
            const RcppParallel::RMatrix<double> trapsR(traps);
            const RcppParallel::RMatrix<double> animalsR(animals);
            
            for (i=0; i<N; i++) {                            // each animal 
                animal.x = animals(i,0);
                animal.y = animals(i,1);
                for (k=0; k<nk; k++) {                        // each transect 
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        before = bswitch (btype, N, i, k, caughtbefore);
                        gsb = getpar (i, s, k, x[i], N, ss, nk, fn, 
                                      before, PIA0, gsb0val, PIA1, gsb1val);
                        n1 = cumk[k];
                        n2 = cumk[k+1]-1;
                        H = hintegral1Ncpp(fn, as<std::vector<double>>(gsb));
                        for (i=0;i<3;i++) gsbval(0,i) = gsb(i);
                        const RcppParallel::RMatrix<double> gsbvalR(gsbval);
                        lambdak = gsb(0) * integral1DNRcpp (fn, i, 0, gsbvalR, trapsR, animalsR,
                                      n1, n2) / H;
                        count = rcount(binomN[s], lambdak, Tski);  // numb detections on transect 
                        maxg = 0;
                        if (count>0) {                       // find maximum - approximate 
                            for (l=0; l<=100; l++) {
                                lx = (cumd[n2]-cumd[n1]) * l/100;
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, gsb, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = R::fmax2(maxg, grx);
                            }
                            for (l=n1; l<=n2; l++) {
                                xyp = line[l];
                                grx = gr (fn, gsb, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = R::fmax2(maxg, grx);
                            }
                            maxg= 1.2 * maxg;                 // safety margin 
                            if (maxg<=0)
                                Rprintf("maxg stop in simsecr\n"); // not found 
                            
                        }
                        for (j=0; j<count; j++) {
                            gotcha = false;
                            l = 0;
                            while (!gotcha) {
                                lx = unif_rand() * (cumd[n2]-cumd[n1]);     // simulate location 
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, gsb, xyp, animal);
                                if (unif_rand() < (grx/maxg))   // rejection sampling 
                                    gotcha = true;
                                l++;
                                if (l % 10000 == 0)
                                    Rcpp::checkUserInterrupt();
                                if (l>1e8) gotcha = true;        // give up and accept anything!!!! 
                            }
                            if (caught[i]==0) {               // first capture of this animal 
                                nc++;
                                caught[i] = nc;
                            }
                            nd++;
                            if (nd >= maxdet) {
                                return(nullresult);  // error
                            }
                            value[i3(i,s,k, N, ss)]++;
                            work[(nd-1)*2] = xyp.x;
                            work[(nd-1)*2+1] = xyp.y;
                            sortkey[nd-1] = (double) (k * N * ss + s * N + caught[i]);
                        }
                    }
                }                                             // end loop over transects 
            }                                                 // end loop over animals 
        }
        
        if ((btype > 0) && (s < (ss-1))) {
            // update record of 'previous-capture' status 
            if (btype == 1) {
                for (i=0; i<N; i++) {
                    if (Markov) 
                        caughtbefore[i] = 0;
                    for (k=0; k<nk; k++)
                        caughtbefore[i] = R::imax2 (value[i3(i, s, k, N, ss)], caughtbefore[i]);
                }
            }
            else if (btype == 2) {
                for (i=0; i<N; i++) {
                    for (k=0; k<nk; k++) {
                        ik = k * (N-1) + i;
                        if (Markov) 
                            caughtbefore[ik] = value[i3(i, s, k, N, ss)];
                        else 
                            caughtbefore[ik] = R::imax2 (value[i3(i, s, k, N, ss)], 
                                                         caughtbefore[ik]);
                    }
                }
            }
            else {
                for (k=0; k<nk;k++) {
                    if (Markov) 
                        caughtbefore[k] = 0;
                    for (i=0; i<N; i++) 
                        caughtbefore[k] = R::imax2 (value[i3(i, s, k, N, ss)], 
                                                    caughtbefore[k]);
                }
            }
        }
        
    }   // loop over s 
    
    for (i=0; i<nd; i++) sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) {
        detectedXY(i,0) = work[sortorder[i]*2];
        detectedXY(i,1) = work[sortorder[i]*2+1];
    }
    
    return (List::create(Named("n") = nc, 
                         Named("caught") = caught,
                         Named("detectedXY") = detectedXY,
                         Named("value") = value,
                         Named("resultcode") = 0));
    
}
//==============================================================================

// [[Rcpp::export]]
List simdetectsignalcpp (
        const int           detect,      // detector -1 single, 0 multi, 1 proximity, 2 count,... 
        const int           nmix,        // number of classes 
        const int           fn,          // code 
        const double        cut,
        const NumericMatrix &gsb0val,    // Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] 
        const IntegerVector &PIA0,       // lookup which g0/sigma/b combination to use for given g, S, K [naive animal] 
        const NumericVector &pmix,
        const IntegerVector &knownclass, // known membership of 'latent' classes 
        const NumericMatrix &animals,    // x,y points of animal range centres (first x, then y) 
        const NumericMatrix &traps,      // x,y locations of traps (first x, then y) 
        const NumericMatrix &dist2,      // distances squared (optional: -1 if unused) 
        const NumericMatrix &Tsk,        // kk x ss array of 0/1 usage codes or effort 
        const NumericVector &miscparm   // detection threshold on transformed scale, etc. 
)
{
    //  detect may take value -
    //  5  signal detectors
    
    int    nk = Tsk.nrow();            // number of traps 
    int    ss = Tsk.ncol();            // number of animals 
    int    N = animals.nrow();         // number of occasions 
    
    int    i,k,s;
    int    nc = 0;
    int    nd = 0;
    double Tski = 1.0;  
    int    maxdet = N * ss * nk;
    
    // work arrays
    double *work = NULL;
    double *noise = NULL;              // detectfn 12,13 only 
    int *sortorder = NULL;
    double *sortkey = NULL;
    NumericVector gsb(3);
    
    IntegerVector x(N, 0);          // mixture class of animal i 
    // return values
    IntegerVector caught(N);        // caught in session 
    NumericVector signal(maxdet);   // return value 
    IntegerVector value (maxdet);   // return value array
    
    //========================================================
    // 'signal-strength only' declarations 
    double muS;
    double muN = 0;
    double sdN = 1;
    double signalvalue;
    double noisevalue;
    //========================================================
    // MAIN LINE 
    
    List nullresult = List::create(Named("n") = 0,
                                   Named("caught") = caught,
                                   Named("signal") = signal,
                                   Named("value") = value,
                                   Named("resultcode") = 2);
    
    RNGScope scope;             // Rcpp initialise and finalise random seed 
    
    if ((detect != 5)) {
        return(nullresult);
    }                                            // signal only 
    if (!((fn == 10) || (fn == 11)))
        Rcpp::stop ("simsecr not implemented for this combination of detector & detectfn");
    work = (double*) R_alloc(maxdet*2, sizeof(double));   /* twice size needed for signal */
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));
    //----------------------------------------------------------------------------
    
    if ((fn==12) || (fn==13)) {
        noise = (double*) R_alloc(maxdet*2, sizeof(double));   /* twice size needed for signal */
    }
    //----------------------------------------------------------------------------
    
    // mixture models 
    if (nmix>1) {
        if (nmix>2)
            Rcpp::stop("simsecr nmix>2 not implemented");
        for (i=0; i<N; i++) {
            if (knownclass[i] > 1) 
                x[i] = knownclass[i] - 2;      // knownclass=2 maps to x=0 etc. 
            else
                x[i] = rdiscrete(nmix, pmix) - 1;
        }
    }
    
    // ------------------------------------------------------------------------- 
    // MAIN LOOP 
    
    // Rprintf("start main loop\n");
    
    for (s=0; s<ss; s++) {
        if ((fn == 12) || (fn == 13)) {
            muN = miscparm[1];
            sdN = miscparm[2];
        }
        for (i=0; i<N; i++) {
            for (k=0; k<nk; k++) {
                Tski = Tsk(k,s);
                if (fabs(Tski) > 1e-10) {
                    // sounds not recaptured 
                    gsb = getpar (i, s, k, x[i], N, ss, nk, fn, 
                                  0, PIA0, gsb0val, PIA0, gsb0val);
                    if ((fn == 10) || (fn == 12))
                        muS  = mufn (i, k, gsb(0), gsb(1), animals, traps, 0);
                    else
                        muS  = mufn (i, k, gsb(0), gsb(1), animals, traps, 1);
                    
                    if ((fn == 10) || (fn == 12))
                        muS  = mufnL (k, i, gsb(0), gsb(1), dist2, 0);
                    else
                        muS  = mufnL (k, i, gsb(0), gsb(1), dist2, 1);
                    
                    signalvalue = norm_rand() * gsb(2) + muS;
                    if ((fn == 10) || (fn == 11)) {
                        if (signalvalue > cut) {
                            if (caught[i]==0) {        // first capture of this animal 
                                nc++;
                                caught[i] = nc;
                            }
                            nd++;
                            value[i3(i, s, k, N, ss)] = 1;
                            work[nd-1] = signalvalue;
                            sortkey[nd-1] = (double) (k * N * ss + s * N + caught[i]);
                        }
                    }
                    else {
                        noisevalue = norm_rand() * sdN + muN;
                        if ((signalvalue - noisevalue) > cut) {
                            if (caught[i]==0) {        // first capture of this animal 
                                nc++;
                                caught[i] = nc;
                            }
                            nd++;
                            value[i3(i,s,k, N,ss)] = 1;
                            work[nd-1] = signalvalue;
                            noise[nd-1] = noisevalue;
                            sortkey[nd-1] = (double) (k * N * ss + s * N + caught[i]);
                        }
                    }
                }
            }
        }
        
    }   // loop over s 
    
    for (i=0; i<nd; i++) sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) signal[i] = work[sortorder[i]];
    if ((fn == 12) || (fn == 13)) {
        for (i=0; i<nd; i++) signal[i+nd] = noise[sortorder[i]];
    }
    
    return (List::create(Named("n") = nc, 
                         Named("caught") = caught,
                         Named("signal") = signal,
                         Named("value") = value,
                         Named("resultcode") = 0));
    
}
//==============================================================================

// expected number of detections
// returns numeric vector with elements of N, ss, kk array of 
// expected detections

// experimental code 2024-07-30

// unfinished:
//   learned responses

double Ey (
        const double p,     // gk[i3(c, k, i, cc, kk)] etc.
        const int binomN, 
        const int detect, 
        const double Tski
)
{
    if (p < -0.1) { 
        return(0.0);
    }  
    if (p > 0) {
        if (detect == 1) {   // binary proximity 
            if (fabs(Tski-1) > 1e-10)   // adjust for usage only if needed
                return(1 - pow(1-p, Tski));
            else
                return(p);
        }
        else if (detect == 2) {       // count proximity 
            if (binomN == 0)                // Poisson
                return(-log(1-p) * Tski);
            else if (binomN == 1)
                return(Tski*p);             // binomial, size = Tsk
            else
                return(binomN * p * Tski);  // binomial, size = binomN
        }
    }
    return(0);
}
//==============================================================================

double adjustp (
        const double p,     // gk[i3(c, k, i, cc, kk)] etc.
        const int binomN, 
        const int detect, 
        const double Tski
)
{
    if (p < -0.1) { 
        return(0.0);
    }  
    if (p>0) {
        if (detect == 1) {   // binary proximity 
            if (fabs(Tski-1) > 1e-10)   // adjust for usage only if needed
                return(1 - pow(1-p, Tski));
            else
                return(p);
        }
        else if (detect == 2) {       // count proximity 
            if (binomN == 0)                // Poisson
                return(-log(1-p) * Tski);
            else if (binomN == 1)
                return(1 - pow(1-p, Tski));          // binomial, size = Tsk
            else
                return(1 - pow(1-p, binomN * Tski));  // binomial, size = binomN
        }
    }
    return(0);
}
//==============================================================================

// [[Rcpp::export]]
NumericVector expdetectpointcpp (
        const int           &detect,     // detector -1 single, 0 multi, 1 proximity, 2 count,... 
        const int           &N, 
        const int           &cc0,
        const int           &cc,
        const NumericVector &gk0, 
        const NumericVector &gk, 
        const NumericVector &hk0, 
        const NumericVector &hk, 
        const IntegerVector &PIA0,       // lookup which g0/sigma/b combination to use for given g, S, K [naive animal] 
        const IntegerVector &PIA1,       // lookup which g0/sigma/b combination to use for given n, S, K  [caught before] 
        const int           &nmix,       // number of classes 
        const IntegerVector &knownclass, // known membership of 'latent' classes 
        const NumericVector &pmix,       // membership probabilities
        const NumericMatrix &Tsk,        // ss x kk array of 0/1 usage codes or effort 
        const int           &btype,      // code for behavioural response  0 none etc. 
        const int           &Markov,     // learned vs transient behavioural response 0 learned 1 Markov 
        const IntegerVector &binomN      // number of trials for 'count' detector modelled with binomial 
)
{
    //  detect may take values -
    //  0  multi-catch traps
    //  1  binary proximity detectors
    //  2  count  proximity detectors
    //            binomN 0 Poisson, 1 Bernoulli, >1 binomial
    
    // 2024-08-02 output order of expdetectpointcpp switched to i,s,k
    
    int    kk = Tsk.nrow();            // number of detectors 
    int    ss = Tsk.ncol();            // number of occasions
    
    double p,ps;
    int    i,k,s;
    int    ik;
    int    wxi   = 0;
    int    c     = 0;
    int    c0    = 0;
    double Tski  = 1.0;  
    double pbefore = 0.0;
    
    std::vector<int> caughtbefore(N * kk, 0);
    std::vector<int> x(N, 0);          // mixture class of animal i 
    
    // return values
    NumericVector value (N*ss*kk);         // return value array

    //========================================================
    // 'multi-catch only' declarations 
    std::vector<double> h(N * kk);        // multi-catch only 
    std::vector<double> hsum(N);          // multi-catch only 
    
    //========================================================
    // MAIN LINE 
    
    if (btype>0) {
        Rcpp::stop("expdetectpointcpp not ready for learned responses");
    }
    
    //----------------------------------------------------------------------------
    // mixture models 
    if (nmix>1) {
        if (nmix>2)
            Rcpp::stop("nmix>2 not implemented");
        for (i=0; i<N; i++) {
            if (knownclass[i] > 1) 
                x[i] = knownclass[i] - 2;          // knownclass=2 maps to x=0 etc. 
            else
                x[i] = rdiscrete(nmix, pmix) - 1;  // assigned at random - a STOPGAP
        }
    }
    
    // ------------------------------------------------------------------------- 
    // MAIN LOOP 
    
    // -------------------------------------------------------------------------- 
    // multi-catch trap; max one site per animal per occasion 
    if (detect == 0) {
        for (i=0; i<N; i++) {
            ps = 0.0;
            pbefore = 0.0;
            for (s=0; s<ss; s++) {
                hsum[i] = 0;
                pbefore += ps;  // general learned response for starters
                for (k=0; k<kk; k++) {
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        wxi =  i4(i, s, k, x[i], N, ss, kk);
                        c0 = PIA0[wxi] - 1;
                        c  = PIA1[wxi] - 1;
                        if (c >= 0) {    // ignore unused detectors 
                            h[k * N + i] = 
                                pbefore       * Tski * hk[i3(c, k, i, cc, kk)] +
                                (1 - pbefore) * Tski * hk0[i3(c0, k, i, cc0, kk)];
                            hsum[i] += h[k * N + i];
                        }
                    }
                }
                // probability newly detected on this occasion:
                ps = (1-pbefore) * (1 - exp(-hsum[i]));   
                for (k=0; k<kk; k++) {
                    value[i3(i, s, k, N, ss)] = ps * h[k * N + i]/hsum[i];
                }
            }
        }
    }
    // -------------------------------------------------------------------------------- 
    // the 'proximity' group of detectors 1:2 - proximity, count
    
    // binomN = 0 — Poisson
    // binomN = 1 — binomial with size determined by usage = Tsk
    // binomN > 1 — binomial with size binomN
    
    else if ((detect == 1) || (detect == 2)) {
        for (i=0; i<N; i++) {
            ps = 0.0;
            pbefore = 0.0;
            for (s=0; s<ss; s++) {
                // general learned response for starters;
                // later depends on btype && Markov
                pbefore += ps; 
                for (k=0; k<kk; k++) {
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        wxi =  i4(i, s, k, x[i], N, ss, kk);
                        c0 = PIA0[wxi] - 1;
                        c  = PIA1[wxi] - 1;
                        if (c >= 0) {    // ignore unused detectors 
                            if (btype>0) {
                                // sum, weighted by probability of previous capture
                                value[i3(i, s, k, N, ss)] = 
                                    pbefore       * Ey(gk[i3(c, k, i, cc, kk)],    binomN[s], detect, Tski) +
                                    (1 - pbefore) * Ey(gk0[i3(c0, k, i, cc0, kk)], binomN[s], detect, Tski);
                                // probability newly detected on this occasion:
                                ps = (1-pbefore) * adjustp(gk0[i3(c0, k, i, cc0, kk)], binomN[s], detect, Tski);
                            }
                            else {
                                value[i3(i, s, k, N, ss)] = Ey(gk[i3(c, k, i, cc, kk)], binomN[s], detect, Tski);
                            }
                        }
                    }
                }  // loop over k
            }      // loop over s
        }          // loop over i
    }
    else {
        Rcpp::stop ("unrecognised or unsupported detector in expdetectpointcpp");
    }
    
    return (value);

}
//==============================================================================
