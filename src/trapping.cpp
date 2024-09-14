#include "poly.h"
using namespace Rcpp;

/*
 trappingXXXX routines perform simulated sampling of 2D popn with various
 detector types
 
 trappingsingle
 trappingmulti
 trappingcapped
 trappingproximity
 trappingcount
 trappingpolygon
 trappingpolygonX
 trappingtransect
 trappingtransectX
 trappingsignal
 trappingtelemetry
 
 2019-07-29 C++
 2019-08-14 slimmed down - pass dist2 etc. 

*/

// [[Rcpp::export]]
List trappingsingle (
    const NumericVector &g0,        // Parameter : detection magnitude  
    const NumericVector &sigma,     // Parameter : detection scale 
    const NumericVector &z,         // Parameter : detection shape (hazard) 
    const NumericMatrix &dist2,     // distances squared (optional: -1 if unused) 
    const NumericMatrix &Tsk,       // ss x kk array of 0/1 usage codes or effort 
    const int    fn,                // code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 uniform 
    const double w2,                // truncation radius 
    const IntegerVector &binomN,    // not used 
    const bool bk
) 
{

  int    ss = Tsk.ncol();    // number of occasions 
  int    kk = dist2.nrow();  // number of traps 
  int    N  = dist2.ncol();  // number of animals 
  
  // return values
  std::vector<int> lastcapt(N);  // occasion of last detection; 0 if none 
  std::vector<int> caught(N);    // caught in session 
  std::vector<int> value(ss*N);  // return value matrix of trap locations n x s 

  int    i,j,k,s,gi,si;
  int    nc         = 0;
  int    tr_an_indx = 0;
  double p;
  int nanimals;        // temporary 
  int ntraps;          // temporary 
  std::vector<int> occupied(kk);    // today 
  std::vector<int> intrap(N);       // today   
  
  //struct  trap_animal *tran;
  std::vector<trap_animal> tran(N*kk);
  IntegerMatrix at(N,kk);
  double event_time;
  int anum = 0;
  int tnum = 0;
  int nextcombo;
  int finished;
  int OK;
  double Tski;
  std::vector<double> gsb(3);
  std::vector<double> miscparm(4);

  // MAIN LINE 
  for (i=0; i<N; i++) {
    caught[i] = 0;   // has animal i been caught in session? 
  }
  for (s=0; s<ss; s++) {
    // initialise day 
    tr_an_indx = 0;
    nanimals = N;
    ntraps   = kk;
    for (i=0; i<N; i++) intrap[i] = 0;
    for (k=0; k<kk; k++) occupied[k] = 0;
    nextcombo = 0;
    
    // make tran 
    for (i=0; i<N; i++) {                     // animals 
      for (k=0; k<kk; k++) {                  // traps 
        Tski = Tsk(k,s);
        if (fabs(Tski) > 1e-10) {         
            if (bk) {
                gi = at(i,k) * ss * kk + k * ss + s;
            }
            else {
                gi = (lastcapt[i]>0) * ss * kk + k * ss + s;
            }
          si = (lastcapt[i]>0) * ss + s;
          gsb[0] = g0[gi];
          gsb[1] = sigma[si];
          gsb[2] = z[s];
          p = pfnS(fn, dist2(k, i), gsb, miscparm, w2); 
          if (fabs(Tski-1) > 1e-10)       
            p = 1 - exp(Tski * log(1 - p));
          
          event_time = randomtime(p);
          if (event_time <= 1) {
            tran[tr_an_indx].time   = event_time;
            tran[tr_an_indx].animal = i;    // 0..*N-1 
            tran[tr_an_indx].trap   = k;    // 0..kk-1 
            tr_an_indx++;
          }
        }
      }
    }
    
    if (tr_an_indx>0) probsort (tr_an_indx, tran);
    
    // make captures 
    while ((nextcombo < tr_an_indx) && (nanimals>0) && (ntraps>0)) {
      finished = 0;
      OK       = 0;
      while ((1-finished)*(1-OK) > 0) {    // until finished or OK 
        if (nextcombo >= (tr_an_indx)) finished = 1;  // no more to process 
        else {
          anum = tran[nextcombo].animal;
          tnum = tran[nextcombo].trap;
          OK = (1-occupied[tnum]) * (1-intrap[anum]); // not occupied and not intrap 
          nextcombo++;
        }
      }
      if (finished==0) {                   // Record this capture 
        occupied[tnum] = 1;
        intrap[anum]   = tnum+1;       // trap = k+1 
        nanimals--;
        ntraps--;
        at(anum,tnum) = 1;
      }
    }
    for (i=0; i<N; i++)
      if (intrap[i]>0) {
        if (caught[i]==0) {                  // first capture of this animal 
          nc++;
          caught[i] = nc;                   // nc-th animal to be captured 
          for (j=0; j<ss; j++)
            value[ss * (nc-1) + j] = 0;
        }
        lastcapt[i] = s+1;
        value[ss * (caught[i]-1) + s] = intrap[i];  // trap = k+1 
      }
  }
  
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("lastcapt") = wrap(lastcapt),
                       Named("caught") = wrap(caught),
                       Named("value") = wrap(value)));
}
//==============================================================================

// [[Rcpp::export]]
List trappingmulti (
    const NumericVector &g0,        // Parameter : detection magnitude  
    const NumericVector &sigma,     // Parameter : detection scale 
    const NumericVector &z,         // Parameter : detection shape (hazard) 
    const NumericMatrix &dist2,     // distances squared (optional: -1 if unused) 
    const NumericMatrix &Tsk,       // ss x kk array of 0/1 usage codes or effort 
    const int    fn,        // code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 uniform 
    const double w2,        // truncation radius 
    const IntegerVector &binomN,  // not used 
    const bool bk
)
  
{
  int    ss = Tsk.ncol();    // number of occasions 
  int    kk = dist2.nrow();  // number of traps 
  int    N  = dist2.ncol();   // number of animals 

  // return values
  std::vector<int> lastcapt(N);  // occasion of last detection; 0 if none 
  std::vector<int> caught(N);    // caught in session 
  std::vector<int> value(ss*N);  // return value matrix of trap locations n x s 

  NumericMatrix h(N, kk);
  NumericVector hsum(N);
  NumericVector cump(kk+1);
  IntegerMatrix at(N,kk);
  std::vector<double> gsb(3);
  std::vector<double> miscparm(4);
  
  double runif;
  int    i,j,k,s,gi,si;
  int    nc;
  double p;
  double Tski;
  
  cump[0] = 0;
  nc = 0;
  for (s=0; s<ss; s++) {
    for (i=0; i<N; i++) {
      hsum(i) = 0;
      for (k=0; k<kk; k++)
      {
        if (bk) {
          gi = at(i,k) * ss * kk + k * ss + s;
        }
        else {
          gi = (lastcapt[i]>0) * ss * kk + k * ss + s;
        }
        si = (lastcapt[i]>0) * ss + s;
        gsb[0] = g0[gi];
        gsb[1] = sigma[si];
        gsb[2] = z[s];
        p = pfnS(fn,  dist2(k, i), gsb, miscparm, w2); 
        Tski = Tsk(k,s);
        //Rprintf("fabs(Tski) %8.6f p %8.6f\n", fabs(Tski), p);
        if (fabs(Tski) > 1e-10) {         
          h(i,k) = -Tski * log(1 - p);
          hsum(i) += h(i,k);
        }
      }
      
      for (k=0; k<kk; k++) {
        cump(k+1) = cump(k) + h(i,k)/hsum(i);
        //Rprintf("cump %8.6f\n", cump(k+1));        
      }
      
      if (unif_rand() < (1-exp(-hsum(i))))
      {
        if (caught[i]==0)           // first capture of this animal 
        {
          nc++;
          caught[i] = nc;
          for (j=0; j<ss; j++)
            value[ss * (nc-1) + j] = 0;
        }
        runif = unif_rand();
        k = 0;
        while ((runif > cump(k)) && (k<kk)) k++;  // pick a trap 
        // Rprintf("k %6d\n", k);        
        lastcapt[i] = s+1;
        at(i,k-1) = 1;
        value[ss * (caught[i]-1) + s] = k;
        
      }
    }
  }
  
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("lastcapt") = wrap(lastcapt),
                       Named("caught") = wrap(caught),
                       Named("value") = wrap(value)));
}
//==============================================================================

// [[Rcpp::export]]
List trappingcapped (    
        const NumericVector &g0,        // Parameter : detection magnitude  
        const NumericVector &sigma,     // Parameter : detection scale 
        const NumericVector &z,         // Parameter : detection shape (hazard) 
        const NumericMatrix &dist2,     // distances squared (optional: -1 if unused) 
        const NumericMatrix &Tsk,       // ss x kk array of 0/1 usage codes or effort 
        const int    fn,                // code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 uniform 
        const double w2,                // truncation radius 
        const IntegerVector &binomN,    // not used 
        const bool bk
)
    
{
    int    ss = Tsk.ncol();    // number of occasions 
    int    kk = dist2.nrow();  // number of traps 
    int    N  = dist2.ncol();   // number of animals 

  // return values
  std::vector<int> caught(N);    // caught in session 
  std::vector<int> value(kk*ss*N);  // return value array of trap locations n x s x k
  
  std::vector<double> h(N);
  double hsum;
  std::vector<double> cumh(N+1);
  double runif;
  int    i,k,s;
  int    gi, si;
  int    nc = 0;
  double p;
  std::vector<double> gsb(3);
  std::vector<double> miscparm(4);
  double Tski;
  
  cumh[0] = 0;
  for (s=0; s<ss; s++) {
    for (k=0; k<kk; k++) {
      hsum = 0;
      Tski = Tsk(k,s);
      if (fabs(Tski) > 1e-10) {    // don't bother if unused 
        for (i=0; i<N; i++) {
          gi = k * ss + s;
          si = s;
          gsb[0] = g0[gi];
          gsb[1] = sigma[si];
          gsb[2] = z[s];
          p = pfnS(fn, dist2(k, i), gsb, miscparm, w2); 
          h[i] = -Tski * log(1 - p);
          hsum += h[i];
        }    
        if (hsum > 0) {       
          // cumulative probability across animals, conditional on one detected 
          for (i=0; i<N; i++) {
            cumh[i+1] = cumh[i] + h[i];
          }
          
          if (unif_rand() < (1-exp(-hsum))) { // success at this detector 
            // pick an animal 
            runif = unif_rand();
            i = 0;
            while ((runif > (cumh[i+1]/hsum)) && (i<N)) i++; 
            
            // first capture of this animal 
            if (caught[i] == 0) {
              nc++;
              caught[i] = nc;
            }
            // Rprintf("i %4d k %4d caught[i]-1, %4d index %10d \n", i, k, caught[i]-1, i3(i, s, k, N, ss)); 
            value[i3(i, s, k, N, ss)] = 1;      
          }
        }
      }
    }
  }
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("caught") = wrap(caught),
                       Named("value") = wrap(value)));
  
}
//==============================================================================

// [[Rcpp::export]]
List trappingproximity (
    const NumericVector &g0,        // Parameter : detection magnitude  
    const NumericVector &sigma,     // Parameter : detection scale 
    const NumericVector &z,         // Parameter : detection shape (hazard) 
    const NumericMatrix &dist2,     // distances squared (optional: -1 if unused) 
    const NumericMatrix &Tsk,       // ss x kk array of 0/1 usage codes or effort 
    const int    fn,            // code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 uniform 
    const double w2,            // truncation radius 
    const IntegerVector &binomN,  // 0 poisson, 1 Bernoulli, or number of trials for 'count' detector modelled with binomial 
    const bool bk
)
{
  double theta;
  int    i,k,s;
  int    gi, si;
  int    nc = 0;
  int    count;
  std::vector<double> gsb(3);
  std::vector<double> miscparm(4);
  double Tski;
  
  int    ss = Tsk.ncol();    // number of occasions 
  int    kk = dist2.nrow();  // number of traps 
  int    N  = dist2.ncol();   // number of animals 

  // return values
  std::vector<int> caught(N);       // caught in session 
  std::vector<int> value(kk*ss*N);  // return value array of trap locations n x s x k
  
  for (s=0; s<ss; s++) {
    for (i=0; i<N; i++) {
      for (k=0; k<kk; k++) {
        Tski =Tsk(k,s);
        if (fabs(Tski) > 1e-10) {       
          gi = k * ss + s;
          si = s;
          gsb[0] = g0[gi];
          gsb[1] = sigma[si];
          gsb[2] = z[s];
          theta = pfnS(fn, dist2(k, i), gsb, miscparm, w2); 

          if (theta>0) {
            count = rcount (1, theta, Tski);
            if (count>0)
            {
              // first capture of this animal 
              if (caught[i]==0)                  
              {
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
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("caught") = wrap(caught),
                       Named("value") = wrap(value)));
}
//===========================================================================

// [[Rcpp::export]]
List trappingcount (
    const NumericVector &g0,        // Parameter : detection magnitude  
    const NumericVector &sigma,     // Parameter : detection scale 
    const NumericVector &z,         // Parameter : detection shape (hazard) 
    const NumericMatrix &dist2,     // distances squared (optional: -1 if unused) 
    const NumericMatrix &Tsk,       // ss x kk array of 0/1 usage codes or effort 
    const int    fn,        // code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 uniform 
    const double w2,        // truncation radius 
    const IntegerVector &binomN,  // 0 poisson, 1 Bernoulli, or number of trials for 'count' detector modelled with binomial 
    const bool bk
)
  
{
  int    ss = Tsk.ncol();    // number of occasions 
  int    kk = dist2.nrow();  // number of traps 
  int    N  = dist2.ncol();  // number of animals 

  // return values
  std::vector<int> caught(N);    // caught in session 
  std::vector<int> value(kk*ss*N);  // return value array of trap locations n x s x k
  
  double theta;
  int    i,k,s;
  int    gi, si;
  int    nc = 0;
  int    count;
  std::vector<double> gsb(3);
  std::vector<double> miscparm(4);
  double Tski;
  
  for (s=0; s<ss; s++) {
    for (i=0; i<N; i++) {
      for (k=0; k<kk; k++) {
        Tski = Tsk(k,s);
        if (fabs(Tski) > 1e-10) {          // nonzero 2012 12 18 
          gi = k * ss + s;
          si = s;
          gsb[0] = g0[gi];
          gsb[1] = sigma[si];
          gsb[2] = z[s];
          theta = pfnS(fn, dist2(k, i), gsb, miscparm, w2); 
          if (theta>0) {
            if (binomN[s] == 1) {
              count = rcount (round(Tski), theta, 1);
            }
            else {
              if (binomN[s] == 0)   
                theta = hazard(theta);
              count = rcount (binomN[s], theta, Tski);
            }
            if (count>0)
            {
              // first capture of this animal 
              if (caught[i]==0)                  
              {
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
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("caught") = wrap(caught),
                       Named("value") = wrap(value)));
}
//===========================================================================

// [[Rcpp::export]]
List trappingpolygon (
    const NumericVector &lambda,  // Parameter : expected detection events per hectare 
    const NumericVector &sigma,   // Parameter : detection scale 
    const NumericVector &z,       // Parameter : detection shape (hazard) 
    const int    npoly,          // number of different polygons 
    const IntegerVector &kk,      // number of vertices + 1 (assumed closed) each polygon 
    const NumericMatrix &animals, // x,y points of animal range centres (first x, then y)  
    const NumericMatrix &traps,   // x,y polygon vertices (first x, then y)  
    const NumericMatrix &Tsk,     // ss x kk array of 0/1 usage codes or effort 
    const int    fn,             // code 0 = halfnormal, 1 = hazard, 2 = exponential 
    const double w2,             // truncation radius 
    const IntegerVector &binomN,  // 0 poisson, 1 Bernoulli, or number of binomial trials 
    const int    maxperpoly      
)
  
{
  int    ss = Tsk.ncol();    // number of occasions 
  int    N = animals.nrow();   // number of animals 

  // return values
  int    maxdet;
  maxdet = N * ss * npoly * maxperpoly;
  std::vector<int> caught(N);          // caught in session 
  std::vector<double> detectedXY(maxdet*2); // x,y locations of detections  
  std::vector<int> value(maxdet*ss);        // return value matrix of trap locations n x s 

  int    i,j,k,l,s,t;
  int    nc = 0;
  int    nd = 0;
  int    count;
  NumericVector gsb(3);
  double w;
  double ws;
  NumericVector xy;
  int cumk[maxnpoly+1];   // limit maxnpoly polygons 
  int n1,n2;
  double *workXY;
  int    *sortorder;
  double *sortkey;
  
  int J;
  double dx,dy,d;
  double Tski = 1.0;
  
  if (npoly>maxnpoly) {
    return (List::create(Named("resultcode") = 2, 
                         Named("n") = nc, 
                         Named("caught") = wrap(caught),
                         Named("detectedXY") = wrap(detectedXY),
                         Named("value") = wrap(value)));
  }
  
  cumk[0] = 0;
  for (k =0; k<npoly; k++) cumk[k+1] = cumk[k] + kk[k];
  workXY = (double*) R_alloc(maxdet*2, sizeof(double));
  sortorder = (int*) R_alloc(maxdet, sizeof(int));
  sortkey = (double*) R_alloc(maxdet, sizeof(double));
  
  // find maximum distance between animal and detector vertex 
  w = 0;
  J = cumk[npoly];
  for (i = 0; i< N; i++) {
    for (j = 0; j < J; j++) {
      dx = animals(i,0) - traps(j,0);
      dy = animals(i,1) - traps(j,1);
      d = std::sqrt(dx*dx + dy*dy);
      if (d > w) w = d;
    }
  } 
  for (i=0; i<N; i++) caught[i] = 0;
  for (s=0; s<ss; s++) {
    if (w > (10 * sigma[s])) 
      ws = 10 * sigma[s];
    else 
      ws = w;
    
    gsb(0) = lambda[s];
    gsb(1) = sigma[s];
    gsb(2) = z[s];
    if (lambda[s]>0) {
      for (i=0; i<N; i++) {
        count = rcount (binomN[s], lambda[s], Tski);   // NOT YET DEFINED !!!!
        // require maximum at r=0 
        if (fn == 6) Rcpp::stop ("annular normal not allowed in trappingpolygon");
        gsb(0) = 1;
        for (j=0; j<count; j++) {
          xy = gxy (fn, gsb, ws);            // simulate location 
          xy[0] = xy[0] + animals[i];
          xy[1] = xy[1] + animals[N + i];
          for (k=0; k<npoly; k++) {             // each polygon 
            Tski = Tsk[s * npoly + k];
            if (fabs(Tski) > 1e-10) {          // 2012 12 18 
              n1 = cumk[k];
              n2 = cumk[k+1]-1;
              // assume polygon closed 
              if (insidecpp(xy, n1, n2, traps)) {
                if (caught[i]==0) {            // first capture of this animal 
                  nc++;
                  caught[i] = nc;
                  for (t=0; t<ss; t++)
                    for (l=0; l<npoly; l++)
                      value[ss * ((nc-1) * npoly + l) + t] = 0;
                }
                nd++;
                if (nd >= maxdet) {
                  return (List::create(Named("resultcode") = 2, 
                                       Named("n") = nc, 
                                       Named("caught") = wrap(caught),
                                       Named("detectedXY") = wrap(detectedXY),
                                       Named("value") = wrap(value)));
                }
                value[ss * ((caught[i]-1) * npoly + k) + s]++;
                workXY[(nd-1)*2] = xy[0];
                workXY[(nd-1)*2+1] = xy[1];
                sortkey[nd-1] = (double) (k * N * ss + s * N + caught[i]);
              }
            }
          }
        }
      }
    }
  }
  for (i=0; i<nd; i++)
    sortorder[i] = i;
  if (nd>0) rsort_with_index (sortkey, sortorder, nd);
  for (i=0; i<nd; i++) {
    detectedXY[i]    = workXY[sortorder[i]*2];
    detectedXY[i+nd] = workXY[sortorder[i]*2+1];
  }
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("caught") = wrap(caught),
                       Named("detectedXY") = wrap(detectedXY),
                       Named("value") = wrap(value)));
}
//==============================================================================

// [[Rcpp::export]]
List trappingtransect (
    const NumericVector &lambda,   // Parameter : expected detection events per hectare 
    const NumericVector &sigma,    // Parameter : detection scale 
    const NumericVector &z,        // Parameter : detection shape (hazard) 
    const int           ntransect, // number of different transects
    const IntegerVector &kk,       // number of vertices each transect
    const NumericMatrix &animals,  // x,y points of animal range centres (first x, then y)  
    const NumericMatrix &traps,    // x,y transect vertices (first x, then y)  
    const NumericMatrix &Tsk,      // ss x kk array of 0/1 usage codes or effort 
    const int           fn,        // code 14 = HHN, 15 = HHR, 16 = HEX
    const double        w2,        // truncation radius 
    const IntegerVector &binomN,   // 0 poisson, 1 Bernoulli, or number of binomial trials 
    const int           maxperpoly  
)
{
  int    ss = Tsk.ncol();                   // number of occasions 
  int    N = animals.nrow();                // number of animals 
  int    maxdet;
  maxdet = N * ss * ntransect * maxperpoly;
  std::vector<int> caught(N);               // caught in session 
  std::vector<double> detectedXY(maxdet*2); // x,y locations of detections  
  std::vector<int> value(maxdet*ss);        // return value matrix n x s 

  int    i,j,k,l,s,t;
  int    nc = 0;
  int    nd = 0;
  int    sumk;
  int    n1,n2;
  int    count;
  bool   gotcha;
  
  std::vector<int> cumk(ntransect+1);
  double *cumd;
  rpoint *line;
  rpoint xy;
  rpoint animal;
  NumericVector gsb(3); 
  NumericMatrix gsbval(1,3); 
  double lx;
  double maxg = 0;
  double lambdak;
  double grx;
  double *workXY;
  int    *sortorder;
  double *sortkey;
  double H;
  double Tski;
  
  cumk[0] = 0;
  for (k =0; k<ntransect; k++) {
    cumk[k+1] = cumk[k] + kk[k];
  }
  sumk = cumk[ntransect];
  line      = (rpoint *) R_alloc(sumk, sizeof(rpoint));
  cumd      = (double *) R_alloc(sumk+1, sizeof(double));
  workXY    = (double*)  R_alloc(maxdet*2, sizeof(double));
  sortorder = (int*)     R_alloc(maxdet, sizeof(int));
  sortkey   = (double*)  R_alloc(maxdet, sizeof(double));

  // coordinates of vertices 
  for (i=0; i<sumk; i++) {
    line[i].x = traps[i];
    line[i].y = traps[i+sumk];
  }
  
  // cumulative distance along lines, transects end on end 
  for (k=0; k<ntransect; k++) {
    cumd[cumk[k]] = 0;
    for (i=cumk[k]; i<(cumk[k+1]-1); i++) {
      cumd[i+1] = cumd[i] + distance1(line[i], line[i+1]);
    }
  }

  for (i=0; i<N; i++) caught[i] = 0;
  for (s=0; s<ss; s++) {                            // each occasion 
    if (lambda[s]>0) {
      for (i=0; i<N; i++) {                     // each animal 
        animal.x = animals(i,0);
        animal.y = animals(i,1);
        for (k=0; k<ntransect; k++) {         // each transect 
          // if (used[s * ntransect + k]) {  
          Tski = Tsk[s * ntransect + k];
          if (fabs(Tski) > 1e-10) {          // 2012 12 18 
            
            n1 = cumk[k];
            n2 = cumk[k+1]-1;
            gsb(0) = lambda[s];
            gsb(1) = sigma[s];
            gsb(2) = z[s];
            H = hintegral1Ncpp(fn, as<std::vector<double>>(gsb)) / gsb(0);
            // H = hintegral1(*fn, par) / par[0];
            
            // flaw in following: integral1D can exceed diameter 
            gsbval(0,0) = gsb(0);
            gsbval(0,1) = gsb(1);
            gsbval(0,2) = gsb(2);
            
            // conversions for RMatrix input to integralxDNRcpp
            const RcppParallel::RMatrix<double> gsbvalR(gsbval);
            const RcppParallel::RMatrix<double> trapsR(traps);
            const RcppParallel::RMatrix<double> animalsR(animals);
            
            lambdak = gsb(0) * integral1DNRcpp (fn, i, 0, gsbvalR, trapsR, animalsR, n1, n2) / H;
            // lambdak = par[0] * integral1D (*fn, i, 0, par, 1, traps, animals, n1, n2, sumk, *N, ex) / H;
            
            count = rcount(binomN[s], lambdak, Tski);
            maxg = 0;
            if (count>0) {    // find maximum - approximate 
              for (l=0; l<=100; l++) {
                lx = cumd[n2] * l/100;
                xy = getxy (lx, cumd, line, sumk, n1);
                grx = gr (fn, gsb, xy, animal);
                // Rprintf("xy %12.6f, %12.6f grx %12.6f\n", xy.x, xy.y, grx); 
                maxg = R::fmax2(maxg, grx);
              }
              for (l=n1; l<=n2; l++) {
                xy = line[l];
                grx = gr (fn, gsb, xy, animal);
                maxg = R::fmax2(maxg, grx);
              }
              maxg= 1.2 * maxg;               // safety margin 
              if (maxg<=0) maxg=0.0001;         // not found 
            }
            for (j=0; j<count; j++) {
              gotcha = false;
              l = 0;
              while (!gotcha) {
                // n2-n1 2011-02-08 here and elsewhere 
                // location along transect 
                lx = unif_rand() * (cumd[n2]-cumd[n1]);   
                xy = getxy (lx, cumd, line, sumk, n1);
                grx = gr (fn, gsb, xy, animal);
                
                // rejection sampling 
                if (fn == 4) {
                  if (grx > 1e-10)
                    gotcha = true;
                }
                else {
                  if (unif_rand() < (grx/maxg))    
                    gotcha = true;
                }
                l++;
                if (l % 10000 == 0)
                  Rcpp::checkUserInterrupt();
                // give up and accept anything!!!! 
                // why why why? 
                if (l>1e6) {
                  Rprintf ("trials exceeded 1e6 in trappingtransect\n"); 
                  gotcha = true;        
                }
              }
              // first capture of this animal 
              if (caught[i]==0) {               
                nc++;
                caught[i] = nc;
                for (t=0; t<ss; t++)
                  for (l=0; l<ntransect; l++)
                    value[ss * ((nc-1) * ntransect + l) 
                    + t] = 0;
              }
              nd++;
              if (nd >= maxdet) {
                return (List::create(Named("resultcode") = 2, 
                                     Named("n") = nc, 
                                     Named("caught") = wrap(caught),
                                     Named("detectedXY") = wrap(detectedXY),
                                     Named("value") = wrap(value)));
              }
              value[ss * ((caught[i]-1) * ntransect + k) + s]++;
              workXY[(nd-1)*2] = xy.x;
              workXY[(nd-1)*2+1] = xy.y;
              sortkey[nd-1] = (double) (k * N * ss + s * N + caught[i]);
              // 2021-05-17 order by occasion, animal, detector
              // sortkey[nd-1] = (double) (ntransect * N * s) + ntransect*(caught[i]-1) + k;
              
            }
          }
        }
      }
    }
  }
  for (i=0; i<nd; i++)
    sortorder[i] = i;
  if (nd>0) rsort_with_index (sortkey, sortorder, nd);
  for (i=0; i<nd; i++) {
    detectedXY[i]    = workXY[sortorder[i]*2];
    detectedXY[i+nd] = workXY[sortorder[i]*2+1];
  }
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("caught") = wrap(caught),
                       Named("detectedXY") = wrap(detectedXY),
                       Named("value") = wrap(value)));
}
//==============================================================================

// [[Rcpp::export]]
List trappingpolygonX (
    const NumericVector &g0,      // Parameter : expected detection events per hectare 
    const NumericVector &sigma,   // Parameter : detection scale 
    const NumericVector &z,       // Parameter : detection shape (hazard) 
    const int    npoly,           // number of different polygons
    const IntegerVector &kk,      // number of vertices + 1 (assumed closed) each polygon 
    const NumericMatrix &animals, // x,y points of animal range centres (first x, then y)  
    const NumericMatrix &traps,   // x,y polygon vertices (first x, then y)  
    const NumericMatrix &Tsk,     // ss x kk array of 0/1 usage codes or effort 
    const int    fn,              // code 14 = HHN, 15 = HHR, 16 = HEX
    const double w2,              // truncation radius 
    const IntegerVector &binomN   // 0 poisson, 1 Bernoulli, or number of binomial trials 
)
  
{
  int    ss = Tsk.ncol();    // number of occasions 
  int    N = animals.nrow(); // number of animals 
  std::vector<int> caught(N);          // caught in session 
  std::vector<double> detectedXY(N * ss *2); // x,y locations of detections  
  std::vector<int> value(N * ss);        // return value matrix of trap locations n x s 

  int    i,j,k,s,t;
  int    nc = 0;
  int    nd = 0;
  NumericVector gsb(3);
  double w;
  double ws;
  NumericVector xy;
  int cumk[maxnpoly+1];   // limit maxnpoly polygons 
  int n1,n2;
  double *workXY;
  int    *sortorder;
  double *sortkey;
  
  int J;
  double dx,dy,d;
  int maybecaught;
  double Tski;
  
  cumk[0] = 0;
  for (k =0; k<npoly; k++) cumk[k+1] = cumk[k] + kk[k];
  workXY = (double*) R_alloc(N*ss*2, sizeof(double));
  sortorder = (int*) R_alloc(N*ss, sizeof(int));
  sortkey = (double*) R_alloc(N*ss, sizeof(double));
  
  // find maximum distance between animal and detector vertex 
  w = 0;
  J = cumk[npoly];
  for (i = 0; i< N; i++) {
    for (j = 0; j < J; j++) {
      dx = animals(i,0) - traps(j,0);
      dy = animals(i,1) - traps(j,1);
      d = std::sqrt(dx*dx + dy*dy);
      if (d > w) w = d;
    }
  } 
  for (i=0; i<N; i++) caught[i] = 0;
  
  for (s=0; s<ss; s++) {
    if (w > (10 * sigma[s])) 
      ws = 10 * sigma[s];
    else 
      ws = w;
    
    gsb(0) = g0[s];
    gsb(1) = sigma[s];
    gsb(2) = z[s];
    
    if (g0[s]>0) {
      for (i=0; i<N; i++) {
        maybecaught = unif_rand() < g0[s];
        // require maximum at r=0 
        if (fn == 6) Rcpp::stop ("annular normal not allowed in trappingpolygonX");
        
        gsb(0) = 1;
        
        // simulate location 
        if (maybecaught) {
          xy = gxy (fn, gsb, ws);
          xy[0] = xy[0] + animals(i,0);
          xy[1] = xy[1] + animals(i,1);
          // which polygon, if any? 
          for (k=0; k < npoly; k++) {             // each polygon 
            // if (used[s * npoly + k]) { 
            Tski = Tsk[s * npoly + k];
            if (fabs(Tski) > 1e-10) {          // 2012 12 18 
              
              n1 = cumk[k];
              n2 = cumk[k+1]-1;
              // assume polygon closed 
              if (insidecpp(xy, n1, n2, traps)) {
                if (caught[i]==0) {            // first capture of this animal 
                  nc++;
                  caught[i] = nc;
                  for (t=0; t<ss; t++)
                    value[ss * (nc-1) + t] = 0;
                }
                nd++;
                value[ss * (caught[i]-1) + s] = k+1;
                workXY[(nd-1)*2] = xy[0];
                workXY[(nd-1)*2+1] = xy[1];
                sortkey[nd-1] = (double) (s * N + caught[i]);
                // 2021-05-17 order by occasion, animal, detector
                // sortkey[nd-1] = (double) (npoly * N * s) + npoly*(caught[i]-1) + k;
                
                break;   // no need to look at more poly 
              }
            }
          }
        }
      }
    }
  }
  for (i=0; i<nd; i++)
    sortorder[i] = i;
  if (nd>0) rsort_with_index (sortkey, sortorder, nd);
  for (i=0; i<nd; i++) {
    detectedXY[i]    = workXY[sortorder[i]*2];
    detectedXY[i+nd] = workXY[sortorder[i]*2+1];
  }
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("caught") = wrap(caught),
                       Named("detectedXY") = wrap(detectedXY),
                       Named("value") = wrap(value)));
}
//==============================================================================

// [[Rcpp::export]]
List trappingtransectX (
    const NumericVector &lambda,  // Parameter : expected detection events per hectare 
    const NumericVector &sigma,   // Parameter : detection scale 
    const NumericVector &z,       // Parameter : detection shape (hazard) 
    const int    ntransect,       // number of different transects
    const IntegerVector &kk,      // number of vertices each transect
    const NumericMatrix &animals, // x,y points of animal range centres (first x, then y)  
    const NumericMatrix &traps,   // x,y polygon vertices (first x, then y)  
    const NumericMatrix &Tsk,     // ss x kk array of 0/1 usage codes or effort 
    const int    fn,              // code 14 = HHN, 15 = HHR, 16 = HEX
    const double w2               // truncation radius 
)
{
  int    ss = Tsk.ncol();    // number of occasions 
  int    N = animals.nrow(); // number of animals 
  int    maxdet;
  maxdet = N * ss;
  std::vector<int> caught(N);          // caught in session 
  std::vector<double> detectedXY(maxdet*2); // x,y locations of detections  
  std::vector<int> value(maxdet);        // return value matrix of trap locations n x s 

  int    i,k,l,s,t;
  int    nc = 0;
  int    nd = 0;
  int    sumk;
  int    n1,n2;
  int    count = 0;
  bool   gotcha;
  
  int    *cumk;
  double *cumd;
  rpoint *line;
  rpoint xy;
  rpoint animal;
  NumericVector gsb(3);
  NumericMatrix gsbval(1,3);
  double lx;
  double maxg = 0;
  double lambdak;
  double grx;
  double *workXY;
  int    *sortorder;
  double *sortkey;
  // double *ex;
  double H;
  double sumhaz;
  double pks;
  double Tski;
  
  cumk = (int *) R_alloc(ntransect+1, sizeof(int));
  cumk[0] = 0;
  for (k =0; k<ntransect; k++)
    cumk[k+1] = cumk[k] + kk[k];
  sumk = cumk[ntransect];
  line = (rpoint *) R_alloc(sumk, sizeof(rpoint));
  cumd = (double *) R_alloc(sumk+1, sizeof(double));
  
  workXY = (double*) R_alloc(maxdet*2, sizeof(double));
  sortorder = (int*) R_alloc(maxdet, sizeof(int));
  sortkey = (double*) R_alloc(maxdet, sizeof(double));
  // ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
  
  // coordinates of vertices 
  for (i=0; i<sumk; i++) {
    line[i].x = traps(i,0);
    line[i].y = traps(i,1);
  }
  
  // cumulative distance along line 
  for (k=0; k<ntransect; k++) {
    cumd[cumk[k]] = 0;
    for (i=cumk[k]; i<(cumk[k+1]-1); i++) {
      cumd[i+1] = cumd[i] + distance1(line[i], line[i+1]);
    }
  }

  for (s=0; s<ss; s++) {                               // each occasion 
    if (lambda[s]>0) {
      
      gsb(0) = lambda[s];
      gsb(1) = sigma[s];
      gsb(2) = z[s];
      H = hintegral1Ncpp(fn, as<std::vector<double>>(gsb));
      // H = hintegral1(*fn, par);
      
      gsbval(0,0) = gsb(0);
      gsbval(0,1) = gsb(1);
      gsbval(0,2) = gsb(2);
      
      // conversions for RMatrix input to integralxDNRcpp
      const RcppParallel::RMatrix<double> gsbvalR(gsbval);
      const RcppParallel::RMatrix<double> trapsR(traps);
      const RcppParallel::RMatrix<double> animalsR(animals);
      
      for (i=0; i<N; i++) {                        // each animal 
        animal.x = animals(i,0);
        animal.y = animals(i,1);
        sumhaz = 0;
        for (k=0; k<ntransect; k++) {            // sum hazard 
          // if (used[s * ntransect + k]) { 
          Tski = Tsk[s * ntransect + k];
          if (fabs(Tski) > 1e-10) {          // 2012 12 18 
            
            n1 = cumk[k];
            n2 = cumk[k+1]-1;
            
            sumhaz += -log(1 - gsb(0) * integral1DNRcpp (fn, i, 0, gsbvalR, trapsR, 
                                                    animalsR, n1, n2) / H);
            // sumhaz += -log(1 - par[0] * integral1D (*fn, i, 0, par, 1, traps, 
            //                                         animals, n1, n2, sumk, *N, ex) / H);
            
          }
        }
        
        for (k=0; k<ntransect; k++) {            // each transect 
          // if (used[s * ntransect + k]) { 
          Tski = Tsk[s * ntransect + k];
          if (fabs(Tski) > 1e-10) {          // 2012 12 18 
            n1 = cumk[k];
            n2 = cumk[k+1]-1;
            lambdak = gsb(0) * integral1DNRcpp (fn, i, 0, gsbvalR, trapsR, 
                                           animalsR, n1, n2) / H;
            // lambdak = par[0] * integral1D (*fn, i, 0, par, 1, traps, 
            //                                animals, n1, n2, sumk, *N, ex) / H;
            pks = (1 - exp(-sumhaz)) * (-log(1-lambdak)) / sumhaz;
            count = unif_rand() < pks;
            maxg = 0;
            
            if (count>0) {                        // find maximum - approximate 
              for (l=0; l<=100; l++) {
                lx = (cumd[n2] - cumd[n1]) * l/100;    
                xy = getxy (lx, cumd, line, sumk, n1);
                grx = gr (fn, gsb, xy, animal);
                maxg = R::fmax2(maxg, grx);
              }
              for (l=n1; l<=n2; l++) {
                xy = line[l];
                grx = gr (fn, gsb, xy, animal);
                maxg = R::fmax2(maxg, grx);
              }
              maxg= 1.2 * maxg;                 // safety margin 
              if (maxg<=0) maxg=0.0001;         // not found 
              
              gotcha = false;
              l = 0;
              while (!gotcha) {
                lx = unif_rand() * (cumd[n2]-cumd[n1]);  // location along transect 
                xy = getxy (lx, cumd, line, sumk, n1);
                grx = gr (fn, gsb, xy, animal);
                gotcha = (unif_rand() < (grx/maxg));    // rejection sampling 
                l++;
                if (l % 10000 == 0)
                  Rcpp::checkUserInterrupt();
                if (l>1e8) gotcha = true;        // give up and accept anything!!!! 
              }
              if (caught[i]==0) {               // first capture of this animal 
                nc++;
                caught[i] = nc;
                for (t=0; t<ss; t++)
                  value[ss * (nc-1) + t] = 0;
              }
              nd++;
              value[ss * (caught[i]-1) + s] = k+1;
              workXY[(nd-1)*2] = xy.x;
              workXY[(nd-1)*2+1] = xy.y;
              sortkey[nd-1] = (double) (s * N + caught[i]);
            }
            if (count>0) break;   // no need to look further 
          }
        }
      }
    }
  }
  for (i=0; i<nd; i++) sortorder[i] = i;
  if (nd>0) rsort_with_index (sortkey, sortorder, nd);
  for (i=0; i<nd; i++) {
    detectedXY[i]    = workXY[sortorder[i]*2];
    detectedXY[i+nd] = workXY[sortorder[i]*2+1];
  }
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("caught") = wrap(caught),
                       Named("detectedXY") = wrap(detectedXY),
                       Named("value") = wrap(value)));
  
}
//==============================================================================

// [[Rcpp::export]]
List trappingsignal (
    const NumericVector &beta0,   // Parameter : intercept 
    const NumericVector &beta1,   // Parameter : slope 
    const NumericVector &sdS,     // Parameter : error sd 
    const double cut,            // detection threshold on transformed scale, etc. 
    const NumericVector &muN,     // noise mean 
    const NumericVector &sdN,     // noise sd 
    const double sdM,            // movement between occasions 
    const NumericMatrix &animals, // x,y points of animal range centres (first x, then y)  
    const NumericMatrix &traps,   // x,y locations of traps (first x, then y)  
    const NumericMatrix &dist2,   // distances squared (optional: -1 if unused) 
    const NumericMatrix &Tsk,     // ss x kk array of 0/1 usage codes or effort 
    const int    fn             
)
{
  // code 10 = signal strength, 11 = signal strength with sph spread 
  // returned signal strength (fn==10) is on transformed scale 
  // limited to Bernoulli count model binomN = 1 
  
  // return value
  
  int    ss = Tsk.ncol();    // number of occasions 
  int    kk = traps.nrow();  // number of detectors
  int    N  = animals.nrow(); // number of animals 
  
  int maxdet = N * ss * kk;
  std::vector<int> caught(N);    // caught in session 
  std::vector<double> signal(maxdet);    // signal strength, one per detection 
  std::vector<double> noise(maxdet);     // noise, one per detection, if signalnoise 
  std::vector<double> value(maxdet);     // return value matrix of trap locations n x s 

  double muS;
  double signalvalue;
  double noisevalue;
  int    i,j,k,l,s;
  int    nc = 0;
  int    nd = 0;
  double *worksignal;
  double *worknoise;
  int    *sortorder;
  double *sortkey;
  NumericMatrix animalss = animals;
  double Tski;
  
  worksignal = (double*) R_alloc(maxdet, sizeof(double));
  worknoise = (double*) R_alloc(maxdet, sizeof(double));
  sortorder = (int*) R_alloc(maxdet, sizeof(int));
  sortkey = (double*) R_alloc(maxdet, sizeof(double));

  for (s=0; s<ss; s++) {
    for (i=0; i<N; i++) {
      if (sdM > fuzz) {
        animalss(i,0) = animals(i,0) + norm_rand() * sdM;
        animalss(i,1) = animals(i,1) + norm_rand() * sdM;
      }
      else {
      }
      for (k=0; k<kk; k++) {
        Tski = Tsk(k,s);
        if (fabs(Tski) > 1e-10) {          // 2012 12 18 
          if ((fn == 10) || (fn == 12))
            muS  = mufnL (k, i, beta0[s], beta1[s], dist2, 0);
          else
            muS  = mufnL (k, i, beta0[s], beta1[s], dist2, 1);
          signalvalue = norm_rand() * sdS[s] + muS;
          if ((fn == 12) || (fn == 13)) {
            noisevalue = norm_rand() * sdN[s] + muN[s];  
            // 2012-09-19 shouldn't this be SminusN? 
            if ((signalvalue-noisevalue) > cut) {
              if (caught[i]==0) {              // first capture of this animal 
                nc++;
                caught[i] = nc;
                for (j=0; j<ss; j++)
                  for (l=0; l<kk; l++)
                    value[ss * ((nc-1) * kk + l) + j] = 0;
              }
              nd++;
              if (nd > maxdet) {
                return (List::create(Named("resultcode") = 2, 
                                     Named("n") = nc, 
                                     Named("caught") = wrap(caught),
                                     Named("signal") = wrap(signal),
                                     Named("noise") = wrap(noise),
                                     Named("value") = wrap(value)));
              }
              value[ss * ((caught[i]-1) * kk + k) + s] = 1;
              worksignal[nd-1] = signalvalue;
              worknoise[nd-1] = noisevalue;
              sortkey[nd-1] = (double) (k * N * ss + s * N + caught[i]);
              // 2021-05-17 order by occasion, animal, detector
              // sortkey[nd-1] = (double) (kk * N * s) + kk*(caught[i]-1) + k;
              
            }
          }
          else {
            if (signalvalue > cut) {
              if (caught[i]==0) {              // first capture of this animal 
                nc++;
                caught[i] = nc;
                for (j=0; j<ss; j++)
                  for (l=0; l<kk; l++)
                    value[ss * ((nc-1) * kk + l) + j] = 0;
              }
              nd++;
              if (nd > maxdet) {
                return (List::create(Named("resultcode") = 2, 
                                     Named("n") = nc, 
                                     Named("caught") = wrap(caught),
                                     Named("signal") = wrap(signal),
                                     Named("noise") = wrap(noise),
                                     Named("value") = wrap(value)));
              }
              value[ss * ((caught[i]-1) * kk + k) + s] = 1;
              worksignal[nd-1] = signalvalue;
              sortkey[nd-1] = (double) (k * N * ss + s * N + caught[i]);
              // sortkey[nd-1] = (double) (kk * N * s) + kk*(caught[i]-1) + k;
            }
          }
        }
      }
    }
  }
  for (i=0; i<nd; i++)
    sortorder[i] = i;
  if (nd>0) rsort_with_index (sortkey, sortorder, nd);
  for (i=0; i<nd; i++)
    signal[i]   = worksignal[sortorder[i]];
  if ((fn == 12) || (fn == 13)) {
    for (i=0; i<nd; i++)
      noise[i]   = worknoise[sortorder[i]];
  }
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("caught") = wrap(caught),
                       Named("signal") = wrap(signal),
                       Named("noise") = wrap(noise),
                       Named("value") = wrap(value)));
}
//==============================================================================


// [[Rcpp::export]]
List trappingtelemetry (
    const NumericVector &lambda,  // Parameter : expected detection events per hectare 
    const NumericVector &sigma,   // Parameter : detection scale 
    const NumericVector &z,       // Parameter : detection shape (hazard) 
    const NumericMatrix &animals, // x,y points of animal range centres (first x, then y)  
    const int    ss,             // number of occasions
    const int    fn,             // code 0 = halfnormal, 1 = hazard, 2 = exponential 
    const double w2,             // truncation radius 
    const IntegerVector &binomN,  // 0 poisson, 1 Bernoulli, or number of binomial trials 
    const int    exactn,         // 0 or a positive integer for the exact number of fixes per animal 
    const int    maxperpoly  
)
{
  int    N = animals.nrow();     // number of animals
  // return values
  int    maxdet;
  maxdet = N * ss * maxperpoly;
  std::vector<int> caught(N);          // caught in session 
  std::vector<double> detectedXY(maxdet*2); // x,y locations of detections  
  std::vector<int> value(maxdet*ss);        // return value matrix of trap locations n x s 

  int    i,j,s,t;
  int    nc = 0;
  int    nd = 0;
  int    count;
  NumericVector gsb(3);
  double ws;
  NumericVector xy;
  double *workXY;
  int    *sortorder;
  double *sortkey;
  double Tski = 1.0;   // not defined !! 2012-12-18 
  
  workXY = (double*) R_alloc(maxdet*2, sizeof(double));
  sortorder = (int*) R_alloc(maxdet, sizeof(int));
  sortkey = (double*) R_alloc(maxdet, sizeof(double));
  
  for (s=0; s<ss; s++) {
    ws = 10 * sigma[s];
    gsb(0) = lambda[s];
    gsb(1) = sigma[s];
    gsb(2) = z[s];
    if (lambda[s]>0) {
      for (i=0; i<N; i++) {
        if (exactn>0)
          count = exactn;
        else
          count = rcount (binomN[s], lambda[s], Tski); 
        // require maximum at r=0 
        if (fn == 6) Rcpp::stop ("annular normal not allowed in trappingtelemetry");
        gsb(0) = 1;
        for (j=0; j<count; j++) {
          xy = gxy (fn, gsb, ws);            // simulate location 
          xy[0] = xy[0] + animals(i,0);
          xy[1] = xy[1] + animals(i,1);
          if (caught[i]==0) {            // first capture of this animal 
            nc++;
            caught[i] = nc;
            for (t=0; t<ss; t++)
              value[ss * (nc-1) + t] = 0;
          }
          nd++;
          if (nd >= maxdet) {
            return (List::create(Named("resultcode") = 2, 
                                 Named("n") = nc, 
                                 Named("caught") = wrap(caught),
                                 Named("detectedXY") = wrap(detectedXY),
                                 Named("value") = wrap(value)));
          }
          value[ss * (caught[i]-1) + s]++;
          workXY[(nd-1)*2] = xy[0];
          workXY[(nd-1)*2+1] = xy[1];
          sortkey[nd-1] = (double) (s * N + caught[i] - 1);
        }
      }
    }
  }
  for (i=0; i<nd; i++)
    sortorder[i] = i;
  if (nd>0) rsort_with_index (sortkey, sortorder, nd);
  for (i=0; i<nd; i++) {
    detectedXY[i]    = workXY[sortorder[i]*2];
    detectedXY[i+nd] = workXY[sortorder[i]*2+1];
  }
  return (List::create(Named("resultcode") = 0, 
                       Named("n") = nc, 
                       Named("caught") = wrap(caught),
                       Named("detectedXY") = wrap(detectedXY),
                       Named("value") = wrap(value)));
}

//==============================================================================

