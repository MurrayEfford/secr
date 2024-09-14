// next two lines must be in order (RcppNumerical precedes secr.h)
// #include <RcppNumerical.h>

#include "secr.h"

using namespace std;
using namespace Rcpp;

double minimumexp = -100;

//--------------------------------------------------------------------------

double expmin (double x)
{
    if (x < minimumexp)
        return(0);
    else
        return(exp(x));
}
//--------------------------------------------------------------------------

// index to vector element corresponding to cell i,j,k in 3D array
// stored in column-major order 

int i3 (int i, int j, int k, int ii, int jj) {
    return(ii * (jj * k + j) + i);
}
//--------------------------------------------------------------------------

// index to vector element corresponding to cell i,j,k,l in 4D array
// stored in column-major order 

int i4 (int i, int j, int k, int l, int ii, int jj, int kk) {
    return (ii *(jj*(kk*l + k) + j) + i);
}
//--------------------------------------------------------------------------

double d2cpp (
        const int k,
        const int m,
        const NumericMatrix &A1,
        const NumericMatrix &A2)
    // return squared distance between two points given by row k in A1
    // and row m in A2, where A1 and A2 have respectively A1rows and A2rows
{
    return(
        (A1(k,0) - A2(m,0)) * (A1(k,0) - A2(m,0)) +
            (A1(k,1) - A2(m,1)) * (A1(k,1) - A2(m,1))
    );
}
//--------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix edist2cpp (
        const NumericMatrix &A1,
        const NumericMatrix &A2)
    // return squared distance between points in A1
    // and A2
{
    int kk = A1.nrow();
    int mm = A2.nrow();
    NumericMatrix d2(kk, mm);
    for (int k=0; k<kk; k++) {
        for (int m=0; m<mm; m++) {
            d2(k,m) =  (A1(k,0) - A2(m,0)) * (A1(k,0) - A2(m,0)) +
                (A1(k,1) - A2(m,1)) * (A1(k,1) - A2(m,1));
        }
    }
    return(d2);
}
//--------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix xydist2cpp (
        const NumericMatrix &A1,
        const NumericMatrix &A2)
    // return squared distance between points in A1
    // and A2
{
    int kk = A1.nrow();
    int mm = A2.nrow();
    NumericMatrix d2(kk, mm);
    for (int k=0; k<kk; k++) {
        for (int m=0; m<mm; m++) {
            d2(k,m) =  std::max((A1(k,0) - A2(m,0)) * (A1(k,0) - A2(m,0)),
                (A1(k,1) - A2(m,1)) * (A1(k,1) - A2(m,1)));
        }
    }
    return(d2);
}
//--------------------------------------------------------------------------

// customised dpois 
double gpois (int count, double lambda)
{
    double x;
    if ((count < 0) || (count>0 && lambda <= 0)) {
        return(0);
    }
    else if (count == 0) {
        return (exp(-lambda));
    }
    else {
        boost::math::poisson_distribution<> pois(lambda);
        x = boost::math::pdf(pois, count);
        return (x);
    }
}
//--------------------------------------------------------------------------

// customised dbinom 
double gbinom(int count, int size, double p)
{
    double x, q;
    int i;
    if ((count < 0) || (count > 0 && p <= 0)) {
        x = 0;
    }
    else if (count == 0) {
        q = 1 - p;
        x = q;
        for (i=1; i< size; i++) x *= q;
    }
    else {
        boost::math::binomial_distribution<> bin(size, p);
        x = boost::math::pdf(bin, count);
    }
    return (x);   
}
//--------------------------------------------------------------------------

// probability of count with distribution specified by binomN 
double countp (int count, int binomN, double lambda) {
    // Poisson 
    if (binomN == 0) {
        if (count == 0) 
            return (exp(-lambda));
        else {
            // return (R::dpois(count, lambda, 0));
            boost::math::poisson_distribution<> pois(lambda);
            return (boost::math::pdf(pois, count));
        }
    }
    
    // Bernoulli 
    else if (binomN == 1) {
        if (count == 0)
            return ( 1 - lambda );
        else
            return ( lambda );
    }
    
    // negative binomial 
    else if (binomN < 0) {
        boost::math::negative_binomial_distribution<> nbin(binomN, lambda);
        return (boost::math::pdf(nbin, count));
    }
    
    // binomial 
    else {
        boost::math::binomial_distribution<> bin(binomN, lambda);
        return (boost::math::pdf(bin, count));
    }
}
//--------------------------------------------------------------------------

// use distance1 rather than distance because of clash with std::distance
double distance1 (const rpoint p1, const rpoint p2) {
    return(std::sqrt ((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y)));
}
//--------------------------------------------------------------------------

double zrcpp (double r, int detectfn, NumericVector par)
{
    if (detectfn == 14) {  // hazard halfnormal
        return (exp(-r*r / 2 / par(1) / par(1)));    
    }
    else {
        if (detectfn == 15) {  // hazard hazard rate
            return (1 - exp(- pow(r /par(1), - par(2))));
        }
        else if (detectfn == 16) {  // hazard exponential
            return (exp(-r / par(1)));
        }
        else if (detectfn == 17) {  // hazard annular normal
            return (exp(-(r-par(2))*(r-par(2)) / 
                    2 / par(1)/ par(1)));
        }
        else if (detectfn == 18) {  // hazard cumulative gamma
            // return (R::pgamma(r,par(2),par(1)/par(2),0,0)); 
            boost::math::gamma_distribution<> gam(par(2),par(1)/par(2));
            return (boost::math::cdf(complement(gam,r))); 
        }
        else if (detectfn == 19) {  // hazard variable power
            return (exp(- pow(r /par(1), par(2))));
        }
        else 
            return (R_NaN);  //Rcpp::stop("unknown or invalid detection function in gxy"));
    }
}

// random points from 2-D radial distribution specified by g function 
// assume intercept par(0) == 1
NumericVector gxy (const int fn, const NumericVector par, const double w) {
    NumericVector xy(2);
    int maxj = 1000000;
    double r;
    for (int j=0; j<maxj; j++) {
        r = w * std::sqrt(unif_rand());
        if (unif_rand() < zrcpp(r, fn, par)) break;
    }
    double theta = unif_rand() * 2 * M_PI;
    xy(0) = r * cos(theta);
    xy(1) = r * sin(theta);
    return(xy);
}
//--------------------------------------------------------------------------

// [[Rcpp::export]]
List nearestcpp (
        const NumericMatrix& xy,       // input points 
        const NumericMatrix& traps,    // input 
        bool non_zero)  // If true, return nearest with non-zero distance.
{
    std::vector<int> p(xy.nrow());        // output indices of nearest point 
    std::vector<double> d (xy.nrow());        // output distances to nearest point 
    int nxy , ntrap;
    int i,j;
    int id=-1;
    double d2;
    double d2min;
    nxy = xy.nrow();
    ntrap = traps.nrow();
    for (j=0; j<nxy; j++) {
        id = -1;
        d2min = 1e100;
        for (i=0; i<ntrap; i++)
        {
            d2 = (traps(i,0) - xy(j,0)) * (traps(i,0) - xy(j,0)) +
                (traps(i,1) - xy(j,1)) * (traps(i,1) - xy(j,1));
            if (d2 < d2min && (!non_zero || d2 > 0)) {
                d2min = d2; id = i;
            }
        }
        d[j] = std::sqrt(d2min);
        p[j] = id+1;
    }
    return (List::create(Named("distance") = wrap(d), Named("index") = wrap(p)));
}
//--------------------------------------------------------------------------

// [[Rcpp::export]]
bool insidecpp (
        const NumericVector &xy,
        const int    n1,
        const int    n2,
        const NumericMatrix &poly)
{
    // Is point xy inside poly?
    // Based on contribution on s-news list by Peter Perkins 23/7/96
    // We assume poly is closed, and in col-major order (x's then y's)
    
    double theta = 0;
    double cutoff = 1e-6;
    int k;
    int ns;
    double N;
    double d;
    ns = n2 - n1 + 1;   // number of selected points 
    std::vector<double> temp((ns+1) * 2);
    
    // get & translate to coords centered at each test point 
    for (k=0; k < ns; k++)
    {
        temp[k]      = poly(k + n1,0) - xy[0];    // x 
        temp[k + ns] = poly(k + n1,1) - xy[1];    // y 
    }
    
    for (k=0; k < (ns-1); k++)
    {
        N = temp[k] * temp[k+1 + ns] - temp[k + ns] * temp[k+1];
        d = temp[k] * temp[k+1]      + temp[k + ns] * temp[k+1 + ns];
        if (fabs(d)>0) { N = N/fabs(d);  d = d/fabs(d); }
        theta += atan2(N, d);
    }
    theta = fabs(theta);
    return (fabs(theta - 2* M_PI) < cutoff);    // M_PI is cmath.h constant 
}
//--------------------------------------------------------------------------

rpoint getxy(
        const double l, 
        double cumd[], 
        const rpoint line[], 
        const int kk, 
        const int offset) {        // double before 2022-01-18
    // return the xy coordinates of point l metres along a transect 
    // offset is the starting position for this transect 
    int k = 1;    // initialised 2022-01-18
    double pr, d, d12;
    rpoint xy;
    for (k=offset+1; k<(offset+kk); k++) {
        if (cumd[k]>l) break;
    }
    k = std::min(k, offset+kk-1);   // prevent overshoot 2022-01-18
    d = l - cumd[k-1];  // distance along leg 
    
    // debug 2022-01-18
    // Rprintf("l %10.6g offset %d k %d cumd[k-1] %10.6g cumd[k] %10.6g \n", l, offset, k, cumd[k-1], cumd[k]);
        
    d12 = cumd[k] - cumd[k-1];
    if (d12>0)
        pr = d / d12;
    else
        pr = 0;
    xy.x = line[k-1].x + (line[k].x - line[k-1].x) * pr;
    xy.y = line[k-1].y + (line[k].y - line[k-1].y) * pr;
    return(xy);
}
//--------------------------------------------------------------------------

// Calculate the length of intersection of a line segment and a circle
// Based on C code of Paul Bourke November 1992
// Line segment is defined from p1 to p2
// Circle is of radius r and centred at sc
// Two potential points of intersection given by
// p = p1 + mu1 (p2-p1)
// p = p1 + mu2 (p2-p1)
// Return 0 if line segment does not intersect circle

double SegCircle2 (
        double p1x, double p1y, 
        double p2x, double p2y, 
        double scx, double scy, 
        double r
) 
{
    double a,b,c;
    double bb4ac;
    double dpx;
    double dpy;
    
    double mu1;
    double mu2;
    int p1in;
    int p2in;
    
    double i1x;
    double i1y;
    double i2x;
    double i2y;
    
    int i1between;
    double d1,d2;
    double seg = 0;
    
    // case where both p1 and p2 inside circle 
    
    // Rprintf ("p1 %6.3f %6.3f\n", p1x, p1y);
    // Rprintf ("p2 %6.3f %6.3f\n", p2x, p2y);
    // Rprintf ("sc %6.3f %6.3f\n", scx, scy);
    // Rprintf ("r %6.3f \n", r);
    
    p1in = ((scx - p1x) * (scx - p1x) + 
        (scy - p1y) * (scy - p1y)) < (r * r);
    p2in = ((scx - p2x) * (scx - p2x) + 
        (scy - p2y) * (scy - p2y)) < (r * r);
    if (p1in && p2in) {        
        seg = std::sqrt ((p1x - p2x) * (p1x - p2x) + 
            (p1y - p2y) * (p1y - p2y));
        return (seg);
    }
    
    dpx = p2x - p1x;
    dpy = p2y - p1y;
    
    a = dpx * dpx + dpy * dpy;
    b = 2 * (dpx * (p1x - scx) + dpy * (p1y - scy));
    c = scx * scx + scy * scy;
    c += p1x * p1x + p1y * p1y;
    c -= 2 * (scx * p1x + scy * p1y);
    c -= r * r;
    bb4ac = b * b - 4 * a * c;
    
    // case of no intersection 
    if ((fabs(a) < 1e-10) || (bb4ac < 0)) {
        return (0);   
    }
    
    mu1 = (-b + std::sqrt(bb4ac)) / (2 * a);
    mu2 = (-b - std::sqrt(bb4ac)) / (2 * a);
    
    i1x = p1x + mu1 * (p2x - p1x);
    i1y = p1y + mu1 * (p2y - p1y);
    i2x = p1x + mu2 * (p2x - p1x);
    i2y = p1y + mu2 * (p2y - p1y);
    
    if (((mu1<0) && (mu2<0)) || ((mu1>1) && (mu2>1))) {
        // no intersection 
        seg = 0;
    }
    else {
        if (((mu1<0) && (mu2>1)) || ((mu1>1) && (mu2<0))) {
            // both inside 
            seg = std::sqrt ((p1x - p2x) * (p1x - p2x) + 
                (p1y - p2y) * (p1y - p2y));
        }
        else {
            if ((mu1>0) && (mu1<1) && (mu2>0) && (mu2<1)) {
                // two intersections 
                seg = std::sqrt ((i1x - i2x) * (i1x - i2x) + 
                    (i1y - i2y) * (i1y - i2y));
            }
            else {
                // one intersection 
                d1 = std::sqrt((i1x - p1x) * (i1x * p1x) + 
                    (i1y - p1y) * (i1y - p1y));
                d2 = std::sqrt((i1x - p2x) * (i1x * p2x) + 
                    (i1y - p2y) * (i1y - p2y));
                i1between = std::sqrt(a) < (d1 + d2 + 1e-10);
                if (p1in) {
                    if (i1between) {
                        i2x = p1x;
                        i2y = p1y;
                    }
                    else {
                        i1x = p1x;
                        i1y = p1y;
                    }
                }
                if (p2in) {
                    if (i1between) {
                        i2x = p2x;
                        i2y = p2y;
                    }
                    else {
                        i1x = p2x;
                        i1y = p2y;
                    }
                }
                seg = std::sqrt ((i1x - i2x) * (i1x - i2x) + 
                    (i1y - i2y) * (i1y - i2y));
            }
        }
    }
    return(seg);    
}

//----------------------------------------------------------------

double randomtime (double p)
    // return random event time for event with probability p 
{
    double minprob = 1e-5;
    double lambda;
    double random_U;
    
    if (p < minprob)
        return(huge);                        // ignore trivial p/lambda 
    else if (p >= 1.0)
        return (-unif_rand());                  // trick to spread P=1 
    else {
        lambda   = -log(1-p);                // rate parameter 
        random_U = unif_rand();
        if (random_U <= 0)                   // trap for zero 
            return(huge);
        else
            return (-log(random_U)/lambda);   // random exponential e.g. Ripley 1987 Algorithm 3.2 p 55 
    }
}
//----------------------------------------------------------------

double randomtimel (double lambda)
    // return random event time for event with hazard lambda
{
    double random_U;
    if (lambda <= 0) {
        return(huge);
    }
    else {
        random_U = unif_rand();
        if (random_U <= 0) {                  // trap for zero 
            return(huge);
        }
        else {
            return (-log(random_U)/lambda);   // random exponential e.g. Ripley 1987 Algorithm 3.2 p 55 
        }
    }
}
//----------------------------------------------------------------

void probsort (
        const int n, 
        std::vector<trap_animal> &tran)
    // Sort using Shell algorithm see Press et al 1989 p 257
    // tran is an array of trap_animal records
    
{
    double aln2i = 1.442695022;
    double tiny  = 1.0e-5;
    int nn,m,lognb2,l,k,j,i;
    trap_animal t;
    lognb2 = trunc(log(n)*aln2i+tiny);
    m = n;
    for (nn=1; nn<=lognb2; nn++)
    {
        m = m / 2;
        k = n-m;
        for (j=1; j<=k; j++)
        {
            i = j;
            lab1:    l = i+m;
            if (tran[l-1].time < tran[i-1].time)
            {
                t = tran[i-1];
                tran[i-1] = tran[l-1];
                tran[l-1] = t;
                i = i-m;
                if (i >= 1)  goto lab1;
            }
        }
    }
}    // end of probsort 
//----------------------------------------------------------------

// random count from different distributions 

double rcount (int binomN, double lambda, const double Tsk) {
    
    // Poisson 
    if (binomN == 0)
        return (R::rpois(lambda * Tsk) );
    
    // negative binomial 
    else if (binomN < 0) {
        // must use 'size, prob' parameters 
        // prob = size / (size + mu) 
        binomN = abs(binomN);
        return (R::rnbinom(binomN, binomN / (binomN+ (lambda * Tsk))) );
    }
    
    else { 
        if (fabs(Tsk-1) > 1e-10)             
            lambda = 1 - pow(1-lambda, Tsk);   // 2012-12-18 
        
        // Bernoulli 
        if (binomN == 1) {
            if (unif_rand() < lambda)
                return (1);
            else
                return (0);
        }
        // binomial 
        else
            return (R::rbinom(binomN, lambda) );
    }
}
//----------------------------------------------------------------

// return probability g(r) for given detection function fn 
// used in simsecr.cpp and trapping.cpp 
double gr (
        const int fn,
        const Rcpp::NumericVector gsb,
        const rpoint xy,
        const rpoint animal) {
    double r;
    fnptrC fnp;
    fnp = getgfns(fn);
    r = distance1 (xy, animal);
    return (fnp(as<std::vector<double>>(gsb),r));
}
//----------------------------------------------------------------


double hazard (double pp) {
    if (pp > (1-fuzz))  // pp close to 1.0 - approx limit 
        pp = huge;      // g0 very large (effecti inf hazard) 
    else {
        if (pp <= 0) 
            pp = 0;
        else 
            pp = -log(1-pp);
    }
    return(pp);
}
//=============================================================

// detect may take values -
// 0  multi-catch traps
// 1  binary proximity detectors
// 2  count  proximity detectors
// 3  exclusive polygon detector
// 4  exclusive transect detector
// 5  signal detector
// 6  polygon detector
// 7  transect detector
// 8  times  (undocumented)
// 9  cue    (undocumented) -- removed in secr 2.10.0
// 12 signalnoise


//--------------------------------------------------------------------------

// Functions to characterize detector type 
// polygon, transect and signal detector types must be constant 
// across occasions 

bool anyexclusive (const IntegerVector detect) {
    bool exclusive = false;
    for (int s=0; s< detect.size(); s++) {
        if ((detect[s]==0) || (detect[s]==3) || (detect[s]==4))
            exclusive = true;
    }
    return exclusive;
}

bool anycapped  (const IntegerVector detect) {
    bool capped = false;
    for (int s=0; s<detect.size(); s++) {
        if (detect[s]==8)
            capped = true;
    }
    return capped;
}

bool anypolygon  (const IntegerVector detect) {
    bool polygon = false;
    for (int s=0; s<detect.size(); s++) {
        if ((detect[s]==3) || (detect[s]==6) )
            polygon = true;
    }
    return polygon;
}

bool anytransect (const IntegerVector detect) {
    bool transect = false;
    for (int s=0; s<detect.size(); s++) {
        if ((detect[s]==4) || (detect[s]==7))
            transect = true;
    }
    return transect;
}

bool anysignal (const IntegerVector detect) {
    bool signal = false;
    for (int s=0; s<detect.size(); s++) {
        if ((detect[s]==5) || (detect[s]==12))
            signal = true;
    }
    return signal;
}

bool anytelemetry (const IntegerVector detect) {
    bool telemetry = false;
    for (int s=0; s<detect.size(); s++) {
        if (detect[s]==13)
            telemetry = true;
    }
    return telemetry;
}

//  check if we need to consider variation among individuals 
// i.e. check if detection parameters constant for given s,k 
bool anyvarying (
        const int    nc,     // number of capture histories (or groups if PIA0 has that dim) 
        const int    ss,     // number of occasions 
        const int    nk,     // number of traps 
        const int    nmix,   // number of mixture classes 
        const IntegerVector &PIA0  // lookup which g0/sigma/b combination to use for given n, S, K [naive] 
) {
    int i,n,s,k,x;
    int wxi;
    bool indiv = false;
    for (s=0; s<ss; s++) {
        for (k=0; k<nk; k++) {
            for (x=0; x<nmix; x++) {
                wxi = i4(0,s,k,x,nc,ss,nk);       
                i = PIA0[wxi];
                for (n=1; n<nc; n++) {
                    wxi = i4(n,s,k,x,nc,ss,nk);    
                    if (i != PIA0[wxi]) {
                        indiv = true; break;
                    }
                }
            }
        }
    }
    return(indiv);
}
//--------------------------------------------------------------------

bool alltelemetry (const IntegerVector detect) {
    bool telemetry = true;
    for (int s=0; s<detect.size(); s++) {
        if ((detect[s]!=13))
            telemetry = false;
    }
    return telemetry;
}

bool allpoint (const IntegerVector detect, bool allowsignal, bool allowtelem) {
    bool point;
    bool OK = true;
    for (int s=0; s<detect.size(); s++) {
        point = (detect[s]==0) || (detect[s]==1) || (detect[s]==2) || detect[s] == 8
        || (detect[s]==10) || (detect[s]==11)
        || (allowsignal && ((detect[s]==5) || (detect[s]==12)))
        || (allowtelem && ((detect[s]==13)));
        OK = OK && point;
    }
    return OK;
}

bool allcapped  (const IntegerVector detect) {
    bool OK = true;
    for (int s=0; s<detect.size(); s++) {
        OK = OK && (detect[s] == 8);
    }
    return OK;
}

bool allmulti (const IntegerVector detect) {
    bool notmulti = false;
    for (int s=0; s<detect.size(); s++) {
        if (detect[s]!=0)
            notmulti = true;
    }
    return (!notmulti);
}
//--------------------------------------------------------------------

// Do parameter values for naive animals differ at all from those for other animals ?
bool anyb (const NumericMatrix &gsbval, const NumericMatrix &gsb0val) {
    bool identical = true;
    for (int i=0; i<gsbval.size(); i++) {
        if (gsbval[i] != gsb0val[i]) identical = false;
    }
    return (!identical);
}
//==============================================================================
// 
// std::vector<int> fillcumkcpp(
//         const IntegerVector detect, 
//         const int ss, 
//         const IntegerVector kk)
// {
//     // determine number of polygons if polygon detector 
//     // for polygon detectors, kk is vector ending in zero
//     // and first nk+1 elements of vector cumk are filled 
//     int i;
//     int nk = 0;
//     std::vector<int> cumk(maxnpoly);
//     if (anypolygon(detect) || anytransect(detect)) {
//         cumk[0] = 0;
//         for (i=0; i<maxnpoly; i++) {
//             if (kk[i]<=0) break;
//             cumk[i+1] = cumk[i] + kk[i];
//             nk++;
//         }
//         for (i=0; i<nk; i++)                // over parts 
//             if ((cumk[i+1] - cumk[i]) > maxvertices)
//                 stop("exceeded maximum number of vertices %d per polygon", maxvertices);
//     }
//     else
//         nk = kk[0];  // return number of detectors unchanged 
//     cumk.resize(nk);
//     return (cumk);
// }

//--------------------------------------------------------------------

// int firstkcpp (const int n,
//                const int x,  
//                const int nc, 
//                const int ss, 
//                const int nk,  
//                const IntegerVector &PIA) {
//     // return index of first detector for which PIA is non-zero 
//     int wxi;
//     int k=-1;
//     do {
//         k++;
//         wxi = i4(n,0,k,x,nc,ss,nk);
//     }
//     while ((PIA[wxi] == 0) && (k<nk));
//     if (k>=nk) return (NAN);  // Rcpp::stop ("no detector used on first occasion? error in getpmix"); 
//     else return(k);
// }

//=============================================================

// void getpmixcpp(int gpar, 
//                 const int nc1, 
//                 const int nmix, 
//                 const IntegerVector &knownclass, 
//                 const int nc, 
//                 const int cc, 
//                 const int ss, 
//                 const int nk, 
//                 const IntegerVector &grp, 
//                 const IntegerVector &PIA, 
//                 const NumericMatrix &gsbval,
//                 std::vector<double> &pmixg, 
//                 std::vector<double> &pmixn) {
//     
//     //-------------------------------
//     // mixture proportions           
//     // by group and by animal        
//     
//     int g, n, c, x, wxi;
//     double pmix;
//     
//     if (nmix>1) {
//         // one extra real parameter for h2, h3 
//         gpar++;
//         for (n=0; n<nc1; n++) {
//             for (x=0; x<nmix; x++) {
//                 wxi = i4(n,0,firstkcpp(n,x,nc,ss,nk,PIA),x,nc,ss,nk);
//                 c = PIA[wxi] - 1;
//                 if (c<0) Rcpp::stop ("c<0 error in getpmix"); 
//                 pmix = gsbval[cc * (gpar-1) + c];  // last column in gsbval 
//                 
//                 // group-specific, and overall pmix by class for knownclass case 
//                 g = grp[n]-1;
//                 pmixg[nmix * g + x] = pmix;
//                 
//                 // individual-specific                 
//                 if (knownclass[n] > 1) {
//                     if (knownclass[n] == (x+2))   // knownclass=2 maps to x=0 
//                         pmixn[nmix * n + x] = 1;
//                     else 
//                         pmixn[nmix * n + x] = 0;
//                 }
//                 else
//                     pmixn[nmix * n + x] = pmix;
//             }
//         }
//     }
// }
//=============================================================

void fillngcpp(
        const int nc, 
        const int gg, 
        const IntegerVector &grp, 
        std::vector<int> &ng) {
    int g, n;
    // Count number per group (not used for CL)                
    // Assume histories sorted by group = individual           
    // CH are numbered 0 <= n < nc in C code                  
    for (g=0; g<gg; g++)
        ng[g] = 0;
    for (n=0; n<nc; n++) { 
        g = grp[n] - 1; 
        ng[g]++;
    }
}

//=============================================================

int nval(int detect0, int nc1, int cc, int ss, int nk) {
    // compute and allocate the space needed for ancillary data 
    int nv;
    if ((detect0==3) || (detect0==4))
        nv = 2 + cc + nc1 * ss;
    else if ((detect0==6) || (detect0==7))
        nv = 2 + cc + nc1 * ss * nk;
    else if (detect0==5)
        nv = 4 + nc1 * ss * nk;
    else if (detect0==12)    // signalnoise 
        nv = 6 + nc1 * ss * nk;
    else
        nv = ss + nc1;    // default for (possibly mixed) 0,1,2,8,13 
    return(nv);
}
//=============================================================

NumericMatrix makedist2cpp (
        const NumericMatrix &traps, 
        const NumericMatrix &mask) 
{
    int kk, mm;
    int k, m;
    kk = traps.nrow();
    mm = mask.nrow();
    NumericMatrix dist2(kk,mm);
    // fill array with squared Euclidean distances 
    for (k=0; k<kk; k++)
        for (m=0; m<mm; m++) 
            dist2(k,m) = d2cpp(k, m, traps, mask);
    return(dist2);
}
//--------------------------------------------------------------------------

void squaredistcpp (
        NumericMatrix &dist2)
{
    int i;
    // just square the inputs 
    for (i=0; i<dist2.size(); i++)
        dist2[i] = dist2[i] * dist2[i];
}

//--------------------------------------------------------------------------

// int getstart(
//         const IntegerVector &detect, 
//         std::vector<int> &start, 
//         const int nc1, 
//         const int nc, 
//         const int ss, 
//         const int nk, 
//         const IntegerVector &w) {
//     
//     //------------------------------------------------------------
//     // Identify start positions of ancillary data for each animal 
//     // Telemetry is incompatible with polygon or signal detectors 
//     
//     int nd = 0;
//     int i,s,k,wi,first;
//     
//     if (anytelemetry(detect)) {
//         // assume all telemetry fixes associate with last detector 
//         for (i=0; i< nc; i++) {
//             first = 1;
//             for (s=0; s<ss; s++) {
//                 if (first) start[i] = nd;
//                 first = 0;
//                 wi = i3(i,s,nk-1,nc,ss);
//                 nd += abs(w[wi]);
//             }
//         }
//     }
//     else {
//         if (anypolygon(detect) || anytransect(detect) || 
//             anysignal(detect)) {
//             // start[z] indexes the first row in xy (or element in signal)
//             //   for each possible count z, where z is w-order (isk) 
//             for (k=0; k<nk; k++) {
//                 for (s=0; s< ss; s++) {
//                     for (i=0; i< nc; i++) {
//                         wi = i3(i,s,k,nc,ss);
//                         start[wi] = nd;
//                         nd += abs(w[wi]);
//                     }
//                 }
//             }
//         }
//     }
//     return(nd);
// }
//=============================================================

void getdetspec (
        const IntegerVector &detect, 
        const int fn, 
        const int nc,  
        const int nc1, 
        const int cc, 
        const int nmix, 
        const int nd, 
        const int nk, 
        const int ss, 
        const int mm, 
        const IntegerVector &PIA, 
        const NumericVector &miscparm, 
        const std::vector<int> &start, 
        std::vector<double> &detspec) {
    
    // detector-type-specific data passed later to prwi functions 
    
    // mixtures are group-specific for full likelihood, and     
    // individual-specific for conditional likelihood           
    
    int i,s;
    
    // default for (possibly mixed) point detectors and telemetry 
    if (allpoint(detect, 0, 1)) {
        for(s=0; s<ss; s++)
            detspec[s] = (double) detect[s];
        
        // start position (in xy) of telemetry fixes of each animal 
        if (anytelemetry(detect)) {
            for (i=0; i< nc; i++)
                detspec[ss+i] = (double) start[i];
        }
    }
    //  polygonX and transectX 
    else if ((detect[0] == 3) || (detect[0] == 4)) {
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        // pass index of first detection of each animal + occasion 
        // maximum 1 per polygon as exclusive 
        for (i=0; i< (nc * ss); i++)
            detspec[2+cc+i] = (double) start[i];
    }
    //  polygon and transect 
    else if ((detect[0] == 6) || (detect[0] == 7)) {
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        // pass index of first detection of each animal + occasion 
        // maximum 1 per polygon as exclusive 
        for (i=0; i< (nc * ss * nk); i++)
            detspec[2+cc+i] = (double) start[i];
    }
    // signal  
    else if (detect[0] == 5) {    
        for (i=0; i<3; i++) detspec[i]= miscparm[i];
        detspec[3]= ((fn == 11) || (fn == 13));     // spherical 
        for (i=0; i< (nc * ss * nk); i++)
            detspec[4+i] = (double) start[i];
    }
    // signal-noise 
    else if (detect[0] == 12) {                           
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;                     // number of detections or detectors??
        detspec[2]= miscparm[0];                      // cut 
        detspec[3]= miscparm[1];                      // noise mean 
        detspec[4]= miscparm[2];                      // noise sd 
        detspec[5]= ((fn == 11) || (fn == 13));     // spherical 
        for (i=0; i< (nc * ss * nk); i++)
            detspec[6+i] = (double) start[i];
    }
}
//=============================================================

void geth2 (
        const int nc1, 
        const int cc, 
        const int nmix, 
        const int mm, 
        const IntegerVector &PIA, 
        const std::vector<double> &hk, 
        const NumericMatrix &Tsk, 
        std::vector<double> &h, 
        std::vector<int> &hindex) 
    
    // This function fills a vector h representing a 4-D (x,m,n,s) array with
    // the total hazard (summed across traps) for animal n on occasion s 
    // wrt mask point m and latent class x
    
    // Computation is limited to combinations of n, s with unique parameter combinations 
    // (values in PIA and Tsk) and the returned n x s matrix 'hindex' contains the index for
    // each n, s to the unique total in h (for given x, m).
    
    // adapted from gethcpp in openCR 2018-11-09
    // replaces defective geth
    
{
    int c,i,m,n,k,x,gi,hi,s;
    int nk = Tsk.nrow();
    int ss = Tsk.ncol(); 
    double Tski;          
    int row, col;
    int nrow, ncol;
    int uniquerows;
    List lookuplist;
    nrow = nc1*ss;
    ncol = nk*(nmix+1);
    NumericMatrix xmat(nrow, ncol);
    
    //---------------------------------------
    // find unique combinations
    for (n=0; n<nc1; n++) {
        for (s=0; s<ss; s++) {
            row = nc1*s+n;
            for (k=0; k<nk; k++) {                         
                for(x=0; x<nmix; x++) {
                    col = x * nk + k;
                    xmat(row,col) = PIA[i4(n,s,k,x, nc1, ss, nk)];
                }
                col = nmix * nk + k;  // append usage for this trap and occasion
                xmat(row,col) = Tsk(k,s);		
            }
        }
    }
    
    lookuplist = makelookupcpp(xmat);
    uniquerows = lookuplist["uniquerows"];
    hindex = as <std::vector<int> >(lookuplist["index"]);
    
    // need zero-based index 
    for (i=0; i< nrow; i++) hindex[i]--;
    
    // zero array for accumulated hazard h
    for (i=0; i<(uniquerows * mm * nmix); i++) h[i] = 0;
    
    // search hindex for each row index in turn, identifying first n,s with the index
    // fill h[] for this row
    hi = 0;
    for (s=0; s < ss; s++) {    // scan by occasion as new hi appear in column order
        for (n=0; n < nc1; n++) {
            if (hindex[s*nc1 + n] == hi) {
                for (k=0; k < nk; k++) {
                    Tski = Tsk[s * nk + k];
                    for (x=0; x<nmix; x++) {
                        c = PIA[i4(n,s,k,x, nc1, ss, nk)]-1;
                        // c<0 (PIA=0) implies detector not used on this occasion
                        if (c >= 0) {
                            for (m=0; m<mm; m++) {
                                gi = i3(c,k,m,cc,nk);
                                h[i3(x,m,hi, nmix, mm)] += Tski * hk[gi];
                            }
                        }
                    }
                }
                hi++;
            } 
            if (hi >= uniquerows) break;
        }
        if (hi >= uniquerows) break;
    }   
}
//=============================================================

// probability of count for session s, detector k, animal i
// The argument 'g' is understood to be a cumulative hazard if binomN=0,
// a probability otherwise

double pski ( int binomN,
              int count,
              double Tski,
              double g,
              double pI) {
    
    double lambda;
    double result = 1.0;
    
    if (binomN == -1) {                              // binary proximity detectors : Bernoulli
        if (abs(Tski-1) > 1e-10) {                   // effort not unity; adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        if (count>0)                                 
            result = g*pI;  
        else 
            result = 1 - g*pI;
    }
    else if (binomN == 0) {                          // count detectors : Poisson 
        lambda = Tski * g * pI;
        if ((count < 0) || (count>0 && lambda<=0)) {         
            result = 0;
        }
        else if (count == 0) {
            result = exp(-lambda);            // routinely apply Tsk adjustment to cum. hazard 
        }
        else {
            boost::math::poisson_distribution<> pois(lambda);
            result = boost::math::pdf(pois,count);
        }
    }
    else if (binomN == 1) {                          // count detectors : Binomial, size from Tsk
        result = gbinom (count, round(Tski), g*pI); 
    }
    else if (binomN > 1) {                           // count detectors : Binomial, specified size 
        if (abs(Tski-1) > 1e-10) {                   // effort not unity, adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        result = gbinom (count, binomN, g*pI);
    }
    else result = NAN; // Rcpp::stop("binomN < -1 not allowed");  // code multi -2 separately
    
    return (result);
}
//--------------------------------------------------------------------------

double classmembership (
        const int n, 
        const int x, 
        const IntegerVector &knownclass, 
        const std::vector<double> &pmixn, 
        const int nmix) {
    
    // Return probability individual n belongs to class x. This may be binary 
    //   (0/1) in the case of known class, or continuous if class is unknown 
    
    double pmixnx = 0;
    int knownx = -1;
    
    if (knownclass[n] == 1) 
        knownx = -1;                         // unknown class 
    else
        knownx = R::imax2(0, knownclass[n]-2);  // known class 
    
    // unknown class : weighted by probability of membership  
    // known class does not require probability of membership 
    if (knownx < 0)
        pmixnx = pmixn[nmix * n + x];
    else if (knownx == x)
        pmixnx = 1.0;
    else 
        pmixnx = 0.0;
    return (pmixnx);
    
}
//--------------------------------------------------------------------------

