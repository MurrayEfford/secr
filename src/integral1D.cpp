#include "poly.h"

#include <R_ext/Utils.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

// The 1-D integration code here is provided for transect simulation in trapping.cpp.
// It does not depend on RcppNumerical, which will lighten the package as a whole 
// once polygonN.cpp is removed to the possible new package 'secrpoly'.

// Calculate the length of intersection of a line segment and a circle
// Based on C code of Paul Bourke November 1992
// Line segment is defined from p1 to p2
// Circle is of radius r and centred at sc
// Two potential points of intersection given by
// p = p1 + mu1 (p2-p1)
// p = p1 + mu2 (p2-p1)
// Return 0 if line segment does not intersect circle

//----------------------------------------------------------------
// integral of 1-D function
// revived old code 2024-09-15
// this code does not depend on RcppNumerical and RcppEigen

// fnptr defined detectfn.cpp, secr.h
// zhnr  defined detectfn.cpp, arguments NumericVector, double
// getzfnr(fn) defined detectfn.cpp, secr.h returns fnptr

void justgr(double *x, int n, void *ex) {
    int i;
    int fn;
    double * p;
    p = (double*) ex;
    NumericVector tmp(4);
    for (i=0; i<4; i++) tmp(i) = p[i];
    fn = tmp(3);
    fnptr fnp;
    fnp = getzfnr(fn);   // fnp expects numeric vector
    for (i=0; i<n; i++) {
        x[i] = fnp(tmp, x[i]);
    }
}

double hintegral1 (int fn, const NumericVector par) {
    double ex[4];
    double a;
    int b;
    double epsabs = 0.0001;
    double epsrel = 0.0001;
    double result = 0;
    double abserr = 0;
    int neval = 0;
    int ier = 0;
    int limit = 100;
    int lenw = 400;
    int last = 0;
    int iwork[100];
    double work[400];
    
    // treat uniform separately
    if (fn == 4) {
        return (2 * par(0) * par(1));  // symmetric about 0
    }
    
    a = 0;
    b = 1;    // signals bounds 0,Inf
    ex[0] = par(0);
    ex[1] = par(1);
    ex[2] = par(2);
    ex[3] = fn;
    
    Rdqagi(justgr, &ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval,
           &ier, &limit, &lenw, &last, iwork, work);
    
    // ignoring ier etc.
    return (result * 2);    // symmetric about 0
}

//----------------------------------------------------------------

void fx1 (double *x, int n, void *ex) {
    int i;
    int ns;
    int fn;
    struct rpoint line[maxvertices * 2];
    struct rpoint mxy;
    struct rpoint xy;
    double * p;
    double cumd[maxvertices * 2];
    double d;
    fnptr fnzr;
    
    // extract parameters passed in void pointer ex
    p = (double*) ex;
    fn = std::lround(p[3]);
    // set detection function - default zhnr
    fnzr = getzfnr(fn);
    
    mxy.x = p[4];
    mxy.y = p[5];
    ns = std::lround(p[9]);
    // coordinates of vertices
    //    line = (struct rpoint *) R_alloc(ns, sizeof(struct rpoint));
    for (i=0; i<ns; i++) {
        line[i].x = p[i+10];
        line[i].y = p[i+ns+10];
    }
    // cumulative distance along line
    // cumd = (double *) R_alloc(ns + 1, sizeof(double));
    cumd[0] = 0;
    for (i=0; i<(ns-1); i++) {
        cumd[i+1] = cumd[i] + distance1 (line[i],line[i+1]);
    }
    NumericVector gsb(4);
    gsb(0) = p[0];
    gsb(1) = p[1];
    gsb(2) = p[2];

    // for each x in x[]
    for (i=0; i<n; i++) {
        xy = getxy (x[i], cumd, line, ns, 0);
        d = distance1 (xy, mxy);
        x[i] = fnzr(gsb, d);   // z(r)
    }
}

double integral1D (int fn, int m, int c,
                   NumericVector gsbval,
                   int cc,
                   NumericMatrix traps,
                   NumericMatrix mask,
                   int n1,
                   int n2,
                   int kk,
                   int mm)
{
    double result = 0;
    double ax=0;
    double bx=0;
    double epsabs = 0.0001;
    double epsrel = 0.0001;
    double abserr = 0;
    int    neval = 0;
    int    ier = 0;
    int    limit = 100;
    int    lenw;
    int    last = 0;
    int    k;
    int    ns = n2-n1+1;
    double *ex;
    int    *iwork;
    double *work;
    lenw  = 4 * limit;
    ex    = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
    iwork = (int *)    R_alloc(limit, sizeof(int));
    work  = (double *) R_alloc(lenw,  sizeof(double));

    // uniform treated separately
    if (fn == 4) {
        // upper bound is length of this transect
        for (k=n1+1; k<=n2; k++) {
            bx += SegCircle2(
                traps(k-1,0), traps(k-1,1),
                traps(k,0), traps(k,1),
                mask(m,0), mask(m,1),
                gsbval(cc + c));
        }
        return (bx);
    }
    // upper bound is length of this transect
    for (k=n1+1; k<=n2; k++) {
        bx += std::sqrt( (traps(k,0) - traps(k-1,0)) * (traps(k,0) - traps(k-1,0)) +
            (traps(k,1) - traps(k-1,1)) * (traps(k,1) - traps(k-1,1)) );
    }

    // pass parameters etc. through pointer
    ex[0] = gsbval(c);
    ex[1] = gsbval(cc + c);
    ex[2] = gsbval(2*cc + c);
    ex[3] = fn;
    ex[4] = mask(m,0);
    ex[5] = mask(m,1);
    ex[6] = 0;
    ex[7] = 0;
    ex[8] = 0;
    ex[9] = ns;
    // pass transect vertices
    for (k=0; k<ns; k++) {
        ex[k+10] = traps(k+n1,0);     // x
        ex[k+ns+10] = traps(k+n1,1);  // y
    }
    Rdqags(fx1, ex, &ax, &bx, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
           &limit, &lenw, &last, iwork, work);
    if (ier != 0) Rprintf("ier error code in integral1D %5d\n", ier);
    return (result);
}

//===============================================================
