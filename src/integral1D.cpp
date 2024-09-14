#include "poly.h"

double hintegral1 (int fn, double par[]) {
    /* integral of 1-D function */
        /* from 2017-03-22 this is strictly a hazard function */
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
    
    /* treat uniform separately */
        if (fn == 4) {
            return (2 * par[0] * par[1]);  /* symmetric about 0 */
        }
    
    a = 0;
    b = 1;    /* signals bounds 0,Inf */
        ex[0] = par[0];
    ex[1] = par[1];
    ex[2] = par[2];
    ex[3] = fn;
    
    Rdqagi(justgr, &ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, 
           &ier, &limit, &lenw, &last, iwork, work);
    
    /* ignoring ier etc. */
    return (result * 2);    /* symmetric about 0 */
}
/*----------------------------------------------------------------*/

void fx1 (double *x, int n, void *ex) {
    int i;
    int ns;
    int fn;
    /*    struct rpoint *line;*/
    struct rpoint line[maxvertices * 2];
    struct rpoint mxy;
    struct rpoint xy;
    double * p;
    double cumd[maxvertices * 2];
    double d;
    fnptr fnzr = zhnr;
    /* extract parameters passed in void pointer ex */
    p = (double*) ex;
    fn = round(p[3]);
    mxy.x = p[4];
    mxy.y = p[5];
    ns = round(p[9]);
    /* coordinates of vertices */
    /*    line = (struct rpoint *) R_alloc(ns, sizeof(struct rpoint));*/
    for (i=0; i<ns; i++) {
        line[i].x = p[i+10];
        line[i].y = p[i+ns+10];
    }
    /* cumulative distance along line */
    /* cumd = (double *) R_alloc(ns + 1, sizeof(double)); */
    cumd[0] = 0;
    for (i=0; i<(ns-1); i++) {
        cumd[i+1] = cumd[i] + distance (line[i],line[i+1]);
    }
    /* set detection function - default zhnr */
    fnzr = getzfnr(fn);  // 2016-01-02, 2017-03-22
    /* for each x in x[] */
    for (i=0; i<n; i++) {
        xy = getxy (x[i], cumd, line, ns, 0);
        d = distance (xy, mxy);
        x[i] = fnzr(p, d);   /* z(r) */
    }
}

double integral1D
(int fn, int m, int c, double gsbval[], int cc, double traps[],
 double mask[], int n1, int n2, int kk, int mm, double ex[])
{
    double ax=0;
    double bx=0;
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
    int k;
    int ns;
    ns = n2-n1+1;
    
    /* 2011-06-21 uniform treated separately */
    if (fn == 4) {
        for (k=n1+1; k<=n2; k++) {  /* upper bound is length of this transect */
    bx += SegCircle2(
        traps[k-1], traps[k-1+kk],
                         traps[k], traps[k+kk],
                                        mask[m], mask[m+mm],
                                                     gsbval[cc + c]);
        }
        return (bx);
    }
    
    for (k=n1+1; k<=n2; k++) {  /* upper bound is length of this transect */
    bx += sqrt( (traps[k] - traps[k-1]) * (traps[k] - traps[k-1]) +
        (traps[k+kk] - traps[k-1+kk]) * (traps[k+kk] - traps[k-1+kk]) );
    }
    /* pass parameters etc. through pointer */
    ex[0] = gsbval[c];
    ex[1] = gsbval[cc + c];
    ex[2] = gsbval[2*cc + c];
    ex[3] = fn;
    ex[4] = mask[m];
    ex[5] = mask[m+mm];
    ex[6] = 0;
    ex[7] = 0;
    ex[8] = 0;
    ex[9] = ns;
    for (k=0; k<ns; k++) {             /* pass transect vertices */
    ex[k+10] = traps[k+n1];        /* x */
    ex[k+ns+10] = traps[k+n1+kk];  /* y */
    }
    Rdqags(fx1, ex, &ax, &bx, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
           &limit, &lenw, &last, iwork, work);
    if (ier != 0) Rprintf("ier error code in integral1D %5d\n", ier);
    return (result);
}

/*===============================================================*/
