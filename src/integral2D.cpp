#include "poly.h"

#include <R_ext/Utils.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

// C++ code for R package 'secr' 
// Murray Efford 

// //---------------------------------------------------------------------------------
// // original 2-D code using Rdqags
// //---------------------------------------------------------------------------------
// 
// // find upper and lower points on perimeter of poly at x-coordinate x
// void yab(double x[], int *i, int *np, double poly[], double *a, double *b) {
//     int k;
//     int nv = 0;
//     double ab[3];
//     // note 'sign' is RMath function
//     for (k=0; k< (*np-1); k++) {
//         if (R::sign(poly[k]- x[*i]) != R::sign(poly[k+1]- x[*i])) {
//             ab[nv] = poly[k+ *np] + (x[*i]-poly[k]) * (poly[k+1+*np]-poly[k+*np]) /
//                 (poly[k+1]-poly[k]);
//             nv++;
//         }
//         if (nv>2) break;
//     }
//     if (ab[0]>ab[1])
//     {*a = ab[1]; *b = ab[0];}
//     else
//     {*a = ab[0]; *b = ab[1];}
// 
// }
// 
// void fy(double *x, int n, void *ex) {
//     int i;
//     int fn;
//     double * p;
//     double mx,my;
//     double xy[2];
//     double d;
//     fnptr fnzr; // = zhnr;
//     NumericVector gsb(3);   // passed to fnzr
//     p = (double*) ex;
//     gsb[0] = p[0];
//     gsb[1] = p[1];
//     gsb[2] = p[2];
//     fn = std::round(p[3]);
//     mx = p[4];
//     my = p[5];
//     xy[0] = p[6];
//     // set detection function
//     fnzr = getzfnr(fn);  // 2016-01-02, 2017-03-22
//     for (i=0; i<n; i++) {
//         xy[1] = x[i];   // set each y value
//         d = std::sqrt ( (xy[1]-my)*(xy[1]-my) + (xy[0]-mx)*(xy[0]-mx) );
//         x[i] = fnzr(gsb, d);   // z(r)
//     }
// }
// 
// void fx(double *x, int n, void *ex) {
//     int i;
//     double * p;
// 
//     double poly[maxvertices * 2];
//     double a;
//     double b;
// 
//     double epsabs = 0.0001;
//     double epsrel = 0.0001;
//     double result = 0;
//     double abserr = 0;
//     int neval = 0;
//     int ier = 0;
//     int limit = 100;
//     int lenw = 400;
//     int last = 0;
//     int iwork[100];
//     double work[400];
//     int kk;
//     p = (double*) ex;
//     kk = std::round(p[9]);
// 
//     for (i=0; i<kk; i++) {
//         poly[i] = p[i+10];
//         poly[i+kk] = p[i+kk+10];
//     }
//     for (i=0; i<n; i++) {
//         yab(x, &i, &kk, poly, &a, &b);  // refine limits here
//         p[6] = x[i];                    // pass forward the value of x; consider &ex etc.
//         Rdqags(fy, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
//                &limit, &lenw, &last, iwork, work);
//         x[i] = result;
//     }
// }
// 
// double integral2Dcpp  (const int &fn,
//                        const int &m,
//                        const int &c,
//                        const RcppParallel::RMatrix<double> &gsbval,
//                        const RcppParallel::RMatrix<double> &traps,
//                        const RcppParallel::RMatrix<double> &mask,
//                        const int &n1,
//                        const int &n2,
//                        double ex[]) {
//     double ax=1e20;
//     double bx=-1e20;
//     double epsabs = 0.0001;
//     double epsrel = 0.0001;
//     double result = 0;
//     double abserr = 0;
//     int neval = 0;
//     int ier = 0;
//     int limit = 100;
//     int lenw = 400;
//     int last = 0;
//     int iwork[100];
//     double work[400];
//     int k;
//     int ns;
//     int reportier = 0;
// 
//     // limits from bounding box of this polygon
//     ns = n2-n1+1;
//     for (k=0; k<ns; k++) {
//         ax = std::min(ax, traps[k+n1]);
//         bx = std::max(bx, traps[k+n1]);
//     }
// 
//     // pass parameters etc. through pointer
//     ex[0] = gsbval(c,0);
//     ex[1] = gsbval(c,1);
//     ex[2] = gsbval(c,2);
//     ex[3] = fn;
//     ex[4] = mask(m,0);
//     ex[5] = mask(m,1);
//     ex[6] = 0;
//     ex[7] = 0;
//     ex[8] = 0;
//     ex[9] = ns;
// 
//     // also pass polygon vertices
//     for (k=0; k<ns; k++) {
//         ex[k+10] = traps(k+n1,0);        // x
//         ex[k+ns+10] = traps(k+n1,1);       // y
//     }
//     Rdqags(fx, ex, &ax, &bx, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
//            &limit, &lenw, &last, iwork, work);
//     if ((ier != 0) & (reportier))
//         Rprintf("ier error code in integral2Dcpp %5d\n", ier);
//     return (result);
// }
// // end original 2-D code using Rdqags
// //---------------------------------------------------------------------------------
// 
