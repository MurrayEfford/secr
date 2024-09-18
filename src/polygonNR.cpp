#include "poly.h"

using namespace Rcpp;

// C++ code for R package 'secr' 
// Murray Efford 

// This file contains functions for integration of home range overlap
// with polygon detectors and similar.

// Uses RcppNumerical = Numer

//------------------------------------------------------------------------------
// original 1-D code using Rdqags moved to integral1D.cpp 2024-09-19
// original 2-D code using Rdqags moved to integral2D.cpp 2024-09-19
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// 1-D code using Numer::integrate Func
//------------------------------------------------------------------------------

// rpoint defined in secr.h

class fx1func: public Numer::Func
{
private:
    std::vector<double> gsb;      
    RcppParallel::RMatrix<double> line;
    int n1,n2;
    int fn;
    double mx, my;
    std::vector<double> cumd; 
    fnptrC fnzr;
    
public:
    
    fx1func(
        const std::vector<double> &gsb,
        const RcppParallel::RMatrix<double> &line,
        const int n1,
        const int n2,
        const int fn,
        const double mx,
        const double my,
        const std::vector<double> &cumd
    ) 
        : gsb(gsb), line(line), n1(n1), n2(n2), fn(fn), mx(mx), my(my), cumd(cumd) 
    {
        fnzr = getzfnrC(fn);   // set detection function 
    }
    
    
    double operator()(const double& x) const
    {
        double d;
        rpoint xy;
        xy = getxycpp (x, cumd, line, n1, n2);  // see transectutils.cpp
        d = std::sqrt (std::pow(xy.x-mx,2) + std::pow(xy.y-my,2));
        return(fnzr(gsb, d));   // z(r) 
    }
};

double integral1DNRcpp (
        const int fn, 
        const int m, 
        const int c, 
        const RcppParallel::RMatrix<double> &gsbval, 
        const RcppParallel::RMatrix<double> &traps,
        const RcppParallel::RMatrix<double> &mask, 
        const int n1, 
        const int n2)
{
    double err_est;
    int err_code;
    double ax=0;
    double bx=0;
    double dx;
    int j;
    std::vector<double> cumd(n2-n1+1,0);
    
    if (gsbval.ncol()>4) Rcpp::stop("bad gsbval matrix");
    std::vector<double> gsb(4);
    
    if (n2<=n1) return (0);
    
    // uniform treated separately 
    if (fn == 4) {
        for (j=n1+1; j<=n2; j++) {  // upper bound is length of this transect 
            bx += SegCircle2(traps(j-1,0), traps(j-1,1), traps(j,0), traps(j,1),
                             mask(m,0), mask(m,1), gsbval(c,1));
        }
        return (bx);
    }
    
    for (j=n1; j<n2; j++) {  
        dx = std::sqrt( 
            std::pow(traps(j+1,0) - traps(j,0), 2) +
                std::pow(traps(j+1,1) - traps(j,1), 2)
        );
        cumd[j-n1+1] = cumd[j-n1] + dx;
    }
    bx = cumd[n2-n1];  // length of this transect 
    
    int npar = gsbval.ncol();
    for (int i=0; i<npar; i++) gsb[i] = gsbval(c,i);
    
    fx1func f(gsb, traps, n1, n2, fn, mask(m,0), mask(m,1), cumd);
    const double res = Numer::integrate(f, ax, bx, err_est, err_code);
    return (res);
}

//---------------------------------------------------------------------------------
// 2-D code using repeated 1-D RcppNumerical Numer::integrate Func
//---------------------------------------------------------------------------------

class yslice: public Numer::Func
{
private:
    std::vector<double> gsb;
    int fn;
    double mx, my;
    double x;
    fnptrC fnzr; // = zhnr;

public:
    yslice(
        const std::vector<double> gsb,
        const int fn,
        const double mx,
        const double my,
        const double x)
        : gsb(gsb), fn(fn), mx(mx), my(my), x(x) {
        // set detection function
        fnzr = getzfnrC(fn);
    }

    double operator()(const double& y) const
    {
        double d;
        d = std::sqrt ((y-my)*(y-my) + (x-mx)*(x-mx));
        return(fnzr(gsb, d));   // z(r)
    }
};

class xfn: public Numer::Func {
    
private:
    std::vector<double> gsb;
    RcppParallel::RMatrix<double> poly;
    const int n1;
    const int n2;
    int fn;
    double mx;
    double my;
    fnptrC fnzr; // = zhnr;
    int np;
    
    class yslicei: public Numer::Func
    {
    private:
        const std::vector<double> gsb;
        const int fn;
              int n1;
              int n2;
        const double mx;
        const double my;
        fnptrC fnzr; // = zhnr;
        
    public:
        
        double x=0.0;
        
        yslicei(const std::vector<double> &gsb,
                const int &fn,
                const double &mx,
                const double &my): gsb(gsb), fn(fn), mx(mx), my(my)
        {fnzr = getzfnrC(fn);}
        
        double operator()(const double& y) const
        {
            double d;
            d = std::sqrt (std::pow(y-my,2) + std::pow(x-mx,2));
            return(fnzr(gsb, d));   // z(r)
        }
    };
    
public:
    
    xfn(const std::vector<double> &gsb,
        const RcppParallel::RMatrix<double> &poly,
        const int &n1,
        const int &n2,
        const int &fn,
        const double &mx,
        const double &my
    )
        : gsb(gsb), poly(poly), n1(n1), n2(n2), fn(fn), mx(mx), my(my) {
        // set detection function
        fnzr = getzfnrC(fn);
        np = poly.nrow();
    }
    
    
    // find upper and lower points on perimeter of poly at x-coordinate x
    std::vector<double> ylim(const double &x) const {
        int nv = 0;
        std::vector<double> ab(2);
        double tmp;
        for (int k=n1; k<n2; k++) {
            // note 'sign' is R function best avoided
            // if (R::sign(poly[k]- x) != R::sign(poly[k+1]- x)) {
            if (((poly[k]< x) && (poly[k+1] > x)) || ((poly[k]> x) && (poly[k+1] < x))) {
                ab[nv] = poly[k+ np] + (x-poly[k]) * (poly[k+1+np]-poly[k+np]) /
                    (poly[k+1]-poly[k]);
                nv++;
            }
            if (nv>2) break;
        }
        if (ab[0]>ab[1])
        {
            tmp = ab[1];
            ab[1] = ab[0];
            ab[0] = tmp;
        }
        return(ab);
    }
    
    double operator()(const double &x) const
    {
        double err_est;
        int err_code;
        std::vector<double> lim;
        yslicei f(gsb, fn, mx, my);
        f.x = x;
        lim = ylim(x);  // refine limits here
        const double res = Numer::integrate(f, lim[0], lim[1], err_est, err_code);
        return res;
    }
};

//---------------------------------------------------------------------------------------
// xfn2 is a version of xfn that uses pointwise insidecppC rather than merely integrating
// between upper and lower bounds in y dimension
//---------------------------------------------------------------------------------------

bool insidecppC (
        const Numer::Constvec &xy,
        const int    &n1,
        const int    &n2,
        const RcppParallel::RMatrix<double> &poly)
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
        if (std::abs(d)>0) { N = N/std::abs(d);  d = d/std::abs(d); }
        theta += std::atan2(N, d);
    }
    theta = std::abs(theta);
    return (std::abs(theta - 2* M_PI) < cutoff);    // M_PI is cmath.h constant 
}

class xfn2: public Numer::Func {
    
private:
    std::vector<double> gsb;
    RcppParallel::RMatrix<double> poly;
    int n1;
    int n2;
    int fn;
    double mx;
    double my;
    double ay;
    double by;
    
    fnptrC fnzr; // = zhnr;
    int np;
    
    class yslicei: public Numer::Func
    {
    private:
        const std::vector<double> gsb;
        const RcppParallel::RMatrix<double> &poly;
        const int n1;
        const int n2;
        const int fn;
        const double mx;
        const double my;
        const double ay;
        const double by;
        fnptrC fnzr; // = zhnr;
        
    public:
        
        double x=0.0;
        
        yslicei(const std::vector<double> &gsb,
                const RcppParallel::RMatrix<double> &poly,
                const int &n1,
                const int &n2,
                const int &fn,
                const double &mx,
                const double &my,
                const double &ay,
                const double &by)
            : gsb(gsb), poly(poly), n1(n1), n2(n2), fn(fn), mx(mx), my(my), ay(ay), by(by)
        {fnzr = getzfnrC(fn);}
        
        double operator()(const double& y) const
        {
            double d;
            Eigen::VectorXd xy(2);
            xy << x,y;
            if (insidecppC(xy, n1, n2, poly)) {
                d = std::sqrt (std::pow(y-my,2) + std::pow(x-mx,2));
                return(fnzr(gsb, d));   // z(r)
            }
            else {
                return(0);
            }
        }
    };
    
public:
    
    xfn2(const std::vector<double> &gsb,
        const RcppParallel::RMatrix<double> &poly,
        const int &n1,
        const int &n2,
        const int &fn,
        const double &mx,
        const double &my,
        const double &ay,
        const double &by
    )
        : gsb(gsb), poly(poly), n1(n1), n2(n2), fn(fn), mx(mx), my(my), ay(ay), by(by) {
        // set detection function
        fnzr = getzfnrC(fn);
        np = poly.nrow();
    }

    double operator()(const double &x) const
    {
        double err_est;
        int err_code;
        std::vector<double> lim;
        yslicei f(gsb, poly, n1, n2, fn, mx, my, ay, by);
        f.x = x;
        const double res = Numer::integrate(f, ay, by, err_est, err_code);
        return res;
    }
};
//---------------------------------------------------------------------------------------

// Apply either xfn or xfn2

double integral2DNRcpp (
        const int &fn,
        const int &m,
        const int &c,
        const RcppParallel::RMatrix<double> &gsbval,
        const RcppParallel::RMatrix<double> &poly,
        const RcppParallel::RMatrix<double> &mask,
        const int &n1,
        const int &n2,
        const bool &convex)
{
    double res;
    double err_est;
    int err_code;
    double ax=1e10;
    double bx=-1e10;
    double ay=1e10;    // used only with xfn2
    double by=-1e10;   // used only with xfn2
    int ns = n2-n1+1;
    for (int i=0; i<ns; i++) {
        ax = std::min(ax, poly(n1+i,0));
        bx = std::max(bx, poly(n1+i,0));
        ay = std::min(ay, poly(n1+i,1));    // used only with xfn2
        by = std::max(by, poly(n1+i,1));    // used only with xfn2
    }

    std::vector<double> gsb(4);
    int npar = gsbval.ncol();
    for (int i=0; i<npar; i++) gsb[i] = gsbval(c,i);

    // "Possible values are GaussKronrod{15, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 201}"
    // GaussKronrod41 default
    // GaussKronrod15 causes R failure
    // this variation seems slower
    // res = Numer::integrate(f, ax, bx, err_est, err_code, 100, 1e-6, 1e-5, Integrator<double>::GaussKronrod21);
    if (convex) {
        xfn f(gsb, poly, n1, n2, fn, mask(m,0), mask(m,1));
        res = Numer::integrate(f, ax, bx, err_est, err_code);
    }
    else {
        xfn2 f2(gsb, poly, n1, n2, fn, mask(m,0), mask(m,1), ay, by);
        res = Numer::integrate(f2, ax, bx, err_est, err_code);
    }
    return (res);
}

//===============================================================

//------------------------------------------------------------------------------
// alternative using 2-D code using RcppNumerical Numer::integrate MFunc
//------------------------------------------------------------------------------

// class fx2func: public Numer::MFunc
// {
// private:
//     std::vector<double> gsb;
//     RcppParallel::RMatrix<double> poly;
//     int n1;
//     int n2;
//     int fn;
//     double mx;
//     double my;
//     fnptrC fnzr;
//     
// public:
//     fx2func(const std::vector<double> &gsb,
//             const RcppParallel::RMatrix<double> &poly,
//             const int &n1,
//             const int &n2,
//             const int &fn,
//             const double &mx,
//             const double &my
//     ) 
//         : gsb(gsb), poly(poly), n1(n1), n2(n2), fn(fn), mx(mx), my(my) {
//         fnzr = getzfnrC(fn);    // set detection function 
//     }
//     
//     double operator()(Numer::Constvec &x) const
//     {
//         double d;
//         if (insidecppC(x, n1, n2, poly)) {
//             d = std::sqrt (std::pow(x[0]-mx,2) + std::pow(x[1]-my,2));
//             return(fnzr(gsb, d));   // z(r) 
//         }
//         else
//             return 0;
//     }
// };

// double integral2DNRcpp  (
//         const int &fn,
//         const int &m,
//         const int &c,
//         const RcppParallel::RMatrix<double> &gsbval,
//         const RcppParallel::RMatrix<double> &poly,
//         const RcppParallel::RMatrix<double> &mask,
//         const int &n1,
//         const int &n2) 
// {
//     Eigen::VectorXd lower(2);
//     Eigen::VectorXd upper(2);
//     double err_est;
//     int err_code;
//     int k;
//     int ns;
//     double xmin,ymin = 1e100;
//     double xmax,ymax = -1e100;
// 
//     // limits from bounding box of this polygon
//     ns = n2-n1+1;
//     for (k=0; k<ns; k++) {
//         xmin = std::min(xmin, poly(k+n1,0));  
//         ymin = std::min(ymin, poly(k+n1,1));  
//         xmax = std::max(xmax, poly(k+n1,0));  
//         ymax = std::max(ymax, poly(k+n1,1));  
//     }
//     lower[0] = xmin;
//     lower[1] = ymin;
//     upper[0] = xmax;
//     upper[1] = ymax;
//     
//     std::vector<double> gsb(4,0);
//     int npar = gsbval.ncol();
//     for (int i=0; i<npar; i++) gsb[i] = gsbval(c,i);
//     
//     fx2func f(gsb, poly, n1, n2, fn, mask(m,0), mask(m,1));
//     
//     const double res = Numer::integrate(f, lower, upper, err_est, err_code);
//     return (res);
// }
//===============================================================

// Miscellaneous 1-D integrations

class grfn: public Numer::Func
{
private:
    std::vector<double> gsb;
    int fn;
    fnptrC fnzr; // = zhnr;
    
public:
    grfn(const std::vector<double> gsb,
         const int fn) 
        : gsb(gsb), fn(fn) {
        fnzr = getzfnrC(fn);  // set detection function 
    }
    double operator()(const double& x) const
    {
        return fnzr(gsb, x);   // z(r) 
    }
};
//===============================================================

class rgrfn: public Numer::Func
{
private:
    std::vector<double> gsb;
    int fn;
    fnptrC fnzr; // = zhnr;
    
public:
    rgrfn(const std::vector<double> gsb,
          const int fn) 
        : gsb(gsb), fn(fn) {
        // set detection function 
        fnzr = getzfnrC(fn);  
    }
    double operator()(const double& r) const
    {
        return( r * fnzr(gsb, r));   // z(r) 
    }
};

//===============================================================

// integral of 1-D function (radial 2-D function)
double hintegral1DNRcpp (
        const int fn, 
        const std::vector<double> &gsb) 
{
    double err_est;
    int err_code;
    grfn f(gsb, fn);
    // infinite limit fails
    // but see https://www.r-bloggers.com/numerical-integration-over-an-infinite-interval-in-rcpp-2/
    // const double res = Numer::integrate(f, 0, R_PosInf, err_est, err_code);
    // ASSUME gsb[1] is scale
    const double res = Numer::integrate(f, 0, 20 * gsb[1], err_est, err_code);
    // Rprintf("res %10.7g err_est %10.7g err_code %4d \n", res, err_est, err_code);
    if (err_code>0) Rcpp::stop ("err_code>0 in hintegral1DNRcpp");
    return (res * 2);
}
//===============================================================

// integral of radial 2-D function 
double hintegral2DNRcpp (
    const int fn, 
    const std::vector<double> &gsb) 
{
  double err_est;
  int err_code;
  rgrfn f(gsb, fn);
  const double res = Numer::integrate(f, 0, 20*gsb[1], err_est, err_code);
  return (res * 2 * M_PI);
}
//===============================================================

