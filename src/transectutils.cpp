#include "poly.h"

using namespace Rcpp;

// exports utility functions ontransectcpp and alongtransectcpp

rpoint getxycpp(
        const double l, 
        const std::vector<double> &cumd, 
        const RcppParallel::RMatrix<double> &line, 
        const int n1, 
        const int n2) {
    // return the xy coordinates of point l metres along a transect 
    // n1 is the starting position for this transect within 'line'
    int j = 1;  // initialised 2022-01-18
    double pr, d, d12;
    rpoint xy;
    // refined by PJ GitHub PR 2023-05-15
    auto upper = std::upper_bound(cumd.begin() + 1,
                                  cumd.begin() + (n2 - n1), l);
    j = std::distance(cumd.begin(), upper);
    d = l - cumd[j-1];  // distance along leg 
    d12 = cumd[j] - cumd[j-1];
    if (d12>0)
        pr = d / d12;
    else
        pr = 0;
    
    j = j+n1;
    xy.x = line(j-1,0) + (line(j,0) - line(j-1,0)) * pr;
    xy.y = line(j-1,1) + (line(j,1) - line(j-1,1)) * pr;
    return(xy);
}

// [[Rcpp::export]]
bool ontransectcpp (
    NumericVector xy,
    NumericMatrix transect,
    int    n1,
    int    n2,
    double tol)
{
    // Is point xy on transect?
        // We assume transect coordinates are in col-major order (x's then y's)
    // http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/ 2009-11-29
    // http://paulbourke.net/geometry/pointlineplane/ 2019-10-13
    int k;
    double r;
    double u;
    rpoint p,p1,p2,p3;
    double minr = 1e20;
    if (n2>=transect.nrow()) Rcpp::stop ("invalid input ontransectcpp");
    
    p3.x = xy(0);
    p3.y = xy(1);
    
    for (k= n1; k < n2; k++)
    {
        p1.x = transect(k,0);
        p1.y = transect(k,1);
        p2.x = transect(k+1,0);
        p2.y = transect(k+1,1);
        if (distance1(p1,p2) > 0) {
            u = ((p3.x-p1.x) * (p2.x-p1.x) + (p3.y-p1.y) * (p2.y-p1.y)) /
                ((p2.x-p1.x) * (p2.x-p1.x) + (p2.y-p1.y) * (p2.y-p1.y));
            if ((u>=0) && (u<=1)) {
                p.x = p1.x + u * (p2.x-p1.x);
                p.y = p1.y + u * (p2.y-p1.y);
                r = distance1 (p,p3);
                minr = std::min(r,minr);
            }
        }
    }
    // check at each vertex to be sure 
    for (k= n1; k <= n2; k++) {
        p1.x = transect(k,0);
        p1.y = transect(k,1);
        r = distance1 (p1,p3);
        minr = std::min(r,minr);
    }
    return(minr < tol);
}
//----------------------------------------------------------------
    
    // [[Rcpp::export]]
double alongtransectcpp (
    NumericVector xy,
    NumericMatrix transect,
    int    n1,
    int    n2,
    double tol
)
{
    // How far is point xy from start of transect? 2011-06-07
    // We assume transect coordinates are in col-major order (x's then y's)
    // http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/ 29/11/09
    int k;
    double r;
    double u;
    rpoint p,p1,p2,p3;
    
    double along=0.0;
    if (n2>=transect.nrow()) Rcpp::stop ("invalid input alongtransectcpp");
    p3.x = xy(0);
    p3.y = xy(1);
    
    for (k= n1; k < n2; k++)
    {
        p1.x = transect(k,0);
        p1.y = transect(k,1);
        r = distance1 (p1,p3);
        if (r < tol) {
            return(along);
        } 
        p2.x = transect(k+1,0);
        p2.y = transect(k+1,1);
        if (distance1(p1,p2) > 0) {
            u = ((p3.x-p1.x) * (p2.x-p1.x) + (p3.y-p1.y) * (p2.y-p1.y)) /
                ((p2.x-p1.x) * (p2.x-p1.x) + (p2.y-p1.y) * (p2.y-p1.y));
            if ((u>=0) && (u<=1)) {
                p.x = p1.x + u * (p2.x-p1.x);
                p.y = p1.y + u * (p2.y-p1.y);
                r = distance1 (p,p3);
                if (r < tol) {
                    along += distance1(p,p1);
                    return(along);
                } 
            }
            along += distance1(p1,p2);
        }
    }
    return(along);
}
//----------------------------------------------------------------
    