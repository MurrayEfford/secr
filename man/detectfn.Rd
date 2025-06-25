\name{detectfn}
\alias{Detection functions}
\alias{detectfn}

\title{ Detection Functions }

\description{

A detection function relates the probability of detection \eqn{g} or the
expected number of detections \eqn{\lambda} for an animal to the
distance of a detector from a point usually thought of as its home-range
centre. In \pkg{secr} only simple 2- or 3-parameter functions are
used. Each type of function is identified by its number or by a 2--3
letter code (version \eqn{\ge} 2.6.0; see below). In most cases the name
may also be used (as a quoted string).

Choice of detection function is usually not critical, and the default
`HN' is usually adequate.

Functions (14)--(19) are parameterised in terms of the expected number
of detections \eqn{\lambda}, or cumulative hazard, rather than
probability. `Exposure' (e.g. Royle and Gardner 2011) is another term
for cumulative hazard. This parameterisation is natural for the `count'
\link{detector} type or if the function is to be interpreted as a
distribution of activity (home range). When one of the functions
(14)--(19) is used to describe detection probability (i.e., for the binary
detectors `single', `multi',`proximity',`polygonX' or
`transectX'), the expected number of detections is internally
transformed to a binomial probability using \eqn{g(d) =
1-\exp(-\lambda(d))}{g(d) = 1-exp(-\lambda(d))}.

The hazard halfnormal (14) is similar to the halfnormal exposure function
used by Royle and Gardner (2011) except they omit the factor of 2 on
\eqn{\sigma^2}, which leads to estimates of \eqn{\sigma} that are larger
by a factor of sqrt(2). The hazard exponential (16) is identical to their
exponential function.

\tabular{llll}{
  Code \tab Name \tab Parameters \tab Function \cr
  
0  HN \tab halfnormal \tab g0, sigma \tab \eqn{ g(d) = g_0 \exp
  \left(\frac{-d^2} {2\sigma^2} \right) }{ g(d) = g0 exp{-d^2 / (2
  \sigma^2) } } \cr

1  HR \tab hazard rate \tab g0, sigma, z \tab \eqn{ g(d) = g_0 [1 - \exp\{
  {-(^d/_\sigma)^{-z}} \}] }{ g(d) = g0 (1 - exp(- (d / \sigma)
  ^(-z) )) } \cr

2  EX \tab exponential  \tab g0, sigma \tab \eqn{ g(d) = g_0  \exp \{
    -(^d/_\sigma) \} }{ g(d) = g0 exp(- d / \sigma) } \cr
  
3 CHN \tab compound halfnormal  \tab g0, sigma, z \tab \eqn{ g(d) = g_0 [1
  - \{1 - \exp \left(\frac{-d^2} {2\sigma^2} \right)\} ^ z] }{ g(d) = g0
   [1 - \{1 - exp(-d^2 / (2 \sigma^2))]^z\} } \cr

4  UN \tab uniform  \tab g0, sigma \tab \eqn{ g(d) = g_0, d <= \sigma;
  g(d) = 0, \mathrm{otherwise}  }{ g(d) = g0, d\le \sigma; g(d) = 0,
  otherwise } \cr

5 WEX \tab w exponential  \tab g0, sigma, w \tab \eqn{ g(d) = g_0, d < w;
  g(d) = g_0  \exp \left( -\frac{d-w}{\sigma} \right), \mathrm{otherwise}
  }{ g(d) = g0, d < w; g(d) = g0 exp(- (d - w) / \sigma), otherwise }
 \cr
 
6 ANN \tab annular normal \tab g0, sigma, w \tab \eqn{ g(d) = g_0 \exp
  \lbrace \frac{-(d-w)^2} {2\sigma^2} \rbrace }{ g(d) = g0 exp(-(d-w)^2 / (2 \sigma^2)) } \cr

7 CLN \tab cumulative lognormal \tab g0, sigma, z \tab \eqn{ g(d) = g_0
  [ 1 - F \lbrace(d-\mu)/s \rbrace ] }{ g(d) = g0 [1 - F\{(d-\mu)/s)\}] }
 \cr

8  CG \tab cumulative gamma \tab g0, sigma, z \tab \eqn{ g(d) = g_0
  \lbrace 1 - G (d; k, \theta)\rbrace }{ g(d) = g0 \{ 1 - G (d; k, \theta) \} }
 \cr

9 BSS \tab binary signal strength  \tab b0, b1 \tab \eqn{ g(d) = 1 - F
    \lbrace - ( b_0 + b_1  d) \rbrace }{ g(d) = 1 - F \{- (b0 + b1 . d) \} } \cr
  
10 SS \tab signal strength  \tab beta0, beta1, sdS \tab \eqn{ g(d) =1 -
  F[\lbrace c - (\beta_0 + \beta_1  d) \rbrace / s] }{ g(d) = 1 - F[ \{c - (\beta0 + \beta1
  . d)\} / sdS] } \cr

11 SSS \tab signal strength spherical \tab beta0, beta1, sdS \tab
    \eqn{ g(d) = 1 - F [ \lbrace c - (\beta_0 + \beta_1 (d-1) -  10 \log
      _{10}  d^2 ) \rbrace / s ]}{ g(d) = 1 - F[\{c - (\beta0 + \beta1 (d-1) - 10 log10
   ( d^2 ) ) \} / sdS ] } \cr

14 HHN \tab hazard halfnormal \tab lambda0, sigma \tab \eqn{ \lambda(d) = \lambda_0 \exp
  \left(\frac{-d^2} {2\sigma^2} \right) }{ \lambda(d) = \lambda0 exp{-d^2 / (2
  \sigma^2) } }; \eqn{g(d) = 1-\exp(-\lambda(d))}{g(d) = 1-exp(-\lambda(d))}\cr

15 HHR \tab hazard hazard rate \tab lambda0, sigma, z \tab \eqn{ \lambda(d)
  = \lambda_0 (1 - \exp \{ -(^d/_\sigma)^{-z} \}) }{ \lambda(d) =
  \lambda0 (1 - exp(- (d / \sigma)^(-z))) }; \eqn{g(d) = 1-\exp(-\lambda(d))}{g(d) = 1-exp(-\lambda(d))} \cr
 
16 HEX \tab hazard exponential \tab lambda0, sigma \tab \eqn{ \lambda(d)
  = \lambda_0 \exp \{ -(^d/_\sigma) \} }{ \lambda(d) = \lambda0 exp(- d
  / \sigma) }; \eqn{g(d) = 1-\exp(-\lambda(d))}{g(d) = 1-exp(-\lambda(d))} \cr
 
17 HAN \tab hazard annular normal \tab lambda0, sigma, w \tab \eqn{ \lambda(d) = \lambda_0 \exp
  \lbrace \frac{-(d-w)^2} {2\sigma^2} \rbrace }{ \lambda(d) = \lambda0 exp(-(d-w)^2 / (2 \sigma^2)) }; \eqn{g(d) = 1-\exp(-\lambda(d))}{g(d) = 1-exp(-\lambda(d))} \cr

18 HCG \tab hazard cumulative gamma \tab lambda0, sigma, z \tab \eqn{ \lambda(d) = \lambda_0
  \lbrace 1 - G (d; k, \theta)\rbrace }{ \lambda(d) = \lambda0 \{ 1 - G
  (d; k, \theta) \} }; \eqn{g(d) = 1-\exp(-\lambda(d))}{g(d) = 1-exp(-\lambda(d))}
 \cr

19 HVP \tab hazard variable power \tab lambda0, sigma, z \tab \eqn{ \lambda(d)
  = \lambda_0 \exp \{ -(^d/_\sigma)^{z} \} }{ \lambda(d) =
  \lambda0 exp(- (d / \sigma)^z) }; \eqn{g(d) = 1-\exp(-\lambda(d))}{g(d) = 1-exp(-\lambda(d))} \cr
  
20 OU \tab Ornstein-Uhlenbeck movement \tab epsilon, sigma, tau \tab See \code{\link{simOU}} \cr
 
}

Functions (1) and (15), the "hazard-rate" detection functions described by Hayes and Buckland
(1983), are not recommended for SECR because of their long tail, and
care is also needed with (2) and (16).

Function (3), the compound halfnormal, follows Efford and Dawson (2009).

Function (4) uniform is defined only for simulation as it poses problems
for likelihood maximisation by gradient methods. Uniform probability
implies uniform hazard, so there is no separate function `HUN'.

For function (7), `F' is the standard normal distribution function and
\eqn{\mu}{\mu} and \eqn{s}{s} are the mean and standard deviation on the
log scale of a latent variable representing a threshold of detection
distance. See Note for the relationship to the fitted parameters sigma
and z.

For functions (8) and (18), `G' is the cumulative distribution function of the
gamma distribution with shape parameter \emph{k} ( = \code{z}) and scale
parameter \eqn{\theta}{\theta} ( = \code{sigma/z}). See R's
\code{\link{pgamma}}.

For functions (9), (10) and (11), `F' is the standard normal
distribution function and \eqn{c} is an arbitrary signal threshold. The two
parameters of (9) are functions of the parameters of (10) and (11):
\eqn{b_0 = (\beta_0 - c) / sdS}{b0 = (\beta0-c)/sdS} and \eqn{b_1 =
\beta_1 / s}{b1 = \beta1/s} (see Efford et al. 2009). Note that (9) does
not require signal-strength data or \eqn{c}.

Function (11) includes an additional `hard-wired' term for sound
attenuation due to spherical spreading. Detection probability at
distances less than 1 m is given by \eqn{g(d) = 1 - F \lbrace(c -
  \beta_0) / sdS \rbrace}{ g(d) = 1 - F\{(c - beta0) / sdS \} }

Functions (12) and (13) are undocumented methods for sound attenuation.

Function (19) has been used in some published papers and is included for 
comparison (see e.g. Ergon and Gardner 2014).

Function (20) defines a detection model based on a bivariate Ornstein-Uhlenbeck 
movement model (mean-reverting random walk). Like the uniform (4), it is used 
in \pkg{secr} only for simulation. The three parameters are epsilon (radius of a 
detector within which detection occurs), sigma (asymptotic bivariate normal 
scale, as in HHN) and tau (serial correlation parameter \eqn{\tau = 1/\beta}). 
See \code{\link{simOU}}.

}

\references{

Efford, M. G. and Dawson, D. K. (2009) Effect of distance-related
heterogeneity on population size estimates from point counts. \emph{Auk}
\bold{126}, 100--111.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

Ergon, T. and Gardner, B. (2014) Separating mortality and emigration: 
modelling space use, dispersal and survival with robust-design spatial 
capture--recapture data. \emph{Methods in Ecology and Evolution} \bold{5}, 
1327--1336.

Hayes, R. J. and Buckland, S. T. (1983) Radial-distance models for the
line-transect method. \emph{Biometrics} \bold{39}, 29--42.

Royle, J. A. and Gardner, B. (2011) Hierarchical spatial
capture--recapture models for estimating density from trapping
arrays. In: A.F. O'Connell, J.D. Nichols & K.U. Karanth (eds)
\emph{Camera Traps in Animal Ecology: Methods and Analyses}. Springer,
Tokyo.  Pp. 163--190.

}

\seealso{ \code{\link{detectfnplot}} }

\note{
    
  The parameters of function (7) are potentially confusing.  The fitted
parameters describe a latent threshold variable on the natural scale:
sigma (mean) = \eqn{\exp(\mu + s^2 / 2)}{exp(\mu + s^2/2)} and z
(standard deviation) = \eqn{\sqrt{\exp(s^2 + 2
\mu)(\exp(s^2)-1)}}{sqrt(exp(s^2+2\mu)(exp(s^2)-1)))}. As with other
detection functions, sigma is a spatial scale parameter, although in
this case it corresponds to the mean of the threshold variable; the
standard deviation of the threshold variable (z) determines the shape
(roughly 1/max(slope)) of the detection function.  }


