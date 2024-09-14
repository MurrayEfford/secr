\name{closedN}
\alias{closedN}

\title{
Closed population estimates
}
\description{

  Estimate N, the size of a closed population, by several conventional
  non-spatial capture--recapture methods. 

}

\usage{

closedN(object, estimator = NULL, level = 0.95, maxN = 1e+07,
    dmax = 10 )

}

\arguments{
  \item{object}{ \code{\link{capthist}} object }
  \item{estimator}{ character; name of estimator (see Details) }
  \item{level}{ confidence level (1 -- alpha) }
  \item{maxN}{ upper bound for population size }
  \item{dmax}{ numeric, the maximum AIC difference for inclusion in confidence set}  
}

\details{

  Data are provided as spatial capture histories, but the spatial
  information (trapping locations) is ignored.

  AIC-based model selection is available for the maximum-likelihood
  estimators \code{null}, \code{zippin}, \code{darroch}, \code{h2}, and
  \code{betabinomial}.

  Model weights are calculated as \deqn{w_i = \frac{\exp(-\Delta_i / 2)}{
  \sum{\exp(-\Delta_i / 2)}}}{w_i = exp(-dAICc_i / 2) / sum{ exp(-dAICc_i
  / 2) }}

  Models for which dAICc > \code{dmax} are given a weight of zero and are
  excluded from the summation, as are non-likelihood models.  

  Computation of \code{null}, \code{zippin} and \code{darroch} estimates
  differs slightly from Otis et al. (1978) in that the likelihood is
  maximized over real values of N between \code{Mt1} and \code{maxN},
  whereas Otis et al. considered only integer values.

  Asymmetric confidence intervals are obtained in the same way for all
  estimators, using a log transformation of \eqn{\hat{N}-Mt1}{Nhat-Mt1}
  following Burnham et al. (1987), Chao (1987) and Rexstad and Burnham
  (1991).

  The available estimators are
  \tabular{llll}{
  Name \tab Model \tab Description \tab Reference \cr
  \code{null} \tab M0 \tab null \tab Otis et al. 1978 p.105 \cr
  \code{zippin} \tab Mb \tab removal \tab Otis et al. 1978 p.108 \cr
  \code{darroch} \tab Mt \tab Darroch \tab Otis et al. 1978 p.106-7 \cr
  \code{h2} \tab Mh \tab 2-part finite mixture \tab Pledger 2000 \cr
  \code{betabinomial} \tab Mh \tab Beta-binomial continuous mixture \tab
  Dorazio and Royle 2003 \cr
  \code{jackknife} \tab Mh \tab jackknife \tab  Burnham and Overton 1978 \cr
  \code{chao} \tab Mh \tab Chao's Mh estimator \tab Chao 1987 \cr
  \code{chaomod} \tab Mh \tab Chao's modified Mh estimator \tab Chao 1987 \cr
  \code{chao.th1} \tab Mth \tab sample coverage estimator 1 \tab Lee and Chao 1994 \cr
  \code{chao.th2} \tab Mth \tab sample coverage estimator 2 \tab Lee and Chao 1994 \cr
  }

}
\value{
  A dataframe with one row per estimator and columns
  \item{model }{ model in the sense of Otis et al. 1978}
  \item{npar }{ number of parameters estimated}
  \item{loglik }{ maximized log likelihood}
  \item{AIC }{ Akaike's information criterion}
  \item{AICc }{AIC with small-sample adjustment of Hurvich & Tsai (1989)}
  \item{dAICc }{ difference between AICc of this model and the one with
    smallest AICc}
  \item{Mt1 }{ number of distinct individuals caught}
  \item{Nhat }{ estimate of population size}
  \item{seNhat }{ estimated standard error of Nhat}
  \item{lclNhat }{ lower 100 x level \% confidence limit }
  \item{uclNhat }{ upper 100 x level \% confidence limit }

}

\note{
  
  Prof. Anne Chao generously allowed me to adapt her code for the
  variance of the `chao.th1' and `chao.th2' estimators.
  
  Chao's estimators have been subject to various improvements not
  included here (e.g., Chao et al. 2016).

}

\section{Warning}{
  
  If your data are from spatial sampling (e.g. grid trapping) it is
  recommended that you do \emph{not} use these methods to estimate
  population size (see Efford and Fewster 2013). Instead, fit a spatial model
  and estimate population size with \code{\link{region.N}}.
  
}

\references{

Burnham, K. P. and Overton, W. S. (1978) Estimating the size of a closed
population when capture probabilities vary among
animals. \emph{Biometrika} \bold{65}, 625--633.

Chao, A. (1987) Estimating the population size for capture--recapture
data with unequal catchability. \emph{Biometrics} \bold{43}, 783--791.

Chao, A., Ma, K. H., Hsieh,  T. C. and Chiu, Chun-Huo (2016) SpadeR: Species-Richness Prediction and
Diversity Estimation with R. R package version 0.1.1. https://CRAN.R-project.org/package=SpadeR

Dorazio, R. M. and Royle, J. A. (2003) Mixture models for estimating the
size of a closed population when capture rates vary among
individuals. \emph{Biometrics} \bold{59}, 351--364.

Efford, M. G. and Fewster, R. M. (2013) Estimating population
  size by spatially explicit capture--recapture. \emph{Oikos}
  \bold{122}, 918--928.

Hurvich, C. M. and Tsai, C. L. (1989) Regression and time series model
selection in small samples. \emph{Biometrika} \bold{76}, 297--307.

Lee, S.-M. and Chao, A. (1994) Estimating population size via sample
coverage for closed capture-recapture models. \emph{Biometrics}
\bold{50}, 88--97.

Otis, D. L., Burnham, K. P., White, G. C. and Anderson, D. R. (1978)
  Statistical inference from capture data on closed animal
  populations. \emph{Wildlife Monographs} \bold{62}, 1--135.

Pledger, S. (2000) Unified maximum likelihood estimates for closed
capture-recapture models using mixtures. \emph{Biometrics} \bold{56},
434--442.

Rexstad, E. and Burnham, K. (1991) User's guide for interactive program
CAPTURE. Colorado Cooperative Fish and Wildlife Research Unit, Fort
Collins, Colorado, USA.

}

\seealso{
  
  \code{\link{capthist}},
  \code{\link{closure.test}},
  \code{\link{region.N}}

}

\examples{
closedN(deermouse.ESG)
}

\keyword{ models }