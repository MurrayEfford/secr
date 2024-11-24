\name{MCgof}
\alias{MCgof}
\alias{MCgof.secr}
\alias{MCgof.secrlist}

\title{Monte Carlo Goodness-of-fit for SECR Models}

\description{

\code{MCgof} implements and extends the Monte Carlo resampling method of Choo et al. (2024) to emulate Bayesian posterior predictive checks (Gelman et al. 1996, Royle et al. 2014). Initial results suggest the approach is more informative than the deviance-based test proposed by Borchers and Efford (2008) and implemented in \code{\link{secr.test}}. 
However, the tests have limited power.

\code{MCgof} is under development. The structure of the output may change and 
bugs may be found. See Warning below for exclusions.

}

\usage{

\method{MCgof}{secr}(object, nsim = 100, statfn = NULL, testfn = NULL, seed = NULL, 
    ncores = 1, clustertype = c("PSOCK", "FORK"), usefxi = TRUE, 
    useMVN = TRUE, Ndist = NULL, quiet = FALSE, debug = FALSE, ...)

\method{MCgof}{secrlist}(object, nsim = 100, statfn = NULL, testfn = NULL, seed = NULL, 
    ncores = 1, clustertype = c("PSOCK", "FORK"), usefxi = TRUE, 
    useMVN = TRUE, Ndist = NULL, quiet = FALSE, debug = FALSE, ...)

}

\arguments{

  \item{object}{secr fitted model or \code{\link{secrlist}} object}
  \item{nsim}{integer number of replicates}
  \item{statfn}{function to extract summary statistics from capture histories}
  \item{testfn}{function to compare observed and expected counts}
  \item{seed}{integer seed}
  \item{ncores}{integer for number of parallel cores}
  \item{clustertype}{character cluster type for parallel::makeCluster}
  \item{usefxi}{logical; if FALSE then AC are simulated de novo from the density 
  process rather than using information on the detected individuals}
  \item{useMVN}{logical; if FALSE parameter values are fixed at the MLE rather
  than drawn from multivariate normal distribution}
  \item{Ndist}{character; distribution of number of unobserved AC (optional) }
  \item{quiet}{logical; if FALSE then a progress bar (ncores=1) and final timing are shown}
  \item{debug}{integer; if >0 then the browser is started at one of 4 points in code}
  \item{\dots}{other arguments passed to testfn}
  
}

\details{

At each replicate parameter values are sampled from the multivariate-normal sampling distribution of the fitted model. The putative location of each detected individual is drawn from the spatial distribution implied by its observations and the resampled parameters (see \code{\link{fxi}}); locations of undetected individuals are simulated from the complement of pdot(x) times D(x). 

New detections are simulated under the model for individuals at the simulated locations, along with the expected numbers. Detections form a capthist object, a 3-D array with dimensions for individual \eqn{i}, occasion \eqn{j} and detector \eqn{k}*. Thus for each replicate and detected individual there are the original observations \eqn{y_{ijk}}, simulated observations \eqn{Y_{ijk}}, and expected counts \eqn{\mbox E (y_{ijk})}. Two discrepancy statistics are calculated for each replicate -- observed vs expected counts, and simulated vs expected counts -- and a record is kept of which of these discrepancy statistics is the larger (indicating poorer fit). 

* Notation differs slightly from Choo et al. (2024), using \eqn{j} for occasion and \eqn{k} for detector to be consistent with usage in \pkg{secr} and elsewhere (e.g., Borchers and Fewster 2016).

The default discrepancy (\code{testfn}) is the Freeman-Tukey statistic as in Choo et al. (2024) and Royle et al. (2014) (see also Brooks, Catchpole and Morgan 2000). The statistic has this general form for \eqn{M} counts \eqn{y_m} with expected value \eqn{\mbox E(y_m)}:
\deqn{T = \sum_{m=1}^{m=M} \left(\sqrt {y_m} - \sqrt{E(y_m)}\right)^2.}

The key output of \code{MCgof} is the proportion of replicates in which the simulated discrepancy exceeds the observed discrepancy. For perfect fit this will be about 0.5, and for poor fit it will approach zero.

By default, tests are performed separately for three types of count: the numbers of detections of each individual (yi), at each detector (yk), and for each individual at each detector (yik) extracted by the default \code{statfn} from the margins of the observed and simulated capture histories.

\tabular{llll}{
\eqn{y_{ik} = \sum_j y_{ijk}} \tab \tab individual x detector \cr 
\eqn{y_{i} = \sum_j \sum_k y_{ijk}} \tab\tab individual \cr
\eqn{y_{k} = \sum_i \sum_j y_{ijk}}  \tab\tab detector \cr 
}

Parallel processing is offered using multiple cores (CPUs) through the package \pkg{parallel} when ncores > 1. This differs from the usual multithreading paradigm in \pkg{secr} and does not rely on the environment variable set by \code{\link{setNumThreads}} \emph{except} that, if ncores = NULL, ncores will be set to the value from \code{setNumThreads}. The cluster type "FORK" is available only on Unix-like systems; it can require large amounts of memory, but is generally fast. A small value of ncores>1 may be optimal, especially with cluster type "PSOCK".

`usefxi' and `useMVN' may be used to drop key elements of the Choo et al. (2024) 
approach - they are provided for demonstration only.

`Ndist' refers to the distribution of the number of unobserved AC, conditional on the expected number \eqn{q = D^*A - n} where \eqn{D^*} is the resampled density, \eqn{A} the mask area, and \eqn{n} the number of detected individuals. By default `Ndist' depends on the distribution component of the `details' argument of the fitted model (``poisson" for Poisson \eqn{n}, ``fixed"" for binomial \eqn{n}). 

`debug' may be used to view intermediate data at certain points in MCgof() numbered 1 to 5. Examine the code of secr:::MCgof.secr or secr:::simfxiAC for these points. Debugging requires `ncores = 1'.

The \code{\link{RNGkind}} of the random number generator is set internally for consistency across platforms.

The \dots argument may be used to pass 'np' and 'verbose' to `Fletcher.chat` used as 'testfn'.
}

\value{

Invisibly returns an object of class `MCgof' with components -

  \item{nsim}{as input}
  \item{statfn}{as input or default}
  \item{testfn}{as input or default}
  \item{all}{list of outputs: for each statistic, a 3 x nsim matrix. Rows
   correspond to Tobs, Tsim, and a binary indicator for Tsim > Tobs}
  \item{proctime}{execution time in seconds}

For secrlist input the value returned is a list of `MCgof' objects.

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially
explicit maximum likelihood methods for capture--recapture studies.
\emph{Biometrics} \bold{64}, 377--385.

Borchers, D. L. and Fewster, R. M. (2016) Spatial capture--recapture models.  
\emph{Statistical Science} \bold{31}, 219--232.

Brooks, S. P., Catchpole, E. A. and Morgan, B. J. T. (2000) Bayesian animal 
survival estimation. \emph{Statistical Science} \bold{15}, 357--376.

Choo, Y. R., Sutherland, C. and Johnston, A. (2024) A Monte Carlo resampling 
framework for implementing goodness-of-fit tests in spatial capture-recapture model
\emph{Methods in Ecology and Evolution} DOI: 10.1111/2041-210X.14386.

Gelman, A., Meng, X.-L., and Stern, H. (1996) Posterior predictive assessment of 
model fitness via realized discrepancies. \emph{Statistica Sinica} \bold{6}, 733--807.

Royle, J. A., Chandler, R. B., Sollmann, R. and Gardner, B. (2014)
\emph{Spatial capture--recapture.} Academic Press.

}

\section{Warning}{

Not all models are covered and some are untested. These models are specifically 
excluded -

\enumerate{ 
  \item multi-session models
  \item models with groups
  \item conditional likelihood
  \item polygon, transect, telemetry or signal detectors
  \item non-binary behavioural responses
}

}

\section{Notes}{

This implementation extends the work of Choo et al. (2024) in these respects -

\enumerate{ 
  \item detector types `multi' and `count' are allowed
  \item the model may include variation among detectors
  \item the model may include behavioural responses
  \item 2-class finite mixture and hybrid mixture models are both allowed.
}

}

\author{

Murray Efford and Yan Ru Choo

}

\seealso{

\link{Parallel},
\code{\link{secr.test}}, 
\code{\link{plot.MCgof}}, 
\code{\link{hist.MCgof}}, 
\code{\link{summary.MCgof}}  

}

\examples{

\donttest{
tmp <- MCgof(secrdemo.0)
summary(tmp)
par(mfrow = c(1,3), pty = 's')
plot(tmp)
}

}
