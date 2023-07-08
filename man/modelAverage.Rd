\name{modelAverage}
\alias{modelAverage}
\alias{modelAverage.secr}
\alias{modelAverage.secrlist}
\alias{collate}
\alias{collate.secr}
\alias{collate.secrlist}

\title{ Averaging of SECR Models Using Akaike's Information Criterion}

\description{
  
AIC- or AICc-weighted average of estimated `real' or `beta' parameters
from multiple fitted secr models, and the tabulation of estimates.

}

\usage{

\method{modelAverage}{secr}(object, ..., realnames = NULL, betanames = NULL, newdata = NULL, 
    alpha = 0.05, dmax = 10, covar = FALSE, average = c("link", "real"),
    criterion = c("AICc","AIC"), CImethod = c("Wald", "MATA"), chat = NULL)

\method{modelAverage}{secrlist}(object, ..., realnames = NULL, betanames = NULL, newdata = NULL, 
    alpha = 0.05, dmax = 10, covar = FALSE, average = c("link", "real"),
    criterion = c("AICc","AIC"), CImethod = c("Wald", "MATA"), chat = NULL)

\method{collate}{secr}(object, ..., realnames = NULL, betanames = NULL, newdata = NULL, 
    alpha = 0.05, perm = 1:4, fields = 1:4)

\method{collate}{secrlist}(object, ..., realnames = NULL, betanames = NULL, newdata = NULL, 
    alpha = 0.05, perm = 1:4, fields = 1:4)

}

\arguments{
  \item{object}{secr or secrlist object}
  \item{\dots}{other secr objects } 
  \item{realnames}{ character vector of real parameter names }
  \item{betanames}{ character vector of beta parameter names }
  \item{newdata}{ optional dataframe of values at which to evaluate models }
  \item{alpha}{ alpha level for confidence intervals }
  \item{dmax}{ numeric, the maximum AIC or AICc difference for inclusion in confidence set }
  \item{covar}{ logical, if TRUE then return variance-covariance matrix }
  \item{average}{ character string for scale on which to average real
    parameters }
  \item{criterion}{ character, information criterion to use for model weights}  
  \item{CImethod}{ character, type of confidence interval (see Details)}
  \item{chat}{numeric optional variance inflation factor for quasi-AIC weights}
  \item{perm}{ permutation of dimensions in output from \code{collate} }
  \item{fields}{ vector to restrict summary fields in output }
}

\details{

Models to be compared must have been fitted to the same data and use the
same likelihood method (full vs conditional). If \code{realnames} =
NULL and \code{betanames} = NULL then all real parameters will be
averaged; in this case all models must use the same real parameters. To
average beta parameters, specify \code{betanames} (this is ignored if a
value is provided for \code{realnames}). See \code{\link{predict.secr}}
for an explanation of the optional argument \code{newdata};
\code{newdata} is ignored when averaging beta parameters. 

Model-averaged estimates for parameter \eqn{\theta}{theta} are given by
\deqn{\hat{\theta} = \sum\limits _k w_k \hat{\theta}_k}{theta-hat = sum(
w_k * theta-hat_k)} where the subscript \eqn{k} refers to a specific
model and the \eqn{w_k} are AIC or AICc weights (see
\code{\link{AIC.secr}} for details). Averaging of real parameters may be
done on the link scale before back-transformation
(\code{average="link"}) or after back-transformation
(\code{average="real"}).

Models for which dAIC > \code{dmax} (or dAICc > \code{dmax}) are given a
weight of zero and effectively are excluded from averaging.

Also, \deqn{\mbox{var} (\hat{\theta}) = \sum\limits _{k} { w_{k}
( \mbox{var}(\hat{\theta}_{k} | \beta _k) + \beta _k ^2)}
}{var(theta-hat) = sum(w_k (var(theta-hat_k) + beta_k^2))}

where \eqn{\hat{\beta} _k = \hat{\theta}_k - \hat{\theta}}{beta-hat_k =
theta-hat_k -- theta-hat} and the variances are asymptotic estimates
from fitting each model \eqn{k}. This follows Burnham and Anderson
(2004) rather than Buckland et al. (1997).

Two methods are offered for confidence intervals. The default `Wald'
uses the above estimate of variance. The alternative `MATA'
(model-averaged tail area) avoids estimating a weighted variance and
is thought to provide better coverage at little cost in increased
interval length (Turek and Fletcher 2012). Turek and Fletcher (2012)
also found averaging with AIC weights (here \code{criterion = 'AIC'})
preferable to using AICc weights, even for small
samples. \code{CImethod} does not affect the reported standard errors.

If 'chat' is provided then quasi-AIC or quasi-AICc weights are used, 
depending on the value of 'criterion'.

\code{collate} extracts parameter estimates from a set of fitted secr
model objects. \code{fields} may be used to select a subset of summary
fields ("estimate","SE.estimate","lcl","ucl") by name or number. 

}

\value{

For \code{modelAverage}, an array of model-averaged estimates, their
standard errors, and a \eqn{100(1-\alpha)}{100(1-alpha)}\% confidence
interval. The interval for real parameters is backtransformed from the
link scale. If there is only one row in \code{newdata} or beta
parameters are averaged or averaging is requested for only one parameter
then the array is collapsed to a matrix. If \code{covar = TRUE} then a
list is returned with separate components for the estimates and the
variance-covariance matrices.

For \code{collate}, a 4-dimensional array of model-specific parameter
estimates. By default, the dimensions correspond respectively to rows in
\code{newdata} (usually sessions), models, statistic fields (estimate, SE.estimate, lcl,
ucl), and parameters ("D", "g0" etc.). For particular comparisons it often helps
to reorder the dimensions with the \code{perm} argument.

}
\references{ 

Buckland S. T., Burnham K. P. and Augustin, N. H. (1997) Model
selection: an integral part of inference. \emph{Biometrics} \bold{53},
603--618.

Burnham, K. P. and Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: A Practical Information-Theoretic Approach}.
Second edition. New York: Springer-Verlag.

Burnham, K. P. and Anderson, D. R. (2004) Multimodel inference -
understanding AIC and BIC in model selection. \emph{Sociological Methods
& Research} \bold{33}, 261--304.

Turek, D. and Fletcher, D. (2012) Model-averaged Wald confidence
intervals. \emph{Computational statistics and data analysis} \bold{56},
2809--2815.

}

\note{\code{modelAverage} replaces the deprecated function \code{model.average} 
whose name conflicted with a method in \pkg{RMark}.}

\seealso{ \code{\link{modelAverage}}, \code{\link{AIC.secr}}, \code{\link{secr.fit}}}

\examples{
## Compare two models fitted previously
## secrdemo.0 is a null model
## secrdemo.b has a learned trap response

modelAverage(secrdemo.0, secrdemo.b)
modelAverage(secrdemo.0, secrdemo.b, betanames = c("D","g0","sigma"))

## In this case we find the difference was actually trivial...
## (subscripting of output is equivalent to setting fields = 1)

collate (secrdemo.0, secrdemo.b, perm = c(4,2,3,1))[,,1,]

}

\keyword{ models }

