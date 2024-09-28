\name{closure.test}
\alias{closure.test}
\title{
Closure tests
}
\description{
  
Perform tests to determine whether a population sampled by
capture-recapture is closed to gains and losses over the period of sampling.

}

\usage{

closure.test(object, SB = FALSE, min.expected = 2)

}

\arguments{

  \item{object}{\code{capthist} object}
  \item{SB}{logical, if TRUE then test of Stanley and Burnham 1999 is
    calculated in addition to that of Otis et al. 1978}
  \item{min.expected}{integer for the minimum expected count in any cell of a component 2x2 table}

}
\details{

  The test of Stanley and Burnham in part uses a sum over 2x2 contingency
tables; any table with a cell whose expected count is less than
min.expected is dropped from the sum. The default value of 2 is that
used by CloseTest (Stanley and Richards 2005, T. Stanley pers. comm.;
see also Stanley and Burnham 1999 p. 203).

}

\value{

In the case of a single-session capthist object, either a vector with
the statistic (z-value) and p-value for the test of Otis et al. (1978
p. 120) or a list whose components are data frames with the statistics
and p-values for various tests and test components as follows --
\item{Otis }{Test of Otis et al. 1978}
\item{Xc }{Overall test of Stanley and Burnham 1999}
\item{NRvsJS }{Stanley and Burnham 1999}
\item{NMvsJS }{Stanley and Burnham 1999}
\item{MtvsNR }{Stanley and Burnham 1999}
\item{MtvsNM }{Stanley and Burnham 1999}
\item{compNRvsJS }{Occasion-specific components of NRvsJS}
\item{compNMvsJS }{Occasion-specific components of NMvsJS}

Check the original papers for an explanation of the components of the Stanley and Burnham test.

In the case of a multi-session object, a list with one component (as
above) for each session.

}

\references{

Otis, D. L., Burnham, K. P., White, G. C. and Anderson, D. R. (1978)
  Statistical inference from capture data on closed animal
  populations. \emph{Wildlife Monographs} \bold{62}, 1--135.

Stanley, T. R. and Burnham, K. P. (1999) A closure test for
  time-specific capture--recapture data. \emph{Environmental and
  Ecological Statistics} \bold{6}, 197--209.

Stanley, T. R. and Richards, J. D. (2005) A program for testing
  capture--recapture data for closure. \emph{Wildlife Society Bulletin}
  \bold{33}, 782--785.

White, G. C., Anderson, D. R., Burnham, K. P. and Otis, D. L. (1982)
\emph{Capture-recapture and removal methods for sampling closed
populations}. Los Alamos National Laboratory, Los Alamos, New Mexico.

}

\note{
  
  No omnibus test exists for closure: the existing tests may indicate
nonclosure even when a population is closed if other effects such as
trap response are present (see White et al. 1982 pp 96--97). The test of
Stanley and Burnham is sensitive to individual heterogeneity which is
inevitable in most spatial sampling, and it should not in general be
used for this sort of data.

}

\seealso{
  \code{\link{capthist}}
}
\examples{
 closure.test(captdata)
}
