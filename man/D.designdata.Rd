\name{D.designdata}
\alias{D.designdata}

\title{ Construct Density Design Data }

\description{
  
  Internal function used by \code{\link{secr.fit}},
  \code{\link{confint.secr}}, and \code{\link{score.test}}. 
  
}

\usage{
D.designdata (mask, Dmodel, grouplevels, sessionlevels, sessioncov =
NULL, meanSD = NULL)
}

\arguments{
  \item{mask}{ \code{\link{mask}} object.}
  \item{Dmodel}{ formula for density model }
  \item{grouplevels}{ vector of group names }
  \item{sessionlevels}{ vector of character values for session names }
  \item{sessioncov}{ optional dataframe of values of session-specific
    covariate(s). }
  \item{meanSD}{optional external values for scaling x- and
  y- coordinates}
}
\details{

This is an internal \pkg{secr} function that you are unlikely ever to
use. Unlike \code{\link{secr.design.MS}}, this function does \emph{not}
call \code{model.matrix}.

}

\value{

Dataframe with one row for each combination of mask point, group and
session. Conceptually, we use a 3-D rectangular array with enough rows
to accommodate the largest mask, so some rows in the output may merely
hold space to enable easy indexing. The dataframe has an attribute
`dimD' that gives the relevant dimensions: \code{ attr(dframe, "dimD") =
c(nmask, ngrp, R)}, where \code{nmask} is the number of mask points,
\code{ngrp} is the number of groups, and \code{R} is the number of
sessions. Columns correspond to predictor variables in Dmodel.

The number of valid rows (points in each session-specific mask) is
stored in the attribute `validMaskRows'.

For a single-session mask, \code{meanSD} is a 2 x 2 matrix of mean and
SD (rows) for x- and y-coordinates. For a multi-session mask, a list of
such objects. Ordinarily these values are from the meanSD attribute of
the mask, but they must be specified when applying a new mask to an
existing model.

}

\seealso{ \code{\link{secr.design.MS}} }

