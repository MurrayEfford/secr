\name{read.mask}
\alias{read.mask}

\title{ Read Habitat Mask From File }

\description{
  Read coordinates of points on a habitat mask from a text file.
}

\usage{

read.mask(file = NULL, data = NULL, spacing = NULL, columns = NULL, ...)

}

\arguments{
  \item{file}{ character string with name of text file }
  \item{data}{ dataframe } 
  \item{spacing}{ spacing of grid points in metres }
  \item{columns}{ character vector naming the columns to save as covariates}
  \item{...}{ other arguments to pass to \code{read.table}}
}
\details{

For file input, the x and y coordinates are usually the first two values
on each line, separated by white space. If the file starts with a line
of column headers and `header = TRUE' is passed to \code{read.table} in
the \dots argument then `x' and `y' need not be the first two fields.

  \code{data} is an alternative input route if the x and y coordinates
  already exist in \R as columns in a dataframe. Only one of \code{data}
  or \code{file} should be specified.

The grid cell size \code{spacing} should be provided if known. If it is
not provided then an attempt is made to infer it from the minimum
spacing of points.  This can be slow and may demand more memory than is
available. In rare cases (highly fragmented masks) it may also yield the
wrong answer.

From 2.3.0, additional columns in the input are saved as covariates. The
default (\code{columns = NULL}) is to save all columns.

}

\value{
object of class \code{mask} with type `user'
}

\seealso{ \code{\link{mask}} }

\note{

  \code{read.mask} creates a single-session mask. If used in
\code{\link{secr.fit}} with a multi-session capthist object a
single-session mask will be replicated to the number of sessions. This
is appropriate if all sessions relate to the same geographical
region. If the `sessions' relate to different regions you will need to
construct a multi-session mask as a list of single-session
masks (e.g. \code{mask <- list(mask1, mask2, mask3)}).

}

\examples{
## Replace file name with a valid local name and remove `#'
# read.mask (file = "c:\\\\myfolder\\\\mask.txt",
# spacing = 3, header = TRUE)
## "mask.txt" should have lines like this
# x   y
# 265 265
# 268 265
# ...
}

\keyword{ IO }

