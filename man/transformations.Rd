\name{transformations}
\alias{transformations}
\alias{flip}
\alias{flip.traps}
\alias{flip.popn}
\alias{rotate}
\alias{rotate.traps}
\alias{rotate.popn}
\alias{shift}
\alias{shift.traps}
\alias{shift.popn}
\alias{shift.mask}

\title{ Transform Point Array }

\description{
  
Flip (reflect), rotate or slide (translate) an array of points. Methods
are provided for `traps' objects that ensure other attributes are
retained. The methods may be used with \code{\link{rbind.traps}} to
create complex geometries.

}

\usage{

flip (object, lr = FALSE, tb = FALSE, ...)
rotate (object, degrees, centrexy = NULL, \dots)
shift (object, shiftxy, ...)

\method{flip}{traps} (object, lr = FALSE, tb = FALSE, ...)
\method{rotate}{traps}(object, degrees, centrexy = NULL, \dots)
\method{shift}{traps}(object, shiftxy, \dots)

\method{flip}{popn} (object, lr = FALSE, tb = FALSE, ...)
\method{rotate}{popn}(object, degrees, centrexy = NULL, \dots)
\method{shift}{popn}(object, shiftxy, \dots)

\method{shift}{mask}(object, shiftxy, \dots)

}

\arguments{
  \item{object}{ a 2-column matrix or object that can be coerced to a
    matrix }
  
  \item{lr}{ either logical for whether array should be flipped
    left-right, or numeric value for x-coordinate of axis about which it
    should be flipped left-right }

  \item{tb}{ either logical for whether array should be flipped
    top-bottom, or numeric value for y-coordinate of axis about which it
    should be flipped top-bottom }

  \item{degrees}{ clockwise angle of rotation in degrees }

  \item{centrexy}{ vector with xy coordinates of rotation centre }
  
  \item{shiftxy}{ vector of x and y displacements }
  
  \item{\dots}{ other arguments (not used) }
}

\details{

\code{flip} reflects points about a vertical or horizontal axis. Logical
values for \code{lr} or \code{tb} indicate that points should be flipped
about the mean on the relevant axis. Numeric values indicate the
particular axis value(s) about which points should be flipped.  The
default arguments result in no change.

\code{shift} shifts the location of each point by the desired amount
on each axis.

\code{rotate} rotates the array about a designated point. If
\code{centrexy} is NULL then rotation is about (0,0)
(\code{rotate.default}), the array centre (\code{rotate.traps}), or the
centre of the bounding box (\code{rotate.popn}).

}

\value{
  
A matrix or object of class `traps' or `popn' with the coordinates of each point
transformed as requested.

}

\seealso{ \code{\link{traps}}, \code{\link{popn}} }

\examples{

temp <- matrix(runif (20) * 2 - 1, nc = 2)

## flip
temp2 <- flip(temp, lr = 1)
plot(temp, xlim=c(-1.5,4), ylim = c(-1.5,1.5), pch = 16)
points (temp2, pch = 1)
arrows (temp[,1], temp[,2], temp2[,1], temp2[,2], length = 0.1)
abline(v = 1, lty = 2)

## rotate
temp2 <- rotate(temp, 25)
plot(temp, xlim=c(-1.5,1.5), ylim = c(-1.5,1.5), pch = 16)
points (0,0, pch=2)
points (temp2, pch = 1)
arrows (temp[,1], temp[,2], temp2[,1], temp2[,2], length = 0.1)	

## shiftxy
temp2 <- shift(temp, c(0.1, 0.1))
plot(temp, xlim=c(-1.5,1.5), ylim = c(-1.5,1.5), pch = 16)
points (0,0, pch=2)
points (temp2, pch = 1)
arrows (temp[,1], temp[,2], temp2[,1], temp2[,2], length = 0.1)	

## flip.traps
par(mfrow = c(1,2), xpd = TRUE)
traps1 <- make.grid(nx = 8, ny = 6, ID = "numxb")
traps2 <- flip (traps1, lr = TRUE)
plot(traps1, border = 5, label = TRUE, offset = 7, gridl = FALSE)
plot(traps2, border = 5, label = TRUE, offset = 7, gridl = FALSE)
par(mfrow = c(1,1), xpd = FALSE)

## rotate.traps
hollow1 <- make.grid(nx = 8, ny = 8, hollow = TRUE)
nested <- rbind (hollow1, rotate(hollow1, 45, c(70, 70)))
plot(nested, gridlines = FALSE)

## shift.traps
hollow1 <- make.grid(nx = 8, ny = 8, hollow = TRUE)
hollow2 <- shift(make.grid(nx = 6, ny = 6, hollow = TRUE), c(20, 20))
nested <- rbind (hollow1, hollow2)
plot(nested, gridlines = FALSE, label = TRUE)

}


