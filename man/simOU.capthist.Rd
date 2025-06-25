\name{simOU.capthist}
\alias{simOU}
\alias{simOU.capthist}

\title{Simulate Detections For Ornstein-Uhlenbeck Movement}

\description{

The usual SECR model ignores the sequential locations of an individual within 
its home range. Movement models predict serial correlation of detections in space. 
The Ornstein-Uhlenbeck (OU) model is a convenient example that over long 
durations leads to a bivariate normal distribution of locations.

Movements of a single animal according to the OU model are simulated in discrete time with \code{simOU}. 

Detections of a population of individuals with pre-defined activity centres are 
simulated with \code{simOU.capthist}. Detection happens when the location of an 
individual at time \eqn{t} (occasion \eqn{t}) is within distance epsilon \eqn{\epsilon} of a detector.

}

\usage{

simOU(xy, tau, sigma, noccasions, start = NULL)

simOU.capthist(traps, popn, detectpar, noccasions, seed = NULL,
    savepopn = FALSE, savepath = FALSE, ...)

}

\arguments{

  \item{xy}{numeric vector of x,y coordinates for one animal}
  \item{tau}{numeric parameter of serial correlation = \eqn{1 / \beta}}
  \item{sigma}{numeric spatial scale parameter}
  \item{noccasions}{integer number of time steps}
  \item{start}{numeric vector of x,y coordinates for initial location (optional)}

  \item{traps}{secr traps object}
  \item{popn}{secr popn object or a 2-column matrix of x-y coordinates of activity centres}
  \item{detectpar}{list with values of detection parameters epsilon, sigma, and tau}
  \item{seed}{either NULL or an integer that will be used in a call to \code{set.seed}}
  \item{savepopn}{logical; if TRUE the population is saved as an attribute}
  \item{savepath}{logical; if TRUE the movement paths are saved as an attribute}
  \item{\dots}{other arguments passed to \code{\link[secr]{reduce.capthist}}}
}

\details{

The first location for \code{simOU} by default is drawn at random from the asymptotic distribution.

The detection parameters are:

\tabular{ll}{
epsilon \tab radius within which individual detected with certainty \cr
sigma \tab spatial scale of asymptotic bivariate normal\cr
tau \tab serial correlation parameter 1/beta \cr
}
In a final step, `simOU.capthist' uses the reduce method for capthist objects to 
coerce the simulated capthist object to the detector type of the traps argument. 

The \dots argument may be used to pass the 'by' argument to reduce.capthist. For
example, 'by = "ALL"' collapses the initially binary data for a single detector 
on noccasions to a single integer. Alternatively, 'by = 10' collapses the 
original occasions in groups of 10. Data will be lost unless the input traps 
object has detector type 'count'.

The default `verify' step in reduce.capthist may also be suppressed 
(set \code{verify = FALSE} in \code{simOU.capthist()}).

The x- and y-dimensions are simulated separately, assuming circularity. 
The distribution of location on the \eqn{x} axis at time \eqn{t+1} 
conditional on the location at time \eqn{t} is then
\deqn{x_{t+1} | x_{t} \sim N(\mu_x + e^{\frac{-1}{\tau}(x_{t} - \mu_x)}, \; \sigma^2 (1 - e^{\frac{-2}{\tau}}),} 
where \eqn{\mu = (\mu_x, \mu_y)} is the long-term activity centre and 
\eqn{\tau} (tau) is a parameter for the strength of serial correlation 
(\eqn{\tau = 1/\beta} in other formulations). 
The scale of the long-term (asymptotic) bivariate normal home range is governed
by \eqn{\sigma} as usual. Steps are implicitly of length 1 occasion so 
\eqn{\Delta t} is omitted.

}

\value{

simOU - matrix of locations dim = c(noccasions, 2)

simOU.capthist - single-session capthist object for \pkg{secr}

}

\seealso{

\code{\link[secr]{sim.capthist}}

}

\examples{

# one animal
locs <- simOU(c(0,0), 20, 1, 100)
par(pty = 's')
plot(locs, type = 'o', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5))
points(0,0, pch = 16, col = 'red')

\donttest{

# simulate some capture data
set.seed(123)
grid <- make.grid(8, 8, spacing = 2)
pop <- sim.popn(D = 1000, core = grid, buffer = 4)
ch <- simOU.capthist(grid, pop, detectpar=list(tau = 50, sigma = 1, 
    epsilon = 0.25), noccasions = 100, savepath = TRUE)

# plot simulated capthist (blue dots), movements (grey), and AC (triangles)
plot(pop, pch = 2, cex = 1, border = 4)
plot(ch, rad = 0.2, tracks = TRUE, varycol = FALSE, add = TRUE,
    gridlines = FALSE)
symbols(grid$x, grid$y, circles = rep(0.25,64), inches = FALSE, fg = 'red', 
    add = TRUE)
sapply(attr(ch, 'path'), lines, col = 'grey')

# fit a misspecified HHN model
# density estimates show a small negative bias 
# sigma estimates show a large negative bias 
# lambda0 cannot be compared
fit <- secr.fit(capthist = ch, buffer = 8, detectfn = 14, trace = FALSE)
predict(fit)

}

}