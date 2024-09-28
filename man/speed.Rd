\name{speed}
\alias{speed}
\alias{Speed tips}

\title{ Speed Tips }

\description{
  
  A list of ways to make \code{\link{secr.fit}} run faster. 
  
}

\section{Use an appropriate mask}{

Check the extent and spacing of the habitat mask that you are using.
Execution time is roughly proportional to the number of mask points
(\code{nrow(mymask)}).  Default settings can lead to very large masks
for detector arrays that are elongated `north-south' because the number
of points in the east-west direction is fixed. Compare results with a
much sparser mask (e.g., nx = 32 instead of nx = 64).

}

\section{Use conditional likelihood}{

If you don't need to model variation in density over space or time then
consider maximizing the conditional likelihood in secr.fit (CL =
TRUE). This reduces the complexity of the optimization problem,
especially where there are several sessions and you want
session-specific density estimates (by default, derived() returns a
separate estimate for each session even if the detection parameters are
constant across sessions).
}

\section{Model selection}{

Do you really need to fit all those complex models? Chasing down small
decrements in AIC is so last-century. Remember that detection parameters
are mostly nuisance parameters, and models with big differences in AIC
may barely differ in their density estimates. This is a good topic for
further research - we seem to need a `focussed information criterion'
(Claeskens and Hjort 2008) to discern the differences that matter. Be
aware of the effects that can really make a difference: learned
responses (b, bk etc.) and massive unmodelled heterogeneity.

Use score.test() to compare nested models. At each stage this
requires only the more simple model to have been fitted in full; further
processing is required to obtain a numerical estimate of the gradient of
the likelihood surface for the more complex model, but this is much
faster than maximizing the likelihood. The tradeoff is that the score
test is only approximate, and you may want to later verify the results
using a full AIC comparison.

}

\section{Break problem down}{

  Suppose you are fitting models to multiple separate datasets that fit
  the general description of `sessions'. If you are fitting separate
  detection parameters to each session (i.e., you do not need to pool
  detection information), and you are not modelling trend in density
  across sessions, then it is much quicker to fit each session
  separately than to try to do it all at once. See Examples.
  
}

\section{Mash replicated clusters of detectors}{

If your detectors are arranged in similar clusters (e.g., small square
grids) then try the function \code{\link{mash}}.

}

\section{Reduce sparse `proximity' data to `multi'}{

  Full data from `proximity' detectors has dimensions n x S x K (n is
  number of individuals, S is number of occasions, K is number of
  traps). If the data are sparse (i.e. multiple detections of an
  individual on one occasion are rare) then it is efficient to treat
  proximity data as multi-catch data (dimension n x S, maximum of one
  detection per occasion). Use \code{reduce(proxCH, outputdetector =
  "multi")}.
  
}

\section{ Use multiple cores when applicable }{
  
  Most computers these days have multiple processors and these will be used by \pkg{secr} if the user sets \code{ncores} greater than one in \code{secr.fit}, \code{sim.secr} and some other functions. If \code{ncores = NULL} then the existing value from the environment variable RCPP_PARALLEL_NUM_THREADS is used (see \code{\link{setNumThreads}}).

}

\section{ Avoid covariates with many levels }{

  Categorical (factor) covariates with many levels and continuous
  covariates that take many values are not handled efficiently in
  secr.fit, and can dramatically slow down analyses and increase memory
  requirements.

}

\section{ Set appropriate typsize }{

Setting typsize manually in the call of `secr.fit` can speed up fitting when 
magnitudes on the link scale are very different (for example, when an identity 
link is used for density and density is very small or very large).

}

\section{ Simulations }{

Model fitting is not needed to assess power. The precision of estimates
  from secr.fit can be predicted without laboriously fitting models to
  simulated datasets. Just use \code{method = "none"} to obtain the asymptotic
  variance at the known parameter values for which data have been
  simulated (e.g. with sim.capthist()).

Suppress computation of standard errors by derived(). For a
  model fitted by conditional likelihood (CL = TRUE) the subsequent
  computation of derived density estimates can take appreciable time. If
  variances are not needed (e.g., when the aim is to predict the bias of
  the estimator across a large number of simulations) it is efficient to
  set se.D = FALSE in derived().

It is tempting to save a list with the entire `secr' object from
each simulated fit, and to later extract summary statistics as
needed. Be aware that with large simulations the overheads associated
with storage of the list can become very large. The solution is to
anticipate the summary statistics you will want and save only these.

}

\references{
  
Claeskens, G. and Hjort N. L. (2008) \emph{Model Selection and Model
Averaging}. Cambridge: Cambridge University Press.

}

\examples{

\dontrun{

## compare timing of combined model with separate single-session models
## for 5-session ovenbird mistnetting data: 2977/78 = 38-fold difference

setNumThreads(7)

system.time(fit1 <- secr.fit(ovenCH, buffer = 300, trace = FALSE, 
    model = list(D ~ session, g0 ~ session, sigma ~ session)))
##    user  system elapsed 
## 1837.71   31.81  730.56 

system.time(fit2 <- lapply (ovenCH, secr.fit, buffer = 300, trace = FALSE))
##   user  system elapsed 
##  43.74    0.46   11.13  

## ratio of density estimates
collate(fit1)[,1,1,"D"] / sapply(fit2, function(x) predict(x)["D","estimate"])
## session=2005 session=2006 session=2007 session=2008 session=2009 
##    1.0000198    1.0000603    0.9999761    0.9999737    0.9999539 

}

}


