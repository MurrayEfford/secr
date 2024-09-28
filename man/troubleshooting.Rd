\name{Troubleshooting}
\alias{Troubleshooting}

\title{ Problems in Fitting SECR Models }

\description{

  Although \code{secr.fit} is quite robust, it does not always
work. Inadequate data or an overambitious model occasionally cause
numerical problems in the algorithms used for fitting the model, or
problems of identifiability, as described for capture--recapture models
in general by Gimenez et al. (2004). Here are some tips that may help
you.

This page has largely been superceded by \href{https://www.otago.ac.nz/density/pdfs/secr-troubleshooting.pdf}{secr-troubleshooting.pdf}.

}

\section{The log-likelihood values shown with trace = TRUE are all NA}{
  Most likely you have infeasible starting values for the
  parameters. try some alternatives, specifying them manually with the
  \code{start} argument.
}

\section{secr.fit finishes, but some or all of the variances are missing}{
  
  This usually means the model did not fit and the estimates should not
  be trusted. Extremely large variances or standard errors also indicate
  problems.
  
  \itemize{

    \item Try another maximization method (\code{method = "Nelder-Mead"}
  is more robust than the default). The same maximum likelihood should
  be found regardless of method, so AIC values are comparable across
  methods.

  \item Repeat the maximization with different starting values. You can use
\code{secr.fit(..., start = last.model)} where \code{last.model} is a
previously fitted secr object.

    \item If you think the estimates are right but there was a problem
    in computing the variances, try re-running secr.fit() with the
    previous model as starting values (see preceding) and set
    \code{method = "none"}. This bypasses maximization and computes the
    variances afresh using \code{fdHess} from \pkg{nlme}.
  
\item Try a finer mask (e.g., vary argument \code{nx} in
\code{\link{make.mask}}). Check that the extent of the mask matches your
data.

\item The maximization algorithms work poorly when the beta coefficients
are of wildly different magnitude. This may happen when using
covariates: ensure beta coefficients are similar (within a factor of
5--10 seems adequate, but this is not based on hard evidence) by scaling
any covariates you provide. This can be achieved by setting the
\code{typsize} argument of \code{nlm} or the \code{parscale}
control argument of \code{optim}.

\item Examine the model. Boundary values (e.g., g0 near 1.0) may give
problems. In the case of more complicated models you may gain insight by
fixing the value of a difficult-to-estimate parameter (argument
\code{fixed}).

}

See also the section `Potential problems' in
\href{https://www.otago.ac.nz/density/pdfs/secr-densitysurfaces.pdf}{secr-densitysurfaces.pdf}.

}

\section{secr.fit finishes with warning nlm code 3}{

  This condition does not invariably indicate a failure of model
  fitting. Proceed with caution, checking as suggested in the preceding
  section.
  
}

\section{secr.fit crashes part of the way through maximization}{

  A feature of the maximization algorithm used by default in \code{nlm}
  is that it takes a large step in the parameter space early on in the
  maximization. The step may be so large that it causes floating point
  underflow or overflow in one or more real parameters.  This can be
  controlled by passing the `stepmax' argument of \code{nlm} in the
  \dots argument of \code{secr.fit} (see first example). See also the
  previous point about scaling of covariates.
  
  }

\section{secr.fit demands more memory than is available}{

  This is a problem particularly when using individual covariates in a
  model fitted by maximizing the conditional likelihood. The memory
  required is then roughly proportional to the product of the number of
  individuals, the number of occasions, the number of detectors and the
  number of latent classes (for finite-mixture models).  When maximizing
  the full-likelihood, substitute `number of groups' for `number of
  individuals'. [The limit is reached in external C used for the
  likelihood calculation, which uses the R function `R_alloc'.]

  The \code{\link{mash}} function may be used to reduce the number of
  detectors when the design uses many identical and independent
  clusters. Otherwise, apply your ingenuity to simplify your model,
  e.g., by casting `groups' as `sessions'. Memory is less often an issue
  on 64-bit systems (see link below).

}

\section{Estimates from mixture models appear unstable}{ These models
  have known problems due to multimodality of the likelihood. See
  \href{https://www.otago.ac.nz/density/pdfs/secr-finitemixtures.pdf}{secr-finitemixtures.pdf}.
  }

\seealso{\code{\link{secr.fit}}, \code{\link[base]{Memory-limits}} }

\references{
  
Gimenez, O., Viallefont, A., Catchpole, E. A., Choquet, R. and Morgan, B. J. T. (2004)
Methods for investigating parameter redundancy. \emph{Animal Biodiversity and
Conservation} \bold{27}, 561--572.

}

