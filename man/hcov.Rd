\name{hcov}
\alias{hcov}

\title{ Hybrid Mixture Model }

\description{

The argument \code{hcov} in \code{secr.fit} is used to fit a hybrid
mixture model. `Hybrid' refers to a flexible combination of latent
classes (as in a finite mixture) and known classes (cf groups or
sessions). A hybrid mixture model includes a parameter `pmix' for the
mixing proportion and optionally allows detection parameters to be
modelled as class-specific ( ~ h2). This is particularly useful for
modelling sex ratio and sex differences in detection, and matches the
Bayesian sex-specific model of Gardner et al. (2010).

For observed animals all of unknown class the model is identical to a
finite mixture (i.e. latent-class) model. For observed animals all of
known class, the classes are no longer `latent' and the model is
equivalent to a grouped model with an additional binomial factor for
class membership.

}

\section{Assumptions}{

  \code{hcov} identifies a single individual covariate (the class
  covariate) that should be a factor with two levels, or contain
  character values that will be coerced to a factor (e.g., `f',
  `m'). Missing values (NA) are used for individuals of unknown
  class. If \code{hcov} has more than two levels, all but the first two
  levels are converted to NA (but see exception for h3 models below).

  It is assumed that the probability of recording a missing value for
  the class covariate is independent of the true class membership (e.g.,
  sex equally likely to be recorded for males and females).

}

\section{Operational details}{

  A hybrid mixture model is fitted whenever \code{hcov} is not
  NULL. Mixture models include a parameter `pmix', the mixing
  proportion. If the covariate identified by \code{hcov} is missing (`'
  or NA) for all individuals \emph{and} a mixture term (h2 or h3)
  appears in the detection model (e.g., g0 ~ h2) then a conventional
  finite mixture model is fitted (cf Pledger 2000, Borchers & Efford
  2008).

  As with finite mixture models, any detection parameter (g0, sigma
  etc.) may be modelled as depending on mixture class by model
  specifications such as (g0 ~ h2, sigma ~ h2). See Examples.

  In general \code{hcov} has been designed for two classes and two
  classes are assumed if neither `h2' nor `h3' appears in the model
  formulae. However, there is a small exception: \code{hcov} may have
  three non-missing levels if `h3' appears in a model formula. Note
  that h2 cannot be combined with h3; h3 is for advanced use only and
  has not been fully tested.

  The number of fitted parameters is the same as the corresponding
  finite mixture model if mixture terms (`h2', `h3') appear in the model
  formulae. Otherwise (no mixture terms) estimating pmix requires a
  single extra parameter. The estimate of pmix then depends solely on
  the observed class proportions in the covariate, and the beta
  variance-covariance matrix will show zero covariance of pmix with
  other detection parameters.

}

\section{Models for pmix}{

  Variation in the parameter pmix may be modelled across sessions i.e.,
  models such as pmix ~ session or pmix ~ Session are valid, as are
  formulae involving session covariates defined in the sessioncov
  argument of secr.fit.

  If no mixture term appears in the formula for pmix then one is added
  automatically (usually `h2'). This serves mostly to keep track of
  values in the output.

  Attempting to model pmix as a function of individual covariates or
  other within-session terms (t, b etc.) will cause an error.
}

\section{Interpreting output}{

  When you display a fitted secr model the parameter estimates are in a
  final section headed `Fitted (real) parameters evaluated at base
  levels of covariates'. The same output may be obtained by calling the
  \code{predict} method directly. Calling \code{predict} has the advantage
  that you can obtain estimates for levels of the covariates other than
  the base levels, by specifying \code{newdata}. An example below shows
  how to specify h2 in \code{newdata}. [Note: \code{predict}
  is generic, and you must consult ?predict.secr to see the help for the
  specific implementation of this method for fitted secr objects].

  The output from \code{predict.secr} for a mixture model is a list with
  one component for each (possibly latent) class. Each row corresponds
  to a fitted real parameter: ordinarily these include the detection
  parameters (e.g., g0, sigma) and the mixing proportion (pmix).

  In the case of a model fitted by maximizing the full likelihood
  (\code{CL = FALSE}), density D will also appear in the output. Note
  that only one parameter for density is estimated, the total density
  across classes. This total density figure appears twice in the
  output, once for each class.

  The standard error (SE.estimate) is shown for each parameter. These
  are asymptotic estimates back-transformed from the link scale. The
  confidence limits are also back-transformed from the link scale (95\%
  CI by default; vary \code{alpha} in \code{predict.secr} if you want
  e.g. 90\% CI).

  The mixing proportion pmix depends on the composition of the sample
  with respect to \code{hcov} and the detection model. For a null
  detection model the mixing proportion is exactly the proportion in the
  sample, with appropriate binomial confidence limits. Otherwise, the
  mixing proportion adjusts for class differences in the probability and
  scale of detection (see Examples).

  The preceding refers to the default behaviour when pmix ~ h2. It is
  possible also to fix the mixing proportion at any arbitrary value
  (e.g., fixed = list(pmix = 0.5) for 1:1 sex ratio).

  On output the classes are tagged with the factor levels of \code{hcov},
  regardless of how few or how many individuals were actually of known
  class. If only a small fraction were of known class, and there is
  cryptic variation unrelated to \code{hcov}, then the association
  between the fitted classes and the nominal classes (i.e. levels of
  \code{hcov}) may be weak, and should not be trusted.
  
}

\section{Limitations}{

  Hybrid mixture models are incompatible with groups as presently
  implemented.

  The hcov likelihood conditions on the number of known-class
  individuals. A model fitted with \code{hcov = NULL} or with a
  different hcov covariate has in effect a different data set, and
  likelihoods, deviances or AICs cannot be compared. AIC can be used to
  compare models provided they all have the same hcov covariate in the
  call to \code{secr.fit}, whether or not h2 appears in the model
  definition.

}

\section{Likelihood}{
  
The likelihood of the hybrid mixture model is detailed in an appendix of
the vignette \href{https://www.otago.ac.nz/density/pdfs/secr-finitemixtures.pdf}{secr-finitemixtures.pdf}.

}

\references{

Borchers, D.L. and Efford, M.G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics} \bold{64},
377--385.

Gardner, B., Royle, J.A., Wegan, M.T., Rainbolt, R. and Curtis,
P. (2010) Estimating black bear density using DNA data from hair
snares. \emph{Journal of Wildlife Management} \bold{74}, 318--325.

Pledger, S. (2000) Unified maximum likelihood estimates for
closed capture--recapture models using mixtures. \emph{Biometrics} \bold{56},
434--442.

}

\seealso{
  \code{\link{secr.fit} }
}

\examples{

\dontrun{

## house mouse dataset, morning trap clearances
## 81 female, 78 male, 1 unknown
morning <- subset(housemouse, occ = c(1,3,5,7,9))
summary(covariates(morning))

## speedy model fitting with coarse mask
mmask <- make.mask(traps(morning), buffer = 20, nx = 32)

## assuming equal detection of males and females
## fitted sex ratio p(female) = 0.509434 = 81 / (81 + 78)
fit.0 <- secr.fit(morning, hcov = "sex", mask = mmask, trace = FALSE)
predict(fit.0)

## allowing sex-specific detection parameters
## this leads to new estimate of sex ratio 
fit.h2 <- secr.fit(morning, hcov = "sex", mask = mmask, trace = FALSE,
    model = list(g0 ~ h2, sigma ~ h2))
predict(fit.h2)

## specifying newdata for h2 - equivalent to predict(fit.h2)
predict(fit.h2, newdata = data.frame(h2 = factor(c('f','m'))))

## conditional likelihood fit of preceding model
## estimate of sex ratio does not change 
fit.CL.h2 <- secr.fit(morning, hcov = "sex", mask = mmask, trace = FALSE,
    CL = TRUE, model = list(g0 ~ h2, sigma ~ h2))
predict(fit.CL.h2)

## did sexes differ in detection parameters?
fit.CL.0 <- secr.fit(morning, hcov = "sex", mask = mmask, trace = FALSE,
    CL = TRUE, model = list(g0 ~ 1, sigma ~ 1))
LR.test(fit.CL.h2, fit.CL.0)

## did sex ratio deviate from 1:1?
fit.CL.h2.50 <- secr.fit(morning, hcov = "sex", mask = mmask, trace = FALSE,
    CL = TRUE, model = list(g0 ~ h2, sigma ~ h2), fixed = list(pmix = 0.5))
LR.test(fit.CL.h2, fit.CL.h2.50)

## did sexes show extra-compensatory variation in lambda0?
## (Efford and Mowat 2014)
fit.CL.a0 <- secr.fit(morning, hcov = "sex", mask = mmask, trace = FALSE,
    CL = TRUE, model = list(a0 ~ 1, sigma ~ h2))
LR.test(fit.CL.h2, fit.CL.a0)

## trend in ovenbird sex ratio, assuming sex-specific detection
omask <- make.mask(traps(ovenCH), buffer = 300, nx = 32)
fit.sextrend <- secr.fit(ovenCH, model = list(g0~h2, sigma~h2, pmix~Session),
    hcov = "Sex", CL = TRUE, mask = omask, trace = FALSE)
predict(fit.sextrend)[1:5]

}

}

