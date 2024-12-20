\name{blackbear}
\docType{data}

\alias{blackbear}
\alias{blackbearCH}
\alias{blackbear.0}
\alias{blackbear.h2bk}
\alias{GSM}

\title{ Black Bear Hair Snag Dataset }

\description{

American black bears \emph{Ursus americanus} were surveyed with baited hair 
snags in the Great Smoky Mountains National Park, Tennessee, in the summer of 2003.

}

\usage{

blackbearCH
GSM
blackbear.0
blackbear.h2bk

}

\details{

American black bears \emph{Ursus americanus} were surveyed in the Tennessee 
sector of Great Smoky Mountains National Park over  9 June--15 August 2003. 
Baited hair snags (barbed wire enclosures) were operated for 10 weeks at 65 sites, 
about 1 km apart and mostly close to trails. Bait consisted of bakery products 
in a small waxed-paper bag. Raspberry extract was used as a scent lure.

Genotyping and non-spatial capture-recapture analysis of a data subset were
described by Settlage et al. (2008). The sex of each genotyped bear was 
determined subsequently and some additional samples were included (J. Laufenberg 
pers. comm. 2012-05-09; Laufenberg et al. 2013). 

The dataset is a single-session capthist object with binary proximity 
detector type. Snags were visited weekly, so there were 10 occasions 
in the raw data. A single covariate `sex' was recorded for each individual.

The dataset comprises 282 detections of 81 females and 58 males. Female 15 
apparently made a long movement (17 km) between occasions 1 and 3.

The hair snag array sampled less than 20\% of the area of the park. The 
unforested area outside the park on the northwestern boundary of the study area 
was not considered to be black bear habitat (F. van Manen pers. comm. 2012-05-18) 
and should be excluded in analyses. The approximate boundary of the park is 
included as a shapefile `GSMboundary.shp' in the `extdata' folder of the package 
and as the \pkg{sf} sfc_POLYGON object \code{GSM}. The latter may be used in 
\code{make.mask} (see Examples).

Two models (blackbear.0 and blackbear.h2bk) were fitted as shown in the Examples.

}

\source{

The data were provided by Jared Laufenberg, Frank van Manen and Joe Clark. 

}

\references{

Laufenberg, J. S., van Manen, F. T., and Clark, J. D.  (2013)  Effects of sampling conditions on DNA-based estimates of American black bear abundance.  \emph{Journal of Wildlife Management} \bold{77}, 1010--1020. 

Settlage, K. E., Van Manen, F. T., Clark, J. D., and King, T. L. (2008) Challenges of DNA-based mark--recapture studies of American black bears. \emph{Journal of Wildlife Management} \bold{72}, 1035--1042.

}

\examples{

summary(blackbearCH)

\donttest{

# GSM is the approximate boundary of Great Smoky Mountains National Park
# Make a habitat mask restricted to the park

tr <- traps(blackbearCH)
msk <- make.mask(tr, buffer = 6000, type = 'trapbuffer', poly = GSM)

# Plot

plot(GSM)
plot(msk, add = TRUE)
plot(blackbearCH, tracks = TRUE, add = TRUE)
plot(tr, add = TRUE)

# Fit models
# suppress fastproximity to allow learned response

setNumThreads()   # as appropriate

# null model
blackbear.0 <- secr.fit(blackbearCH,  detectfn = 'EX', hcov = 'sex', 
    mask = msk, details = list(fastproximity = FALSE), trace = FALSE)

# sex differences and site-specific behavioural response
blackbear.h2bk <- secr.fit(blackbearCH,  detectfn = 'EX', hcov = 'sex', 
    model = list(g0~bk+h2, sigma~h2), mask = msk, 
    details = list(fastproximity = FALSE), trace = FALSE)
    
AIC(blackbear.0, blackbear.h2bk)
summary(blackbear.h2bk)   

# How many if we extrapolate to GSM NP?
region.N(blackbear.h2bk, region = GSM)

}


}

