\name{secr-defunct}
\alias{secr-defunct}
\alias{model.average}
\alias{ip.secr}
\alias{pfn}
\alias{make.newdata}
\alias{secr.make.newdata}
\alias{read.SPACECAP}
\alias{write.SPACECAP}
\alias{derivedSystematic}

\title{Defunct Functions in Package \pkg{secr}}

\description{

These functions are no longer available in \pkg{secr}.

}

\usage{

# Defunct in 5.0.0

derivedSystematic()

# Defunct in 4.6.2 (2023-09-30)

model.average()

# Defunct in 4.5.10 (2023-03-10)

ip.secr()
pfn()

# Defunct in 4.4.2 (2021-05-04)

make.newdata()

# Defunct in 4.4.0 (2021-05-01)

secr.make.newdata()

# Defunct in 4.0.0 (2019-10-27)

read.SPACECAP()
write.SPACECAP()

}

\details{

Some of these have stubs which report that they are defunct, but most have been 
removed completely (apart from being documented here). 

\code{model.average} is replaced by a method for `secr' and `secrlist' objects of the generic \code{\link{modelAverage}}. The internal code is essentially the same for \code{model.average} and \code{modelAverage.secrlist}. The generic avoids a name conflict with RMark and is also used in \pkg{openCR}.

\code{ip.secr} and \code{pfn} have been superceded by \code{ipsecr.fit} and 
\code{proxy.ms} in package \pkg{ipsecr}.

Internal functions \code{secr.make.newdata} and \code{make.newdata} were 
replaced with \code{makeNewData} generic from 4.4.2.

SPACECAP was removed from the CRAN archive on 2019-08-31.

}

\seealso{

\code{\link{secr-deprecated}}

}
