\name{secr-defunct}
\alias{secr-defunct}
\alias{ip.secr}
\alias{pfn}
\alias{make.newdata}
\alias{secr.make.newdata}
\alias{read.SPACECAP}
\alias{write.SPACECAP}

\title{Defunct Functions in Package \pkg{secr}}

\description{

These functions are no longer available in \pkg{secr}.

}

\usage{

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

\code{ip.secr} and \code{pfn} have been superceded by \code{ipsecr.fit} and 
\code{proxy.ms} in package \pkg{ipsecr}.

Internal functions \code{secr.make.newdata} and \code{make.newdata} were 
replaced with \code{makeNewData} generic from 4.4.2.

SPACECAP was removed from the CRAN archive on 2019-08-31.

}

\seealso{

\code{\link{secr-deprecated}}

}
