\name{newdata}
\alias{newdata}
\alias{makeNewData}
\alias{makeNewData.secr}
\alias{makeNewData.default}

\title{
Create Default Design Data
}

\description{
  
Internal function used to generate a dataframe containing design data
for the base levels of all predictors in an secr object.

}

\usage{
\method{makeNewData}{default}(object, all.levels = FALSE, ...)
\method{makeNewData}{secr}(object, all.levels = FALSE, ...)

}

\arguments{
  \item{object}{fitted secr model object}
  \item{all.levels}{logical; if TRUE then all levels of factors are included}
  \item{\dots}{other arguments (not used)}
}

\details{

\code{makeNewData} is used by \code{predict} in lieu of
user-specified `newdata'. There is seldom any need to call the function 
\code{makeNewData} directly.

}

\value{

A dataframe with one row for each session and group, and columns for the
predictors used by \code{object$model}.

}

\seealso{

\code{\link{predict.secr}}, \code{\link{secr.fit}}

}

\examples{

## from previously fitted model
makeNewData(secrdemo.b)

}


