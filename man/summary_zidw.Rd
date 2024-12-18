\name{summary.zidw}
\title{Object Summaries}
\alias{summary.zidw}
\alias{print.summary.zidw}
\usage{
\method{summary}{zidw}(object, ...)
\method{print}{summary.zidw}(x, digits = max(3, getOption("digits") - 3), ...)
}
\description{
summary is a generic function used to produce result summaries of the results of various model fitting functions. The function invokes particular methods which depend on the class of the first argument.
}

\arguments{
  \item{object}{an object for which a summary is desired.}
  \item{...}{additional arguments affecting the summary produced.}
  \item{x}{a \code{"summary.zidw”} object.}
  \item{digits}{the number of digits in the output.}
} 
\details{
 Additional information about the ZIDW fit represented by \code{object} is extracted and included as components of
  \code{object}.  The returned object has a \code{\link{print}} method.
}





\examples{ 
## data
data("rwm1984", package = 'COUNT')

## model with covariates
zidw_out <- zidw_reg(docvis ~ outwork + female + age, ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
                     
summary(zidw_out)
}

\keyword{file}


