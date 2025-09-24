\name{summary.zidw}
\title{Object Summaries}
\alias{summary.zidw}
\alias{print.summary.zidw}
\usage{
\method{summary}{zidw}(object, ...)
\method{print}{summary.zidw}(x, digits = max(3, getOption("digits") - 3), ...)
}
\description{
\code{summary} method to produce results for objects of class \code{"zidw"}. 
}

\arguments{
  \item{object}{An object of class \code{"zidw"} for which a summary is desired.}
  \item{...}{Additional arguments affecting the summary produced.}
  \item{x}{A \code{summary.zidw} object.}
  \item{digits}{The number of digits in the output.}
} 
\details{
 Additional information about the ZIDW fit represented by \code{object} is extracted and included as components of
  \code{object}.  The returned object has a \code{\link{print}} method.
}
\value{
  Currently, limited results from an object of class \code{"zidw"} is returned.  The output is of class \code{summary.zidw}.
  Future versions of this code will attempt to mimic \code{summary.glm}.
}

 \seealso{
  \code{\link{summary.glm}}
} 


\examples{ 
## data
data("rwm1984", package = 'COUNT')

## model with covariates
zidw_out <- zidw_reg(docvis ~ outwork + female + age, 
                     ziformula = ~ age + female, betaformula = ~ 1, 
                     data = rwm1984[1:100, ], beta = .1)
                     
summary(zidw_out)
}

\keyword{file}


