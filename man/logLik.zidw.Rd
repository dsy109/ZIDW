\name{logLik.zidw}
\title{Extract Log-Likelihood for ZIDW Model Fits}
\alias{logLik.zidw}
\usage{
\method{logLik}{zidw}(object, ...)
}
\description{
  Extract the log-likelihood of a ZIDW model fit.
}

\arguments{
  \item{object}{A fitted ZIDW model object for which there exists a \code{logLik} method to extract the corresponding log-likelihood, or an object inheriting from class \code{logLik}.}
  \item{...}{Some methods for this generic function require additional arguments.}
} 
  \value{
  A numeric value with the corresponding log-likelihood.
  
} \seealso{
  \code{\link{logLik}} \code{\link{AIC}}
} 




\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, 
                     ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)

logLik(zidw_out)
}

\keyword{file}


