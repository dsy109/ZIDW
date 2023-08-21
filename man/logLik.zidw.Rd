\name{logLik.zidw}
\title{Extract Log-Likelihood for ZIDW}
\alias{logLik.zidw}
\usage{
\method{logLik}{zidw}(object, ...)
}
\description{
  Extract the log-likelihood of ZIDW
}

\arguments{
  \item{object}{a fitted ZIDW model object for which there exists a logLik method to extract the corresponding log-likelihood, or an object inheriting from class logLik.}
  \item{...}{some methods for this generic function require additional arguments.}
} 
  \value{
  a numeric value with the corresponding log-likelihood
  
} \seealso{
  \code{\link{logLik}} \code{\link{AIC}}
} 




\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)

logLik(zidw_out)
}

\keyword{file}


