\name{AIC.zidw}
\title{AIC and BIC for ZIDW Model Fits}
\alias{AIC.zidw}
\alias{BIC.zidw}
\usage{
\method{AIC}{zidw}(object, ..., k = 2)
\method{BIC}{zidw}(object, ...)
}
\description{
  Generic function calculating AIC or BIC for ZIDW model fits.
}

\arguments{
  \item{object}{A fitted ZIDW model object for which there exists a \code{logLik} method to extract the corresponding log-likelihood, or an object inheriting from class \code{logLik}.}
  \item{...}{Optionally more fitted model objects.}
  \item{k}{Numeric value for the penalty per parameter to be used; default = 2.}
} 
  \value{
  A numeric value with the corresponding AIC or BIC.
  
} \seealso{
  \code{\link{logLik}} \code{\link{AIC}}
} 


\references{
Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986), \emph{Akaike Information Criterion Statistics}. D. Reidel Publishing Company.
} 

\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, 
                     ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)

AIC(zidw_out)
BIC(zidw_out)
}

\keyword{file}


