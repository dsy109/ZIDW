\name{AIC.zidw}
\title{Akaike's An Information Criterion for ZIDW}
\alias{AIC.zidw}
\alias{BIC.zidw}
\usage{
\method{AIC}{zidw}(object, ..., k = 2)
\method{BIC}{zidw}(object, ...)
}
\description{
  Generic function calculating Akaike's 'An Information Criterion' for one or several fitted model objects for which a log-likelihood value can be obtained, according to the formula -2*log-likelihood + k*npar, where npar represents the number of parameters in the fitted model, and k = 2 for the usual AIC, or k = log(n) (n being the number of observations) for the so-called BIC or SBC (Schwarz's Bayesian criterion).
}

\arguments{
  \item{object}{a fitted ZIDW model object for which there exists a logLik method to extract the corresponding log-likelihood, or an object inheriting from class logLik.}
  \item{...}{optionally more fitted model objects}
  \item{k}{numeric, the penalty per parameter to be used; default = 2}
} 
  \value{
  a numeric value with the corresponding AIC
  
} \seealso{
  \code{\link{logLik}} \code{\link{AIC}}
} 


\references{
Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986). Akaike Information Criterion Statistics. D. Reidel Publishing Company.
} 

\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)

AIC(zidw_out)
BIC(zidw_out)
}

\keyword{file}


