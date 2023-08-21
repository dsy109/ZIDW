\name{zidw_r_squared}
\title{Calculating the pseudo-\eqn{R^2}}
\alias{zidw_r_squared}
\usage{
zidw_r_squared(object, adj = TRUE)
}
\description{
 Calculating the pseudo-\eqn{R^2} for ZIDW regression. Modeling q and \eqn{\lambda}.
}

\arguments{
  \item{object}{an object from ZIDW regression}
  \item{adj}{adjustment for pseudo-\eqn{R^2}. Assume to be TRUE.}
 
} \value{
  a numerical value with corresponding pseudo-\eqn{R^2}
} 

\references{
  Jacob Martin and Daniel B. Hall. R2 measures for zero-inflated regression models for count data with excess zeros. Journal of Statistical Computation and Simulation, 86:18, 3777-3790, DOI: 10.1080/00949655.2016.1186166
} 

\examples{ 
## data
data("rwm1984", package = 'COUNT')

## model with covariates
zidw_out <- zidw_reg(docvis ~ outwork + female + age, ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
                     
zidw_r_squared(zidw_out)
}

\keyword{file}


