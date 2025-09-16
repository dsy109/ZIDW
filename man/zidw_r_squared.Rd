\name{zidw_r_squared}
\title{Pseudo-\eqn{R^2} for ZIDW Regression}
\alias{zidw_r_squared}
\usage{
zidw_r_squared(object, adj = TRUE)
}
\description{
 Calculates the pseudo-\eqn{R^2} for ZIDW regression fits. Modeling \eqn{q} and \eqn{\lambda}.
}

\arguments{
  \item{object}{An object from a ZIDW regression fit of class \code{"zidw"}.}
  \item{adj}{An adjustment for calculating the pseudo-\eqn{R^2}. Default is \code{adj = TRUE}.}
 
} \value{
  A numerical value with the corresponding pseudo-\eqn{R^2}.
} 

\references{
  Martin, J. and Hall, D. B. (2016), \eqn{R^2} Measures for Zero-Inflated Regression Models for Count Data with Excess Zeros. \emph{Journal of Statistical Computation and Simulation}, \bold{84}:18, 3777--3790.
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


