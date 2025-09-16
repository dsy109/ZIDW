\name{rqres_zidw_reg}
\title{Randomized Quantile Residuals for a ZIDW Regression Fit}
\alias{rqres_zidw_reg}
\usage{
rqres_zidw_reg(test, plot = FALSE)
}
\description{
  Calculate randomized quantile residuals from a ZIDW regression fit.
}

\arguments{
  \item{test}{Zero-inflated discrete Weibull regresssion output.}
  \item{plot}{If \code{plot = TRUE}, produce the randomized quantile residuals plot.}
} 
  \value{
  Return a vector of randomized quantile resuduals or a randomized quantile residuals plot:
  \item{rqr}{Randomized quantile resuduals.}
  
} \seealso{
  \code{\link{glm}}, \code{\link{lm}}
} 

\references{
  Dunn, P. K. and Smyth, G. K. (1996), Randomzied Quantile Residuals. \emph{Journal of Computational and Graphical Statistics}, \bold{5}:3, 236--244.

  Yeh, P. and Young, D. S. (2025), Some Estimation and Inference Considerations for the Zero-Inflated Discrete Weibull Distribution. \emph{Communications in Statistics - Simulation and Computation (in press)}, 1--22.
} 

\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, 
                     ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
rqr <- rqres_zidw_reg(zidw_out)
}

\keyword{file}


