\name{rqres_zidw_reg}
\title{Produce randomized quantile residuals}
\alias{rqres_zidw_reg}
\usage{
rqres_zidw_reg(test, plot = FALSE)
}
\description{
  Calculate randomized quantile resuduals for residuals diagnostic
}

\arguments{
  \item{test}{zero-inflated discrete Weibull regresssion output}
  \item{plot}{if TRUE, produce the randomized quantile residuals plot}
} 
  \value{
  Return a vector of randomized quantile resuduals or rqr plot
  \item{rqr}{randomized quantile resuduals}
  
} \seealso{
  \code{\link{glm}}, \code{\link{lm}}
} 



\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
rqr <- rqres_zidw_reg(zidw_out)
}

\keyword{file}


