\name{coef.zidw}
\title{Print Coefficients from a ZIDW Model Fit}
\alias{coef.zidw}
\usage{
\method{coef}{zidw}(object, ...)
}
\description{
  Print the coefficients of a ZIDW object.
}

\arguments{
  \item{object}{A ZIDW object to extract the model coefficients.}
  \item{...}{Further arguments.}
} 

  \value{
  Coefficients extracted from the ZIDW object \code{object}. This will be a named numeric vector. 
}



\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, 
                     ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)

coef.zidw(zidw_out)
}

\keyword{file}


