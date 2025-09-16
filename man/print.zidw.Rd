\name{print.zidw}
\title{Print Values from a ZIDW Object}
\alias{print.zidw}
\usage{
\method{print}{zidw}(x, digits = max(3, getOption("digits") - 3), ...)
}
\description{
  Print the output of a ZIDW object.
}

\arguments{
  \item{x}{A ZIDW object used to select a method.}
  \item{digits}{Minimal number of significant digits.}
  \item{...}{Further arguments passed to or from other methods.}
} 





\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, 
                     ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)

print(zidw_out)
}

\keyword{file}


