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

\value{
  Returns an object of class \code{zidw} based on the object \code{x} with mostly the same output when printing output from
  an object of class \code{lm}.  Whereas an \code{lm} object has one \code{model.matrix} in its output,
  use of \code{print.zidw} returns three model matrices, one each for the model parameters \eqn{q} and \eqn{beta} as well as for the
  zero-inflation component.
}
 \seealso{
  \code{\link{print.lm}}
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


