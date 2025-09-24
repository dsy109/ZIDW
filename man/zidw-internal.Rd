\name{zidw-internal}
\alias{dw.parest2}
\alias{ddw2}
\alias{pdw2}


\title{Internal 'ZIDW' Functions}
\description{
  Internal functions for the package \code{ZIDW}.
}
\usage{
dw.parest2(data, method = "likelihood", method.opt = "NR")
ddw2(x, q = exp(-1), beta = 1)
pdw2(x, q = exp(-1), beta = 1)
}

\arguments{
  \item{data}{Data that is a vector.}
  \item{method}{Either \code{"likelihood"} or \code{"proportion"}.}
  \item{method.opt}{Optimization criterion used in \code{maxLik} (default is \code{"NR"}).}
  \item{x}{Vector of quantiles.}
  \item{q, beta}{Parameters of the distribution.}
}

\value{
  These functions yield different outputs, which are called for their desired effects within
  other functions.  Thus, they are usually not to be called by the user.
}

\details{
  These are usually not to be called by the user.
}


\keyword{internal}
