\name{zidw-internal}
\alias{dw.parest2}
\alias{ddw2}
\alias{pdw2}


\title{Internal 'ZIDW' Functions}
\description{
  Internal functions for the package \code{zidw}.
}
\usage{
dw.parest2(data, method = "likelihood", method.opt = "NR")
ddw2(x, q = exp(-1), beta = 1)
pdw2(x, q = exp(-1), beta = 1)
}

\arguments{
  \item{data}{Data that is a vector.}
  \item{method}{Either "likelihood" or "proportion"}
  \item{method.opt}{Optimization criterion used in maxLik (default is "NR")}
  \item{x}{vector of quantile}
  \item{q, beta}{Parameters of the distribution}
}

\details{
  These are usually not to be called by the user.
}


\keyword{internal}
