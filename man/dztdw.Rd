\name{ZTDW}
\title{Zero-Truncated Discrete Weibull Distribution}
\alias{ZTDW}
\alias{rztdw}
\alias{qztdw}
\alias{pztdw}
\alias{dztdw}
\usage{
dztdw(x, q_par, beta, log = FALSE)
pztdw(q, q_par, beta, lower.tail = TRUE, log.p = FALSE)
qztdw(p, q_par, beta, lower.tail = TRUE, log.p = FALSE)
rztdw(n, q_par, beta)
}
\description{
Density, distribution function, quantile function and random generation for the zero-truncated Discrete Weibull distribution with parameter \eqn{q_{par}, \beta}.
}

\arguments{
  \item{x, q}{vector of quantiles}
  \item{p}{vector of probabilities}
  \item{n}{number of observation}
  \item{q_par}{shape parameter}
  \item{beta}{shape parameter}
  \item{log, log.p}{logical; if TRUE, probabilities are returned on log-scale}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}}
} 
\details{
The zero-truncated discrete Weibull has the mass function

\eqn{f(y) = \frac{q^{y^\beta} - q^{(y + 1)^\beta}}{q}}

for \eqn{y = 1, 2, 3,...} where \eqn{q \in (0, 1)}, and \eqn{\beta > 0} 
}
  \value{
    \item{dztdw}{gives the density}
    \item{pztdw}{gives the cumulative probability}
    \item{qztdw}{gives the quantile value}
    \item{rztdw}{generates random numbers}
  
} 




\examples{ 
dztdw(1, .6, 1)
pztdw(1, .6, 1)
qztdw(.8, .6, 1)
rztdw(100, .6, 1)
}

\keyword{file}


