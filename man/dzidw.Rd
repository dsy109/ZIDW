\name{ZIDW}
\title{Zero-Inflated Discrete Weibull Distribution}
\alias{ZIDW}
\alias{rzidw}
\alias{qzidw}
\alias{pzidw}
\alias{dzidw}
\usage{
dzidw(x, q_par, beta, lam, log = FALSE)
pzidw(q, q_par, beta, lam, lower.tail = TRUE, log.p = FALSE)
qzidw(p, q_par, beta, lam, lower.tail = TRUE, log.p = FALSE)
rzidw(n, q_par, beta, lam)
}
\description{
Density, distribution function, quantile function and random generation for the zero-inflated Discrete Weibull distribution with parameter \eqn{q_{par}, \beta, \lambda}
}

\arguments{
  \item{x, q}{vector of quantiles}
  \item{p}{vector of probabilities}
  \item{n}{number of observation}
  \item{q_par}{shape parameter}
  \item{beta}{shape parameter}
  \item{lam}{zero-infation parameter}
  \item{log, log.p}{logical; if TRUE, probabilities are returned on log-scale}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}}
} 
\details{
The zero-inflated discrete Weibull has the mass function

\eqn{f(y) = \lambda + (1 - \lambda)(q^{y^\beta} - q^{(y + 1)^\beta})}

for \eqn{y = 0, 1, 2,...} where \eqn{\lambda \in (0, 1)} is zero-inflated parameter, \eqn{q \in (0, 1)}, and \eqn{\beta > 0} 
}
  \value{
    \item{dzidw}{gives the density}
    \item{pzidw}{gives the cumulative probability}
    \item{qzidw}{gives the quantile value}
    \item{rzidw}{generates random numbers}
  
} 




\examples{ 
dzidw(.5, .6, 1, .4)
pzidw(.5, .6, 1, .4)
qzidw(.8, .6, 1, .4)
rzidw(100, .6, 1, .4)
}

\keyword{file}


