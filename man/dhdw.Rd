\name{HDW}
\title{Hurdle Discrete Weibull Distribution}
\alias{HDW}
\alias{rhdw}
\alias{qhdw}
\alias{phdw}
\alias{dhdw}
\usage{
dhdw(x, q_par, beta, lam, log = FALSE)
phdw(q, q_par, beta, lam, lower.tail = TRUE, log.p = FALSE)
qhdw(p, q_par, beta, lam, lower.tail = TRUE, log.p = FALSE)
rhdw(n, q_par, beta, lam)
}
\description{
Density, distribution function, quantile function and random generation for the hurdle Discrete Weibull distribution with parameter \eqn{q_{par}, \beta} and hurdle crossing probability \eqn{\lambda} (i.e., \eqn{1 - \lambda} is the probability for observed zeros).
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
The hurdle discrete Weibull has the mass function

\eqn{f(y) = \lambda + (1 - \lambda)\frac{q^{y^\beta} - q^{(y + 1)^\beta}}{q}}

for \eqn{y=0, 1, 2,...} where \eqn{\lambda \in (0, 1)} is zero-inflated parameter, \eqn{q \in (0, 1)}, and \eqn{\beta > 0} 
}
  \value{
    \item{dhdw}{gives the density}
    \item{phdw}{gives the cumulative probability}
    \item{qhdw}{gives the quantile value}
    \item{rhdw}{generates random numbers}
  
} 




\examples{ 
dhdw(.5, .6, 1, .4)
phdw(.5, .6, 1, .4)
qhdw(.8, .6, 1, .4)
rhdw(100, .6, 1, .4)
}

\keyword{file}


