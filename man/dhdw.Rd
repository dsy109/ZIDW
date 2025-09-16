\name{hdw}
\title{Hurdle Discrete Weibull Distribution}
\alias{hdw}
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
Density, distribution function, quantile function and random generation for the hurdle Discrete Weibull distribution with parameter \eqn{q_{par}}, \eqn{\beta} and hurdle crossing probability \eqn{\lambda} (i.e., \eqn{1 - \lambda} is the probability for observed zeros).
}

\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observation.}
  \item{q_par}{Shape parameter.}
  \item{beta}{Shape parameter.}
  \item{lam}{Zero-infation parameter.}
  \item{log, log.p}{Logical; if \code{TRUE}, probabilities are returned on log-scale.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
} 
\details{
The hurdle discrete Weibull has the mass function

\eqn{f(y) = \lambda + (1 - \lambda)\frac{q^{y^\beta} - q^{(y + 1)^\beta}}{q}},

for \eqn{y=0, 1, 2, ...} where \eqn{\lambda \in (0, 1)} is the zero-inflation parameter, \eqn{q \in (0, 1)}, and \eqn{\beta > 0}.
}
  \value{
    \item{dhdw}{Gives the density.}
    \item{phdw}{Gives the cumulative probability.}
    \item{qhdw}{Gives the quantile value.}
    \item{rhdw}{Generates random numbers.}
  
} 

\references{
  Kalktawi, H. S. (2017), \emph{Discrete Weibull Regression Model for Count Data}. Ph.D. Thesis, Brunel University London.

  Yeh, P. and Young, D. S. (2025), Some Estimation and Inference Considerations for the Zero-Inflated Discrete Weibull Distribution. \emph{Communications in Statistics - Simulation and Computation (in press)}, 1--22.
} 



\examples{ 
dhdw(.5, .6, 1, .4)
phdw(.5, .6, 1, .4)
qhdw(.8, .6, 1, .4)
rhdw(100, .6, 1, .4)
}

\keyword{file}


