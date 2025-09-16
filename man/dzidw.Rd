\name{zidw}
\title{Zero-Inflated Discrete Weibull Distribution}
\alias{zidw}
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
Density, distribution function, quantile function and random generation for the zero-inflated Discrete Weibull distribution with parameters \eqn{q_{par}}, \eqn{\beta}, \eqn{\lambda}.
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
The zero-inflated discrete Weibull has the mass function

\eqn{f(y) = \lambda + (1 - \lambda)(q^{y^\beta} - q^{(y + 1)^\beta})},

for \eqn{y = 0, 1, 2, ...} where \eqn{\lambda \in (0, 1)} is zero-inflation parameter, \eqn{q \in (0, 1)}, and \eqn{\beta > 0}. 
}
  \value{
    \item{dzidw}{Gives the density.}
    \item{pzidw}{Gives the cumulative probability.}
    \item{qzidw}{Gives the quantile value.}
    \item{rzidw}{Generates random numbers.}
  
} 

\references{
  Kalktawi, H. S. (2017), \emph{Discrete Weibull Regression Model for Count Data}. Ph.D. Thesis, Brunel University London.

  Yeh, P. and Young, D. S. (2025), Some Estimation and Inference Considerations for the Zero-Inflated Discrete Weibull Distribution. \emph{Communications in Statistics - Simulation and Computation (in press)}, 1--22.
} 


\examples{ 
dzidw(.5, .6, 1, .4)
pzidw(.5, .6, 1, .4)
qzidw(.8, .6, 1, .4)
rzidw(100, .6, 1, .4)
}

\keyword{file}


