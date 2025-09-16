\name{ztdw}
\title{Zero-Truncated Discrete Weibull Distribution}
\alias{ztdw}
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
Density, distribution function, quantile function and random generation for the zero-truncated Discrete Weibull distribution with parameter \eqn{q_{par}}, \eqn{\beta}.
}

\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observation.}
  \item{q_par}{Shape parameter.}
  \item{beta}{Shape parameter.}
  \item{log, log.p}{Logical; if \code{TRUE}, probabilities are returned on log-scale.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
} 
\details{
The zero-truncated discrete Weibull has the mass function

\eqn{f(y) = \frac{q^{y^\beta} - q^{(y + 1)^\beta}}{q}},

for \eqn{y = 1, 2, 3, ...} where \eqn{q \in (0, 1)}, and \eqn{\beta > 0}.
}
  \value{
    \item{dztdw}{Gives the density.}
    \item{pztdw}{Gives the cumulative probability.}
    \item{qztdw}{Gives the quantile value.}
    \item{rztdw}{Generates random numbers.}
  
} 

\references{
  Kalktawi, H. S. (2017), \emph{Discrete Weibull Regression Model for Count Data}. Ph.D. Thesis, Brunel University London.

  Yeh, P. and Young, D. S. (2025), Some Estimation and Inference Considerations for the Zero-Inflated Discrete Weibull Distribution. \emph{Communications in Statistics - Simulation and Computation (in press)}, 1--22.
} 


\examples{ 
dztdw(1, .6, 1)
pztdw(1, .6, 1)
qztdw(.8, .6, 1)
rztdw(100, .6, 1)
}

\keyword{file}


