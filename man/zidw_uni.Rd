\name{zidw_uni}
\title{Parameter Estimation for the Zero-Inflated Discrete Weibull}
\alias{zidw_uni}
\usage{
zidw_uni(y, par = NULL, method = c('MLE', 'mde'), B = NULL, max_method = NULL, 
         constraint = TRUE)
}
\description{
  Parameter estimation for the univariate zero-inflated discrete Weibull distribution using maximum likelihood or minimum distance estimation.
}

\arguments{

  \item{y}{Vector of univariate counts.}
  \item{par}{Vector of starting values in the order of \eqn{\lambda}, \eqn{\beta}, and \eqn{q}.}
  \item{method}{Parameter estimation method.  Methods include \code{uni_method = 'MLE'} (default) 
  for maximum likelihood estimation and \code{uni_method = 'MDE'} for minimum distance estimation.}
  \item{B}{Number of bootstrap replciates for standard error estimation when modeling univariate data using \code{uni_method = 'MDE'}.}
  \item{max_method}{Maximum likelihood estimation method for the univariate setting. See \link[maxLik:maxLik]{maxLik} for the different methods.}
  \item{constraint}{Logical value to constrain the parameter to be within the parameter space.  Default is \code{constraint = TRUE}.}
} \value{
  Output for for maximum likelihood estimation (\code{uni_method = 'MLE'}):
  \item{MLE}{Maximum likelihood estimates.}
  \item{coefficients}{Coefficients for \eqn{q}, \eqn{\beta}, and \eqn{\lambda}.}
  \item{convergence}{Convergence code of \code{maxLik}.}
  \item{iteration}{Number of iterations.}
  \item{SE}{Estimated standard errors.}
  Output for for minimum distance estimation (\code{uni_method = 'MDE'}):
  \item{lam}{Parameter estimate for \eqn{\lambda}.}
  \item{Beta}{Parameter estimate for \eqn{\beta}.}
  \item{q}{Parameter estimate for \eqn{q}.}
  \item{distance}{The final distance value.}
  \item{SE}{Bootstrap standard error estimates for \eqn{\lambda}, \eqn{\beta}, and \eqn{q}.}
} \seealso{
  \code{\link{glm}}, \code{\link{lm}}
} 


\references{
  Taconeli, C. A. and Rodrigues de Lara, I. A. (2022), Discrete Weibull Distribution: Different Estimation Methods Under Ranked Set Sampling and Simple Random Sampling. \emph{Journal of Statistical Computation and Simulation}, \bold{92}:8, 1740--1762.

  Yeh, P. and Young, D. S. (2025), Some Estimation and Inference Considerations for the Zero-Inflated Discrete Weibull Distribution. \emph{Communications in Statistics - Simulation and Computation (in press)}, 1--22.
} 



\examples{ 
## data
data("rwm1984", package = 'COUNT')

## MLE method
zidw_out <- zidw_uni(rwm1984$docvis)
                     
\dontrun{                     
## mde without standard error estimation
zidw_out2 <- zidw_uni(rwm1984$docvis, method = 'mde')

## mde with bootstrap standard error estimation
zidw_out2 <- zidw_uni(rwm1984$docvis, method = 'mde', B = 5)
}

}

\keyword{file}


