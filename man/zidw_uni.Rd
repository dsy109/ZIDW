\name{zidw_uni}
\title{Fitting Zero-Inflated Discrete Weibull}
\alias{zidw_uni}
\usage{
zidw_uni(y, par = NULL, method = c('MLE', 'mde'), B = NULL, max_method = NULL, 
         constraint = TRUE)
}
\description{
  Fit zero-inflated univariate model for count data via \code{maxLik} or \code{mde}.
}

\arguments{

  \item{y}{univariate response data.}
  \item{par}{starting values in the order of \eqn{\lambda}, \eqn{\beta}, q.}
  \item{method}{optimaziation method. MLE or mde. MLE by defalut.}
  \item{B}{number of bootstrap sample size for standard error estimation when mde is used. Ignore by deafult.}
  \item{max_method}{algorithm to use for MLE method. See \link{maxLik}}
  \item{constraint}{logical value which constraint the parameter to be within the parameter space}
} \value{
  Output for MLE
  \item{MLE}{maximum likelihood estimate}
  \item{coefficients}{coefficients for q, \eqn{\beta}, and \eqn{\lambda}}
  \item{convergence}{convergence code of \code{maxLik}}
  \item{iteration}{number of iteration}
  \item{SE}{estimated standard error}
  Output for mde
  \item{lam}{parameter estimation for \eqn{\lambda}}
  \item{Beta}{parameter estimation for \eqn{\beta}}
  \item{q}{parameter estimation for q}
  \item{distance}{the distance}
  \item{SE}{bootstrap standard error estimates for \eqn{\lambda}, \eqn{\beta}, and q}
} \seealso{
  \code{\link{glm}}, \code{\link{lm}}
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


