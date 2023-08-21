\name{zidw_reg}
\title{Fitting Zero-Inflated Discrete Weibull Regression}
\alias{zidw_reg}
\usage{
zidw_reg(qformula, betaformula = ~ 1, ziformula = ~ 1, data, lam = NULL, beta = NULL, 
         q = NULL, k = 1000, uni_method =  c('MLE', 'mde'), 
         max_method = NULL, constraint = TRUE, B = NULL)
}
\description{
  Fit zero-inflated regression model for count data via \code{constrOptim}. Fit univariate model for count data via maximum likelihood and minimum distance.
}

\arguments{
  \item{qformula}{a symbolic description of the model to be fitted for parameter q including the response variable.}
  \item{betaformula}{a symbolic description of the model to be fitted for parameter \eqn{\beta}.}
  \item{ziformula}{a symbolic description of the model to be fitted for parameter \eqn{\lambda}.}
  \item{data}{data frame containg model variables.}
  \item{lam}{starting value for parameter \eqn{\lambda}.}
  \item{beta}{starting value for parameter \eqn{\beta}}
  \item{q}{starting value for parameter q.}
  \item{k}{cutofff value to calculate residuals.}
  \item{uni_method}{univariate fitting method. MLE by default.}
  \item{max_method}{maximum likelihood estimate method for univariate setting. NM by default. See \code{maxLik}.}
  \item{constraint}{FALSE for unconstrained maximization. TRUE by default.}
  \item{B}{number of bootstrap sample size for standard error estimates when modeling univariate data using MDE.}
} \value{
  \code{zidw_reg} An object of class "zidw", i.e., a list with components including
  \item{call}{the original function call.}
  \item{coefficients}{a list with elements "zero", "beta", and "q" containing the coefficients from the respective models}
  \item{loglik}{log-likehood of the fitted model}
  \item{SE}{standard error}
  \item{convergence}{convergence code of \code{optim}}
  \item{nall}{number of observations}
  \item{res}{a vector of raw residuals (observed - fitted)}
  \item{fitted_values}{fitted values}
  \item{model_matrix_q}{model matrix of parameter q}
  \item{model_matrix_beta}{model matrix of parameter \eqn{\beta}}
  \item{model_matrix_zi}{model matrix of parameter \eqn{\lambda}}
  \item{response}{vector of response variable}
  \item{model}{dataset}
  \item{formula}{a list contains three formulas for q, \eqn{\beta}, and \eqn{\lambda}}
} \seealso{
  \code{\link{glm}}, \code{\link{lm}}
} 



\examples{ 
## data
data("rwm1984", package = 'COUNT')

## model with covariates
zidw_out <- zidw_reg(docvis ~ outwork + female + age, ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
                     
## model without covariates
zidw_out2 <- zidw_reg(docvis ~ 1, data = rwm1984)

}

\keyword{file}


