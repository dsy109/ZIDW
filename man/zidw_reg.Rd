\name{zidw_reg}
\title{Fitting of Zero-Inflated Discrete Weibull Regression Models}
\alias{zidw_reg}
\usage{
zidw_reg(qformula, betaformula = ~ 1, ziformula = ~ 1, data, lam = NULL, 
         beta = NULL, q = NULL, k = 1000, uni_method =  c('MLE', 'mde'), 
         max_method = NULL, constraint = TRUE, B = NULL)
}
\description{
  Fitting of the zero-inflated discrete Weibull regression model is done via \code{constrOptim}. Fitting of the univariate model can be accomplished via maximum likelihood or minimum distance estimation.
}

\arguments{
  \item{qformula}{A symbolic description of the model to be fitted for the parameter \eqn{q}, including the response variable.}
  \item{betaformula}{A symbolic description of the model to be fitted for the parameter \eqn{\beta}.}
  \item{ziformula}{a symbolic description of the model to be fitted for the parameter \eqn{\lambda}.}
  \item{data}{A data frame containg model variables.}
  \item{lam}{The starting value for \eqn{\lambda}.}
  \item{beta}{The starting value for parameter \eqn{\beta}.}
  \item{q}{The starting value for parameter \eqn{q}.}
  \item{k}{A cutofff value used to calculate residuals.}
  \item{uni_method}{The univariate fitting method. Methods include \code{uni_method = 'MLE'} (default) 
  for maximum likelihood estimation and \code{uni_method = 'MDE'} for minimum distance estimation.}
  \item{max_method}{Maximum likelihood estimation method for the univariate setting. See \link[maxLik:maxLik]{maxLik} for the different methods.}
  \item{constraint}{Set \code{constraint = FALSE} for unconstrained maximization.  \code{constraint = TRUE} is the default.}
  \item{B}{Number of bootstrap replciates for standard error estimation when modeling univariate data using \code{uni_method = 'MDE'}.}
} \value{
  \code{zidw_reg} An object of class \code{"zidw"}, i.e., a list with components including the following:
  \item{call}{The original function call.}
  \item{coefficients}{A list with elements "zero", "beta", and "q" containing the coefficients from the respective models.}
  \item{loglik}{The log-likehood of the fitted model.}
  \item{SE}{Estimated standard errors.}
  \item{convergence}{The convergence code of \code{optim}.}
  \item{nall}{The number of observations.}
  \item{res}{A vector of raw residuals (observed - fitted).}
  \item{fitted_values}{The fitted values.}
  \item{model_matrix_q}{The model matrix of parameter \eqn{q}.}
  \item{model_matrix_beta}{The model matrix of parameter \eqn{\beta}.}
  \item{model_matrix_zi}{The model matrix of parameter \eqn{\lambda}.}
  \item{response}{Vector of response values.}
  \item{model}{The dataset.}
  \item{formula}{A list containing three formulas for \eqn{q}, \eqn{\beta}, and \eqn{\lambda}.}
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

## model with covariates
zidw_out <- zidw_reg(docvis ~ outwork + female + age, ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
                     
## model without covariates
zidw_out2 <- zidw_reg(docvis ~ 1, data = rwm1984)

}

\keyword{file}


