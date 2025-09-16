\name{predict.zidw}
\title{Prediction for ZIDW Model Fits}
\alias{predict.zidw}
\usage{
\method{predict}{zidw}(object, newdata, 
                       type = c("response", "prob", "count", "zero"), 
                       at = NULL, ...)
}
\description{
Obtains predictions from the fitted ZIDW model object.
}

\arguments{
  \item{object}{A fitted object of class inheriting from \code{"zidw"}.}
  \item{newdata}{Optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted predictors are used.}
  \item{type}{The type of prediction required. For details see below.}
  \item{at}{Optionally, if \code{type = "prob"}, a numeric vector at which the probabilities are evaluated. By default 0:max(y) is used where y is the original observed response.}
  \item{...}{Currently not used.}
} 
\details{
 The default is \code{type = "response"}, which on the scale of the response variable. 
 \code{type = "prob"} provides the predicted density (i.e., probabilities for the observed counts). 
 \code{type = "count"} offers the predicted mean from the count component without zero-inflation and \code{type = "zero"} predicts the probability for the zero-component.
}





\examples{ 
## data
data("rwm1984", package = 'COUNT')

## model with covariates
zidw_out <- zidw_reg(docvis ~ outwork + female + age, 
                     ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
                     
predict.zidw(zidw_out)
}

\keyword{file}


