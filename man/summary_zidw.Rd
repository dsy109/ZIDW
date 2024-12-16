\name{summary.zidw}
\title{Object Summaries}
\alias{summary.zidw}
\usage{
\method{summary}{zidw}(object, ...)
}
\description{

}

\arguments{
  \item{object}{an object for which a summary is desired.}
  \item{...}{additional arguments affecting the summary produced.}
} 
\details{
 The functions \code{summary.lm} and \code{summary.glm} are examples of particular methods which summarize the results produced by lm and glm.
}





\examples{ 
## data
data("rwm1984", package = 'COUNT')

## model with covariates
zidw_out <- zidw_reg(docvis ~ outwork + female + age, ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
                     
summary(zidw_out)
}

\keyword{file}


