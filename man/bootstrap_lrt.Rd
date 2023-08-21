\name{bootstrap_lrt}
\title{Perform bootstrap likelihood ratio test}
\alias{bootstrap_lrt}
\usage{
bootstrap_lrt(data, B, tol = -1)
}
\description{
  Perform bootstrap likelihood ratio test comparing DW and ZIDW
}

\arguments{
  \item{data}{data to test}
  \item{B}{number of bootstrap sample size}
  \item{tol}{threshold of bootstrap likelihood ratio test statistics. See details}
} 
\details{
  Theoretically, LR is always positive. In practice, however, if a negative (but close to zero) LR occurred, then the calculations are probably just numerically unstable, or the fitting procedure failed to converge to the MLEs. We, therefore, set the threshold to -1 by default. If the bootstrap LR is less than \code{tol}, we drop this sample and resample. If it is between \code{tol} and 0, we truncated it at 0
}
  \value{
  Return a list containing four vectors
  \item{pvalue}{p-value of the test}
  \item{Observe likelihood ratio test statistics}{likelihood ratio test statistics from the data}
  \item{Bootstrap likelihood ratio test statistics}{Bootstrap likelihood ratio test statistics}
  \item{count}{number of times that bootstrap LRT is less than \code{tol}}
  
} \seealso{
  \code{\link{glm}}, \code{\link{lm}}
} 


\examples{ 
\dontrun{
## data
data("rwm1984", package = 'COUNT')

test <- bootstrap_lrt(rwm1984, B = 1000)
}
}

\keyword{file}


