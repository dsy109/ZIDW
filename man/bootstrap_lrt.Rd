\name{bootstrap_lrt}
\title{Bootstrap Likelihood Ratio Test for Zero-Inflation}
\alias{bootstrap_lrt}
\usage{
bootstrap_lrt(data, B, tol = -1)
}
\description{
  Perform the bootstrap likelihood ratio test for comparing DW and ZIDW model fits.
}

\arguments{
  \item{data}{Data to test.}
  \item{B}{Number of bootstrap samples to draw.}
  \item{tol}{Threshold of bootstrap likelihood ratio test statistics. See details.}
} 
\details{
  Theoretically, the likelihood ratio is always positive. In practice, however, if a negative (but close to zero) likelihood ratio occurrs, then the calculations are typically numerically unstable or the fitting procedure failed to converge to the MLEs. We, therefore, set the threshold to -1 by default. 
  If the bootstrap likelihood ratio is less than \code{tol}, we drop this sample and resample. If it is between \code{tol} and 0, we truncate it at 0.
}
  \value{
  Return a list containing four vectors:
  \item{pvalue}{P-value of the test.}
  \item{Observe likelihood ratio test statistics}{Likelihood ratio test statistics from the data.}
  \item{Bootstrap likelihood ratio test statistics}{Bootstrap likelihood ratio test statistics.}
  \item{count}{Number of times that bootstrap LRT is less than \code{tol}.}
  
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


