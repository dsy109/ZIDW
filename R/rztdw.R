rztdw <- function (n, q_par, beta) 
{
  if (any(q_par <= 0 | q_par >= 1)) {
    stop(paste("q must be between 0 and 1!", "\n"))
  }
  if (any(beta <= 0)){
    stop(paste("beta must be greater than 0!", "\n"))
  }

  n <- ceiling(n)
  p <- runif(n)
  r <- qztdw(p, q_par = q_par, beta = beta)
  r
}
