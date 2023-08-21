rhdw <- function (n, q_par, beta, lam) 
{
  if (q_par <= 0 | q_par >= 1) {
    stop(paste("q must be between 0 and 1!", "\n"))
  }
  if(beta <= 0){
    stop(paste("beta must be greater than 0!", "\n"))
  }
  if(lam < 0 | lam >= 1){
    stop(paste("lambda must be between 0 and 1!", "\n"))
  }
  n <- ceiling(n)
  p <- runif(n)
  r <- qhdw(p, q_par = q_par, beta = beta, lam = lam)
  r
}


