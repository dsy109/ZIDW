pdw2 <- function (x, q = exp(-1), beta = 1) 
{
  maxl <- max(length(x), length(q), length(beta))
  x <- rep(x, length.out = maxl)
  q <- rep(q, length.out = maxl)
  beta <- rep(beta, length.out = maxl)
  
  if (any(q <= 0 | q >= 1)) {
    stop(paste("q must be between 0 and 1!", "\n"))
  }
  if (any(beta <= 0)) {
    stop(paste("beta must be greater than 0!", "\n"))
  }
  
  res <- 1 - q^((floor(x) + 1)^(beta))
  return(res)
}
