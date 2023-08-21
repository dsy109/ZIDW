ddw2 <- function (x, q = exp(-1), beta = 1) 
{
  LLL <- max(length(x), length(q), length(beta))
  x <- rep(x, length.out = LLL)
  q <- rep(q, length.out = LLL)
  beta <- rep(beta, length.out = LLL)
  
  if (any(q > 1 | q < 0)) 
    stop("q must be between 0 and 1", call. = FALSE)
  if (any(beta <= 0))
    stop("beta must be positive", call. = FALSE)
  # is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - 
  #                                                                    round(x)) < tol

    
  res <- q^(x^(beta)) - q^((x + 1)^(beta))
  
  return(res)
}
