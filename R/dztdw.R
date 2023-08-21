dztdw <- function (x, q_par, beta, log = FALSE) 
{
  
  maxl <- max(length(x), length(q), length(beta))
  x <- rep(x, length.out = maxl)
  q <- rep(q_par, length.out = maxl)
  beta <- rep(beta, length.out = maxl)

  if (any(q <= 0 | q >= 1)) {
    stop(paste("q must be between 0 and 1!", "\n"))
  }
  if (any(beta <= 0)){
    stop(paste("beta must be greater than 0!", "\n"))
  }

  x <- floor(x)
  # delta <- ifelse(x == 0, 1, 0)
  p <- rep(0, length.out = maxl)
  ind <- which(x > 0)
  if(length(ind) > 0){
    p[ind] <- (q[ind]^(x[ind]^beta[ind]) - q[ind]^((x[ind] + 1)^beta[ind])) / q[ind]
    }
  

  
  
  if (log) 
    p <- log(p)
  p[is.nan(p)] <- 0
  if (!log) 
    p <- pmin(pmax(p, 0), 1)
  p
}





