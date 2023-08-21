dhdw <- function (x, q_par, beta, lam, log = FALSE) 
{
  
  maxl <- max(length(x), length(q_par), length(beta), length(lam))
  x <- rep(x, length.out = maxl)
  q <- rep(q_par, length.out = maxl)
  beta <- rep(beta, length.out = maxl)
  lam <- rep(lam, length.out = maxl)
  if (any(q <= 0 | q >= 1)) {
    stop(paste("q must be between 0 and 1!", "\n"))
  }
  if (any(beta <= 0)){
    stop(paste("beta must be greater than 0!", "\n"))
  }
  if (any(lam < 0 | lam >= 1)){
    stop(paste("lambda must be between 0 and 1!", "\n"))
  }
  x <- floor(x)
  delta <- ifelse(x == 0, 1, 0)
  # p <- (lam^delta) * ((1 - lam) * ((q^x^beta - q^(x + 1)^beta) / q))^(1 - delta)
  p <- (ddw2(x, q = q, beta = beta) / (1 - pdw2(0, q = q, beta = beta))) * lam
  p[x == 0] <- 1 - rep_len(lam, length.out = length(p))[x == 0]
  if (log) 
    p <- log(p)
  p[is.nan(p)] <- 0
  if (!log) 
    p <- pmin(pmax(p, 0), 1)
  p
}



