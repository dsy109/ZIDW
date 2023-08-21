qhdw <- function(p, q_par, beta, lam, lower.tail = TRUE, log.p = FALSE) 
{
  ly <- max(length(p), length(q_par), length(beta), length(lam))
  p <- rep(p, length = ly)
  lam <- rep(lam, length = ly)
  q <- rep(q_par, length = ly)
  beta <- rep(beta, length = ly)
  
  if (any(q <= 0 | q >= 1)) {
    stop(paste("q must be between 0 and 1!", "\n"))
  }
  if (any(beta <= 0)){
    stop(paste("beta must be greater than 0!", "\n"))
  }
  if (any(lam < 0 | lam >= 1)){
    stop(paste("lambda must be between 0 and 1!", "\n"))
  }
  if (log.p) 
    p <- exp(p)
  
  if (lower.tail == TRUE){
    p <- p
  } else{
    p <- 1 - p
  }
  
  pnew <- log(pmax(0, (p + lam - 1))) - log(lam)
  # print(pnew)
  val <- qztdw(pnew, q_par = q, beta = beta, log.p = TRUE)
  #val <- sapply(1:length(pnew), function(i) qztdw(pnew[i], q = q[i], beta = beta[i]))
  val[!is.finite(pnew)] <- 0
  val[lam < 0 | lam > 1] <- NaN
  val
}

