qzidw <- function(p, q_par, beta, lam, lower.tail = TRUE, log.p = FALSE) 
{
  if(any(q_par <= 0 | q_par >= 1)) {
    stop(paste("q must be between 0 and 1!", "\n"))
  }
  if(any(beta <= 0)){
    stop(paste("beta must be greater than 0!", "\n"))
  }
  if(any(lam < 0 | lam >= 1)){
    stop(paste("lambda must be between 0 and 1!", "\n"))
  }
  if (log.p) 
    p <- exp(p)

  if (lower.tail == TRUE){
    p <- p
  } else{
    p <- 1 - p
  }
  ly <- max(length(p), length(q), length(beta), length(lam))
  p <- rep(p, length = ly)
  lam <- rep(lam, length = ly)
  q <- rep(q_par, length = ly)
  beta <- rep(beta, length = ly)
  pnew <- ((p - lam)/(1 - lam)) - (1e-07)
  pnew <- ifelse(pnew > 0, pnew, 0)
  q <- sapply(1:length(pnew), function(i) qdw(pnew[i], q = q[i], beta = beta[i]))
  q
}







