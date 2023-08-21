

qztdw <- function(p, q_par, beta, lower.tail = TRUE, log.p = FALSE){
  
  if (lower.tail == FALSE){
    p <- 1 - p
  }
  
  if (any(q_par <= 0 | q_par >= 1)) {
    stop(paste("q must be between 0 and 1!", "\n"))
  }
  if (any(beta <= 0)){
    stop(paste("beta must be greater than 0!", "\n"))
  }
  
  ly <- max(length(p), length(q), length(beta))
  p <- rep(p, length = ly)
  q <- rep(q_par, length = ly)
  beta <- rep(beta, length = ly)
  
  p_orig <- p
  p <- if (log.p) 
    p
  else log(p)
  
  p <- p + log(1 - pdw2(0, q = q, beta = beta))
  p <- exp(p) + ddw2(0, q = q, beta = beta)
  
  ## lower.tail=TRUE
  # if(lower.tail){
  #   rval <- sapply(1:length(p), function(i) qdw(p[i], q = q[i], beta = beta[i]))
  # }else{
  #   rval <- sapply(1:length(p), function(i) qdw(1 - p[i], q = q[i], beta = beta[i]))
  # }
  rval <- sapply(1:length(p), function(i) qdw(p[i], q = q[i], beta = beta[i]))

  if(lower.tail){
    rval[p_orig < dztdw(1, q_par = q, beta = beta, log = log.p)] <- 1
  }
  rval
}
