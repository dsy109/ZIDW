zidw_r_squared <- function(object, adj = TRUE){
  
  d_zidw <- function(x, q, beta, lam) {
    lam * (x == 0) + (1 - lam) * ddw2(x, q = q, beta = beta)
  }
  
  dw_out <- dw.parest(object$response, method = 'likelihood')
  
  ll_dw <- sum(log(ddw2(object$response, q = dw_out$q, beta = exp(object$coefficients$beta))))
  
  # ll_full <- logLik(zidw_out)[1]
  # ll_red <- ll_dw
  
  obj <- function(q,x,N,beta,lambda) x-(1-lambda)*sum(q^((1:N)^beta))
  
  tmp <- sapply(1:object$nall,function(i) uniroot(obj,interval=c(0,1),x=object$response[i],N=300,beta=exp(object$coefficients$beta),lambda=as.numeric(object$response==0)[i])$root)
  
  ll_zidw_sat <- sum(log(d_zidw(object$response, q = tmp, beta = exp(object$coefficients$beta), as.numeric(object$response == 0))))
  
  
  ##
  if(adj == TRUE){
    R2_zidw_adj <- 1 - ((ll_zidw_sat - logLik(object)[1] + length(unlist(object$coefficients)) - 3 + 1.5) / (ll_zidw_sat - ll_dw))
  }else{
    R2_zidw_adj <- 1 - ((ll_zidw_sat - logLik(object)[1]) / (ll_zidw_sat - ll_dw))
  }
  return(R2_zidw_adj)
  
}
