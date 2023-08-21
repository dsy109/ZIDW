BIC.zidw <- function(object, ...)
{
  n <- object$nall
  -2 * object$loglik + log(n) * length(unlist(coef(object)))
}


