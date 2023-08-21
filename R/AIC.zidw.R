AIC.zidw = function(object, ..., k = 2)
{
  -2 * object$loglik + 2 * length(unlist(coef(object)))
}


