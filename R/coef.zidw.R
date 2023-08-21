coef.zidw <- function(object, ...)
{
  c(object$coefficients$q, object$coefficients$beta, object$coefficients$zero)
}
