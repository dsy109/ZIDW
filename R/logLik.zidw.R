logLik.zidw <- function (object, ...) 
{
  p <- sum(lengths(object$coefficients))
  val <- object$loglik
  attr(val, "nall") <- object$nall
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}

