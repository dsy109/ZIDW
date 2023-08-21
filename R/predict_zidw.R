predict.zidw <- function (object, newdata, type = c("response", "prob", "count", 
                                    "zero"), at = NULL, ...) 
{
  type <- match.arg(type)
  
  
  if (missing(newdata)) {
    rval <- object$fitted_values
    if (type != "response") {
      
      X <- object$model_matrix_q
      W <- object$model_matrix_beta
      Z <- object$model_matrix_zi

      q <- inv.logit(X %*% object$coefficients$q)[, 1]
      beta <- exp(W %*% object$coefficients$beta)[, 1]
      lam <- inv.logit(Z %*% object$coefficients$zero)[, 1]
      
      # mu <- exp(X %*% object$coefficients$count)[, 1]
      
      y <- seq(1, 1000, 1)
      mu <- sapply(1:length(q), function(i) sum(q[i]^y^beta[i]))
      phi <- lam
    }
  }

  else {
    mf <- model.frame(object$model, newdata) ## full model
    
    X <- model.matrix(object$formula$qformula, mf)
    W <- model.matrix(object$formula$betaformula, mf)
    Z <- model.matrix(object$formula$ziformula, mf)
    
    q <- inv.logit(X %*% object$coefficients$q)[, 1]
    beta <- exp(W %*% object$coefficients$beta)[, 1]
    lam <- inv.logit(Z %*% object$coefficients$zero)[, 1]
    
    y <- seq(1, 1000, 1)
    mu <- sapply(1:length(q), function(i) sum(q[i]^y^beta[i]))
    phi <- lam
    rval <- (1 - phi) * mu
  }
  

  
  if (type == "count") 
    rval <- mu
  
  if (type == "zero") 
    rval <- phi
  
  if (type == "prob") {
    if (!is.null(object$response)) 
      y <- object$response
    else if (!is.null(object$model))   ## extract dataset
      y <- model.response(object$model)
    else stop("predicted probabilities cannot be computed for fits with y = FALSE and model = FALSE")
    yUnique <- if (is.null(at))  
      0:max(y)
    else at
    nUnique <- length(yUnique)
    rval <- matrix(NA, nrow = length(rval), ncol = nUnique)
    dimnames(rval) <- list(rownames(X), yUnique)
    
    
    rval[, 1] <- phi + (1 - phi) * exp(-mu)
    for (i in 2:nUnique) rval[, i] <- (1 - phi) * dzidw(yUnique[i], q_par = q, beta = beta, lam = lam)
    
  }
  rval
}
