
rqres_zidw_reg <- function(test, plot = FALSE)
{
  beta0 <- test$coefficients$beta
  q0 <- test$coefficients$q
  lam0 <- test$coefficients$zero
  y <- test$response
  n <- length(y)
  X <- test$model_matrix_q
  Z <- test$model_matrix_beta
  W <- test$model_matrix_zi
  q <- inv.logit(X %*% q0)
  beta <- exp(Z %*% beta0)
  lam <- inv.logit(W %*% lam0)
  
  # rqres <- function(y, F, eps = 1e-6)
  # {
  #   n <- length(y)
  #   FL <- F(y - eps)
  #   FU <- F(y)
  #   u <- runif(n, min = FL, max = FU)
  #   qres <- qnorm(u)
  #   return(qres)
  # }
  
  rqres <- function(y, F, eps = 1e-6, ...)
  {
    n <- length(y)
    FL <- pmin(F(pmax(y - eps,0), ...),1-eps)
    FU <- F(y, ...)
    u <- runif(n, min = FL, max = FU)
    qres <- qnorm(u)
    return(qres)
  }
  
  F <- function(y) {
    ret <- numeric(n)
    for (i in 1:n) {
      ret[i] <- pzidw(y[i], q_par = q[i], beta = beta[i], lam = lam[i])
    }
    return(ret)
  }
  
  rval <- rqres(y, F)
  
  if(plot == TRUE){
    df_zidw <- data.frame('qqzidw' = rval)
    ggplot(df_zidw, aes(sample = qqzidw))+stat_qq()+
      #  geom_abline(intercept = 0, slope = 1,col=2)+
      stat_qq_line(col=2)+
      ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
      ggtitle("Randomized Quantile Residuals (ZIDW Fit)")+
      theme(text = element_text(size = 15))
  }else{
    rval
  }
  

}
