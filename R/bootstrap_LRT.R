bootstrap_lrt <- function(data, B, tol = -1){
  
  y <- data[, 1]
  colnames(data) <- c('y')
  ## fit dw
  dw_fit <- dw.parest(y, method = 'likelihood')
  ## fit zidw
  zidw_fit <- zidw_reg(qformula = y ~ 1, data = data)
  
  ## calculate loglik
  dw_loglik <- sum(log(ddw2(y, q = dw_fit$q, beta = dw_fit$beta)))
  lrs <- max(-2 * (dw_loglik - zidw_fit$loglik), 0)
  
  ##### bootstrap
  lrs_bs <- c()
  count <- 0
  i = 0
  while (i <= B) {
    i = i + 1
    bs <- rdweibull(dim(data)[1], q = dw_fit$q, beta = dw_fit$beta, zero = TRUE)
    
    bs_dat <- data.frame(y = bs)
    ## fit dw
    dw_fit_bs <- dw.parest2(bs, method = 'likelihood')
    
    ## fit zidw
    zidw_fit_bs <- zidw_reg(qformula = y ~ 1, data = bs_dat, 
                                q = inv.logit(zidw_fit$coefficients$q), beta = exp(zidw_fit$coefficients$beta), lam = inv.logit(zidw_fit$coefficients$zero))
    
    
    h0 <- try(ddw2(bs, q = dw_fit_bs$q, beta = dw_fit_bs$beta), silent = TRUE)
    ha <- try(zidw_fit_bs$loglik, silent = TRUE)
    
    
    if(inherits(zidw_fit_bs, "try-error", which = TRUE)){
      i = i - 1
    }else if(any(is.na(h0)) | !inherits(ha, "numeric", which = FALSE) | is.na(ha)){
      i = i - 1
    }else{
      dw_loglik <- sum(log(h0))
      lrs_bs[i] <- -2 * (dw_loglik - zidw_fit_bs$loglik)
    
    }
    
    if(lrs_bs[i] < tol){
      i = i - 1
      count <- count + 1
    }
    
    if(between(lrs_bs[i], tol, 0)){
      lrs_bs[i] = 0
    }
    
  }
  ## p-value
  pval <- mean(lrs_bs > lrs)
  
  return_list <- list(pvalue = pval, 'Observe likelihood ratio test statistics' = lrs, 'Bootstrap likelihood ratio test statistics' = lrs_bs, 
                      count = count)
  return(return_list)
}







