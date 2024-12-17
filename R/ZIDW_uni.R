

zidw_uni <- function(y, par = NULL, method = c('MLE', 'mde'), B = NULL, max_method = NULL, constraint = TRUE){
  method <- match.arg(method)
  n <- length(y)
  
  if(is.null(max_method)){
    max_method <- 'NM'
  }else if(is.null(max_method) & constraint == FALSE){
    max_method <- 'NR'
  }
  
  if(constraint == TRUE & max_method %in% c('NR', 'BHHH')){
    stop('Wrong maximization method')
  }else if(constraint == FALSE & max_method %in% c('BFGS', 'BFGSR', 'CG', "SANN", 'NM')){
    stop('Wrong maximization method')
  }

  if(is.null(par)){
    data1 <- y[which(y > 0)] - 1
    q <- dw.parest2(data1)$q
    lambda <- mean(y == 0)
    beta <- dw.parest2(data1)$beta
  }else{
    lambda <- par[1]
    beta <- par[2]
    q <- par[3] 
  }

  p <- rbinom(n, size = 1, prob = lambda)
  
  lam <- mean(p)
  if(method == 'mde'){
    ddw2 <- Vectorize(ddw)

    d2dw <- function(y, lam, Beta, q){
      
      return((y == 0) * (lam + (1 - q) * (1 - lam))  + (1 - (y == 0)) * (1 - lam) * ddw2(y, q, Beta))
    }
    
    p2dw <- function(y, lam, Beta, q){
      return(sapply(1:length(y),function(i) sum(d2dw(y = 0:y[i], Beta = Beta, q = q, lam = lam))))
    }
    
    p2dw(y = y, lam = lam, Beta = beta, q = q)
    
    out <- try(mde(y, p2dw, start = list(lam = lam, Beta = beta, q = q), 
               measure = "CvM", method = 'L-BFGS-B', lower = c(0, 0, 0), upper = c(1, Inf, 1)), silent = TRUE)
    
    ## bootstrap
    if(is.null(B)){
      return(out)
    } else{
      B <- B
      se <- matrix(NA, nrow = B, ncol = 3)

      i <- 1
      while (i < B+1) {
        
        y1 <- sample(y, n, replace = TRUE)
        
        out_se <- try(mde(y1, p2dw, start = list(lam = lam, Beta = beta, q = q),
                       measure = "CvM", method = 'L-BFGS-B', lower = c(0, 0, 0), upper = c(1, Inf, 1)), silent = TRUE)

        
        if(class(out_se)[1] == 'try-error') {
          i = i - 1
        } else{
          if(any(is.na(out_se$estimate))){
          i = i - 1
        } else {
          se[i, 1] <- out_se$estimate[1]
          se[i, 2] <- out_se$estimate[2]
          se[i, 3] <- out_se$estimate[3]
        }} 
        i <- i+1
        
      }
      
      return(list('out' = out, 'SE' = sqrt(diag(var(se)))))
      
    }


  }
  
  else if(method == 'MLE'){
    ll <- function(par, y){
      delta <- ifelse(y == 0, 1, 0)
      
      sum(delta * log(par[1] + (1 - par[1]) * (1 - par[3])) + (1 - delta) * log((1 - par[1]) * (par[3]^y^par[2] - par[3]^(y + 1)^par[2])))
      
    }
    
    ui <- rbind(c(1, 0, 0), c(-1, 0, 0), c(0, 1, 0), c(0, -1, 0), c(0, 0, 1), c(0, 0, -1))
    ci <- c(0, -1, 0, -100, 0, -1)
    
    ineqA = ui
    ineqB = -ci
    
    if(constraint == TRUE){
      out <- maxLik::maxLik(logLik = ll, y = y, start = c(lam, beta, q), constraints = list(ineqA = ineqA, ineqB = ineqB), method = max_method)
    }else if(constraint == FALSE){
      out <- maxLik::maxLik(logLik = ll, y = y, start = c(lam, beta, q), method = max_method)
    }

    fit <- list('MLE' = out$maximum,
                'coefficients' = list(q = out$estimate[3], beta = out$estimate[2], zero = out$estimate[1]),
                'convergence' = out$code,
                'iteration' = paste(out$iteration, 'iterations'),
                'SE' = sqrt(diag(solve(-out$hessian))),
                'hessian' = out$hessian)
    
    return(fit)
  }
}
