

zidw_reg <- function(qformula, betaformula = ~ 1, ziformula = ~ 1, data, lam = NULL, beta = NULL, q = NULL, k = 1000, uni_method =  c('MLE', 'mde'), max_method = NULL, constraint = TRUE, B = NULL){
  uni_method <- match.arg(uni_method)
  
  n <- nrow(data)
  
  outq <- lm(qformula,data=data)
  mf <- outq$model
  X <- model.matrix(outq)
  # ind <- which(colnames(data)==rownames(attr(outq$terms,"factors"))[1])
  ind <- which(colnames(data)==names(attr(outq$terms,"dataClasses"))[1])
  Y <- data[,ind]
  

    
  old.Y <- colnames(data)[ind] 
  colnames(data)[ind] <- "Y"
  new.form <- as.formula(paste("Y",paste(format(terms(betaformula)), collapse = '')))
  Z <- model.matrix(lm(new.form,data=data))
  
  new.form.zi<- as.formula(paste("Y",paste(format(terms(ziformula)), collapse = '')))
  W <- model.matrix(lm(new.form.zi,data=data))
    
  
  data_beta <- data.frame(Z, Y)
  data_q <- data.frame(X, Y)
  data_lam <- data.frame(W, Y = as.numeric(Y == 0))
  para.q1 <- ncol(X) > 1
  para.beta <- ncol(Z) > 1
  para.lam <- ncol(W) > 1

  
  kx <- ncol(X)
  kz <- ncol(Z)
  kw <- ncol(W)

  
  # par=c(lam, beta, q)
  ll <- function(par, y, X, W, Z){
    delta <- ifelse(y == 0, 1, 0)
    
    # lam <- par[1:ncol(W)]
    
    kx <- ncol(X)
    kz <- ncol(Z)
    kw <- ncol(W)
    p <- sum(kx, kz, kw)

    
    q <- c(inv.logit(X %*% cbind(par[-c(1:(kw + kz))])))
    beta <- c(exp(Z %*% cbind(par[-c(1:kw, (kw + kz + 1):p)])))
    lam <- c(inv.logit(W %*% cbind(par[-c((kw + 1) : p)])))
    
    loglike <- sum(delta * log(lam + (1 - lam) * (1 - q)) + (1 - delta) * log((1 - lam) * (exp(y^beta * log(q)) - exp((y + 1)^beta * log(q)))  ))
    # loglike <- sum(delta * log(lam + (1 - lam) * (1 - q)) + (1 - delta) * ll_tran)
    return(loglike)
    
    
    
  }
  
  
  

  ui <- function(p){
    
    temp <- rbind(c(1,0,rep(0,p)),c(-1,0,rep(0,p)),c(0,1,rep(0,p)),c(0,-1,rep(0,p)))
    
    for(i in 1:p){
      
      temp <- rbind(temp,c(0,0,rep(0,i-1),1,rep(0,p-i)),c(0,0,rep(0,i-1),-1,rep(0,p-i)))
      
    }
    
    temp
    
  }
  
  ## constraint for q
  if(para.q1){
    ui_kx <- ui(kx)
    ci_kx <- rep(c(-50, -50), kx)
  }else{
    ui_kx <- ui(kx)
    ci_kx <- c(0, -1)
  }

  
  ## constraint for beta
  if(para.beta){
    ui_kz <- ui(kz)
    ci_kz <- rep(c(-50, -50), kz)
  }else{
    ui_kz <- ui(kz)
    ci_kz <- c(0, -50)
  }

  ## constraint for lambda
  if(para.lam){
    ui_kw <- ui(kw)
    ci_kw <- rep(c(-50, -50), kw)
  }else{
    ui_kw <- ui(kw)
    ci_kw <- c(0, -1)
  }

  ## combine constraint
  # fuzz <- - 1e-6
  ui_all <- direct.sum(direct.sum(ui_kw, ui_kz), ui_kx)
  ci_all <- c(ci_kw, ci_kz, ci_kx)

  ####### initial values for q
  data1 <- data_q[which(data_q$Y > 0), ]

  
  if(is.null(q)){
    if(para.q1){
      out.q <- suppressWarnings(dw.reg(I(Y - 1) ~ . -1, data = data1, para.q1 = TRUE))
      
      q0 <- as.numeric(coef(out.q)[1:kx])
    }else{
      q0 <- dw.parest2(data1$Y - 1)$q
    }
  }else{
    q0 <- q
  }

  
  ### initial values for lambda

  
  if(is.null(lam)){
    if(para.lam){
      lam0 <- as.numeric(coef(glm(Y ~ . -1, data = data_lam, family = 'binomial')))
    }else{
      lam0 <- mean(Y == 0)
    }
  }else{
    lam0 <- lam
  }

  
  
  ##  beta
  data2 <- data_beta[which(data_beta$Y > 0), ]

  
  if(is.null(beta)){
    if(para.beta){
      out.beta <- suppressWarnings(dw.reg(I(Y - 1) ~ . -1, data = data2, para.beta = TRUE))
      
      beta0 <- as.numeric(coef(out.beta)[1:kz])
    }else{
      beta0 <- dw.parest2(data2$Y - 1)$beta
    }
  }else{
    beta0 <- beta
  }

  d <- length(c(lam0, beta0, q0))
  ui_all <- ui_all[1:(2*d), 1:d]



  X <- as.matrix(X)


  

  Z <- as.matrix(Z)


  W <- as.matrix(W)

 
  if(!para.beta && !para.lam && !para.q1){
    if(uni_method == 'MLE'){
      out <- zidw_uni(Y, method = uni_method, max_method = max_method, constraint = constraint, par = c(lam0, beta0, q0))
      coefq <- gtools::logit(out$coefficients$q)
      coefb <- log(out$coefficients$beta)
      coefzi <- gtools::logit(out$coefficients$zero)
      names(coefq) <- names(coefb) <- names(coefzi) <- '(Intercept)'
      names(out)[1] <- 'value'
      SE <- out$SE
    }else{
      if(!is.null(B)){
        out <- zidw_uni(Y, method = uni_method, par = c(lam0, beta0, q0), B = B)
        coefq <- gtools::logit(out$out$estimate[3])
        coefb <- log(out$out$estimate[2])
        coefzi <- gtools::logit(out$out$estimate[1])
        names(coefq) <- names(coefb) <- names(coefzi) <- '(Intercept)'
        out[3] <- 'value'
        names(out)[3] <- 'value'
        SE <- out$SE
        
        out$value <- sum(dzidw(Y, q_par = out$out$estimate[3], beta = out$out$estimate[2], lam = out$out$estimate[1], log = TRUE))
      }else{
        out <- zidw_uni(Y, method = uni_method, par = c(lam0, beta0, q0))
        coefq <- gtools::logit(out$estimate[3])
        coefb <- log(out$estimate[2])
        coefzi <- gtools::logit(out$estimate[1])
        names(coefq) <- names(coefb) <- names(coefzi) <- '(Intercept)'
        out[3] <- 'value'
        names(out)[3] <- 'value'
        SE <- NULL
        out$value <- sum(dzidw(Y, q_par = out$estimate[3], beta = out$estimate[2], lam = out$estimate[1], log = TRUE))
      }
      
    }
    
    
  }else{
    out <- constrOptim(theta = c(lam0, beta0, q0), f = ll, y = Y, X = X, Z = Z, W = W, ui = ui_all, ci = ci_all, grad = NULL, control=list(fnscale=-1))
    
    ## SE
    ineqA = ui_all
    ineqB = -ci_all
    out1 <- maxLik::maxLik(logLik = ll, y = Y, X = X, Z = Z, W = W, start = out$par, constraints = list(ineqA = ineqA, ineqB = ineqB), method = 'NM')
    
    
    SE <- sqrt(abs(diag(solve(-out1$hessian))))

    ## coefficients and covariances
    coefq <- out$par[(length(out$par) - kx + 1): length(out$par)]
    names(coefq) <- colnames(X)
    coefzi <- out$par[1 : kw]
    names(coefzi) <- colnames(W)
    coefb <- out$par[(kw + 1) : (kw + kz)]
    names(coefb) <- colnames(Z)
    
  }
    cl <- match.call()
    
    ## fitted and residuals
    q <- c(inv.logit(X %*% cbind(coefq)))
    beta <- c(exp(Z %*% cbind(coefb)))
    lam <- c(inv.logit(W %*% cbind(coefzi)))
    # dim 801
  
    success <- FALSE
    while(!success){
      y <- seq(1, k, 1)
      mu <- sapply(1:length(q), function(i) sum(q[i]^y^beta[i]))
      Yhat <- (1 - lam) * mu
      res <- Y - Yhat
      success <- max(round(q^k^beta, 8)) == 0
      k <- k * 10
    }
    
    
    
    # print(out$convergence)
    fit <- list(call = cl,
                coefficients = list(zero = coefzi, beta = coefb, q = coefq),
                loglik = out$value,
                SE = SE,
                convergence = out$convergence,
                nall = n,
                res = res,
                fitted_values = Yhat,
                response = Y,
                model_matrix_q = X,
                model_matrix_beta = Z,
                model_matrix_zi = W,
                model = mf,
                formula = list(qformula = qformula, betaformula = betaformula, ziformula = ziformula))
    class(fit) <- 'zidw'
    # return(fit)
    fit

}











