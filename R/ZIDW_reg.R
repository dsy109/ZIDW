

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
   out1 <- out 
    
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

    vc <- tryCatch(abs(-solve(as.matrix(out1$hessian))),
                   error = function(e){
                     warning(e$message, call = FALSE)
                     k <- nrow(as.matrix(out1$hessian))
                     return(matrix(NA, k, k))
                   })

    colnames(vc) <- rownames(vc) <- c(paste("q", colnames(X), sep = "_"),
                                      paste('beta', colnames(Z), sep = '_'),
                                      paste("zero",  colnames(W), sep = "_"))
    
    
    
    # print(out$convergence)
    fit <- list(call = cl,
                coefficients = list(zero = coefzi, beta = coefb, q = coefq),
                loglik = out$value,
                SE = SE,
                convergence = out$convergence,
                nall = n,
                residuals = res,
                df.residual = n - (kx + kz + kw),
                fitted_values = Yhat,
                response = Y,
                vcov = vc,
                model_matrix_q = X,
                model_matrix_beta = Z,
                model_matrix_zi = W,
                model = mf,
                formula = list(qformula = qformula, betaformula = betaformula, ziformula = ziformula))
    class(fit) <- 'zidw'
    # return(fit)
    fit

}





summary.zidw <- function(object,...)
{
  ## residuals
  object$residuals <- residuals(object, type = "pearson")
  
  ## compute z statistics
  kb <- length(object$coefficients$beta)
  kq <- length(object$coefficients$q)
  kz <- length(object$coefficients$zero)
  
  se <- sqrt(diag(object$vcov))
  coef <- c(object$coefficients$q, object$coefficients$beta, object$coefficients$zero)  
  # if(object$dist == "negbin") {
  #   coef <- c(coef[1:kc], "Log(theta)" = log(object$theta), coef[(kc+1):(kc+kz)])
  #   se <- c(se[1:kc], object$SE.logtheta, se[(kc+1):(kc+kz)])
  #   kc <- kc+1
  # }
  zstat <- coef / se
  pval <- 2 * pnorm(-abs(zstat))
  coef <- cbind(coef, se, zstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients$q <- coef[1:kq,,drop = FALSE]
  object$coefficients$beta <- coef[(kq + 1) : (kq + kb),, drop = FALSE]
  object$coefficients$zero <- coef[(kq + kb + 1) : (kq + kb + kz),,drop = FALSE]
  
  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL
  
  ## return
  class(object) <- "summary.zidw"
  object
}

print.summary.zidw <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!as.logical(x$convergence)) {
    cat("model did not converge\n")
  } else {
    
    cat("Pearson residuals:\n")
    print(structure(quantile(x$residuals),
                    names = c("Min", "1Q", "Median", "3Q", "Max")), digits = digits, ...)  
    
    cat(paste("\nq model coefficients (with logit link):\n", sep = ""))
    printCoefmat(x$coefficients$q, digits = digits, signif.legend = FALSE)
    
    cat(paste("\nBeta model coefficients (with log link):\n", sep = ""))
    printCoefmat(x$coefficients$beta, digits = digits, signif.legend = FALSE)
    
    cat(paste("\nZero-inflation model coefficients (binomial with logit link):\n", sep = ""))
    printCoefmat(x$coefficients$zero, digits = digits, signif.legend = FALSE)
    
    if(getOption("show.signif.stars") & any(rbind(x$coefficients$q, x$coefficients$beta, x$coefficients$zero)[,4] < 0.1, na.rm=TRUE))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    
    # if(x$dist == "negbin") cat(paste("\nTheta =", round(x$theta, digits), "\n")) else cat("\n")
    # cat(paste("Number of iterations in", x$method, "optimization:", tail(na.omit(x$optim$count), 1), "\n"))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits), "on", x$nall - x$df.residual, "Df\n")
  }
  
  invisible(x)
}






