print.zidw <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 
                                                        0.85)), "", sep = "\n")
  # if (x$convergence != 0) {
  #   cat("model did not converge\n")
  # }

    cat(paste("q model coefficients with logit link:\n", 
              sep = ""))
    print.default(format(x$coefficients$q, digits = digits), 
                  print.gap = 2, quote = FALSE)

    cat(paste("\nZero-inflation model coefficients (binomial with logit link):\n", sep = ""))
    print.default(format(x$coefficients$zero, digits = digits), 
                  print.gap = 2, quote = FALSE)
    
    cat("\n")
    
    cat(paste("Beta model coefficients with log link:\n", 
              sep = ""))
    print.default(format(x$coefficients$beta, digits = digits), 
                  print.gap = 2, quote = FALSE)
    
    cat("\n")
  
  invisible(x)
}

