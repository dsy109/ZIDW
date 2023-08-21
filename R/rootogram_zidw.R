rootogram_zidw <- function (object, type = c("hanging", "standing", "suspended"),
                            sqrt = TRUE, ref_line = TRUE, warn_limits = TRUE, fitted_colour = "steelblue",
                            bar_colour = NA, bar_fill = "grey", ref_line_colour = "black",
                            warn_line_colour = "black", ylab = NULL, xlab = NULL, ...)
  
{
  
  type <- match.arg(type)
  
  mf <- model.frame(object)
  y <- model.response(mf)
  
  max_count <- max(y)
  bin <- seq(0, max_count, 1)
  names(bin) <- as.character(bin)
  obs <- table(factor(y, levels = bin))
  
  mu <- predict(object, newdata = mf, type = 'response')
  
  q_par <- inv.logit(object$model_matrix_q %*% cbind(object$coefficients$q))
  beta <- exp(object$model_matrix_beta %*% cbind(object$coefficients$beta))
  lam <- inv.logit(object$model_matrix_zi %*% cbind(object$coefficients$zero))
  
  fitted <- purrr::map_dfc(bin, .f = dzidw, q_par = q_par, beta = beta, lam = lam)
  
  fitted <- colSums(fitted)
  
  df <- tibble::tibble(bin = bin, observed = as.integer(obs[seq_along(bin)]), fitted = fitted)
  
  if (is.null(ylab)) {
    
    ylab <- if (as.logical(sqrt)) {
      
      df <- mutate(df, observed = sqrt(.data$observed),
                       
                       fitted = sqrt(.data$fitted))
      
      expression(sqrt(Frequency))
      
    }
    
    else {
      
      "Frequency"
      
    }
    
  }
  
  if (is.null(xlab)) {
    
    xlab <- "Response"
    
  }
  
  #  distr <- attr(object, "distribution")
  
  width <- 0.45
  
  #  width <- if (distr %in% c("Poisson", "Negative Binomial")) {
  
  #    0.45
  
  #  }
  
  #  else {
  
  #    w <- attr(object, "width")
  
  #    (w * 0.9)/2
  
  #  }
  
  df <- mutate(df, x_low = .data$bin - width, x_high = .data$bin +
                     
                     width)
  
  nr <- nrow(df)
  
  df <- if (type == "hanging") {
    
    mutate(df, y_bot = .data$fitted - .data$observed,
           
           y_top = .data$y_bot + .data$observed)
    
  }
  
  else if (type == "suspended") {
    
    mutate(df, y_bot = rep(0, nr), y_top = .data$fitted -
             
             .data$observed)
    
  }
  
  else {
    
    mutate(df, y_bot = rep(0, nr), y_top = .data$y_bot +
             
             .data$observed)
    
  }
  
  plt <- ggplot(df) + geom_rect(aes(xmin = .data$x_low,
                                        
                                        xmax = .data$x_high, ymin = .data$y_bot, ymax = .data$y_top),
                                    
                                    fill = bar_fill, col = bar_colour) + geom_line(aes(x = .data$bin,
                                                                                       
                                                                                       y = .data$fitted), colour = fitted_colour, linewidth = 1) +
    
    geom_point(aes(x = .data$bin, y = .data$fitted), colour = fitted_colour,
               
               size = 2.5)
  
  if (as.logical(ref_line) & type != "standing") {
    
    plt <- plt + geom_hline(yintercept = 0, colour = ref_line_colour)
    
  }
  
  if (as.logical(warn_limits) & type != "standing") {
    
    plt <- plt + geom_hline(yintercept = c(-1, 1), lty = "dashed",
                            
                            colour = warn_line_colour)
    
  }
  
  plt <- plt + labs(y = ylab, x = xlab)
  
  plt
  
}
