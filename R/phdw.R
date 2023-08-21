phdw <- function (q, q_par, beta, lam, lower.tail = TRUE, log.p = FALSE) 
{
  LLL <- max(length(q), length(q_par), length(beta), length(lam))
  if (length(q) != LLL) 
    q <- rep_len(q, LLL)
  if (length(q_par) != LLL) 
    q_par <- rep_len(q_par, LLL)
  if (length(beta) != LLL) 
    beta <- rep_len(beta, LLL)
  if (length(lam) != LLL) 
    lam <- rep_len(lam, LLL)
  
  q <- floor(q)
  ans <- lam * (pdw2(q, q = q_par, beta = beta) - ddw2(0, q = q_par, beta = beta)) / (1 - pdw2(0, q = q_par, beta = beta))
  
  ans <- ans + (1 - lam)
  ans <- ifelse(q < 0, 0, ans)

  if (lower.tail == FALSE) 
    ans <- 1 - ans
  
  if (log.p) 
    ans <- log(ans)
  if (!log.p) 
    ans <- pmin(pmax(ans, 0), 1)
  ans
}







