pztdw <- function (q, q_par, beta, lower.tail = TRUE, log.p = FALSE) 
{
  LLL <- max(length(q), length(q_par), length(beta))
  if (length(q) != LLL) 
    q <- rep_len(q, LLL)
  if (length(q_par) != LLL) 
    q_par <- rep_len(q_par, LLL)
  if (length(beta) != LLL) 
    beta <- rep_len(beta, LLL)

  q <- floor(q)
  ans <- (pdw2(q, q_par, beta) - ddw2(0, q_par, beta)) / (1 - pdw2(0, q_par, beta))
  ans <- ifelse(q < 0, 0, ans)

  if (lower.tail == FALSE) 
    ans <- 1 - ans
  
  if (log.p) 
    ans <- log(ans)
  if (!log.p) 
    ans <- pmin(pmax(ans, 0), 1)
  ans
}




