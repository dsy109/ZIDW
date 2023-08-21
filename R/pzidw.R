pzidw <- function (q, q_par, beta, lam, lower.tail = TRUE, log.p = FALSE) 
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
  ans <- pdw2(q, q_par, beta)
  ans <- ifelse(q < 0, 0, lam + (1 - lam) * ans)
  # deflat.limit <- -1/expm1(lambda)
  # ans[pstr0 < deflat.limit] <- NaN
  # ans[pstr0 > 1] <- NaN
  if (lower.tail == FALSE) 
    ans <- 1 - ans

  if (log.p) 
    ans <- log(ans)
  if (!log.p) 
    ans <- pmin(pmax(ans, 0), 1)
  ans
}










