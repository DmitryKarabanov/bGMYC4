bgmyc.lik <- function(params, data) {
  # Explicit cast to integer for safe list indexing
  t_idx <- as.integer(params[3])
  n <- data$n[[t_idx]]
  p <- c(rep(params[2], n), params[1])
  mat <- data$list.i.mat[[t_idx]]
  s_nod <- data$list.s.nod[[t_idx]]
  internod <- data$internod
  
  # OPT: Replace ifelse() with Boolean indexing
  log_mat <- log(mat)
  mat_p <- mat^p
  pos <- mat > 0
  mat_p[pos] <- exp(p * log_mat)[pos]
  
  # Lambda (Yule process)
  denom1 <- sum(mat_p[1:n, ] %*% internod)
  lambda1 <- sum(s_nod[1:n, ]) / denom1
  
  # Lambda (Coalescent process)
  denom2 <- sum(mat_p[n + 1, ] * internod)
  lambda2 <- sum(s_nod[n + 1, ]) / denom2
  
  # Combine into a vector
  lambda <- c(rep(lambda1, n), lambda2)
  
  b <- t(mat_p) %*% lambda
  
  lik <- b * exp(-b * internod)
  out <- sum(log(lik))
  
  # Protection against numerical artifacts
  if (!is.finite(out)) return(-Inf)
  return(out)
}
