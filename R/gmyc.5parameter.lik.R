#' GMYC likelihood with all five parameters free (experimental)
#'
#' @param params Numeric vector `c(py, pc, lambda_y, lambda_c, t)` — rate-change
#'   parameters first, then branching rates, then the threshold index.
#' @param data Output of [bgmyc.dataprep()].
#' @details Not used by the default sampler. The Pons et al. implementation and
#'   [bgmyc.lik()] profile the branching rates at their ML estimates (3 free
#'   parameters); letting all 5 vary can produce a flat likelihood surface with
#'   a peak far from the true value. Reid & Carstens (2012) found the profile
#'   version empirically better; this function is kept for experimentation.
#' @return Log-likelihood (scalar).
#' @export
gmyc.5parameter.lik <- function(params, data) {
  t_idx <- as.integer(params[5])
  n     <- data$n[[t_idx]]
  p     <- c(rep(params[2], n), params[1])
  lambda <- c(rep(params[4], n), params[3])
  mat_p <- data$list.i.mat[[t_idx]]^p
  mat_p[data$list.i.mat[[t_idx]] == 0] <- 0   # empty intervals contribute zero rate, even at p = 0
  b <- as.numeric(t(mat_p) %*% lambda)
  if (any(b <= 0)) return(-Inf)
  sum(log(b) - b * data$internod)             # stable form of log(b * exp(-b*x))
}
