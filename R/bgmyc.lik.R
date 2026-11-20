#' GMYC log-likelihood with profiled branching rates
#'
#' @param params Numeric vector `c(py, pc, t)`.
#' @param data Output of [bgmyc.dataprep()].
#' @details Implements `L(x_i) = b_i exp(-b_i x_i)` with
#'   `b_i = lambda_y n_i^py + sum_j lambda_c (n_ij (n_ij - 1))^pc`.
#'   `lambda_y` and `lambda_c` are profile maximum-likelihood estimates
#'   (Nee 2001), shared across coalescent processes, as in Reid & Carstens
#'   (2012). Empty intervals contribute zero rate even when `p = 0`
#'   (the `0^0` guard is mathematically required).
#' @return Log-likelihood (scalar); `-Inf` for invalid states.
#' @references Reid & Carstens (2012) \doi{10.1186/1471-2148-12-196}.
#' @export
bgmyc.lik <- function(params, data) {
  t_idx <- as.integer(params[3])
  if (is.na(t_idx) || t_idx < 2 || t_idx > length(data$n)) return(-Inf)
  n     <- data$n[[t_idx]]
  mat   <- data$list.i.mat[[t_idx]]
  s_nod <- data$list.s.nod[[t_idx]]
  x     <- data$internod
  p     <- c(rep(params[2], n), params[1])

  mat_p <- mat^p                 # Direct ^: accurate and fast
  mat_p[mat == 0] <- 0           # Empty intervals do not contribute even when p=0.

  ## ML-based lambda estimates (Nee 2001)
  denom_c <- sum(mat_p[1:n, , drop = FALSE] %*% x)
  denom_y <- sum(mat_p[n + 1, ] %*% x)
  if (denom_c <= 0 || denom_y <= 0) return(-Inf)
  lambda_c <- sum(s_nod[1:n, , drop = FALSE]) / denom_c
  lambda_y <- sum(s_nod[n + 1, ]) / denom_y
  lambda   <- c(rep(lambda_c, n), lambda_y)

  b <- as.numeric(t(mat_p) %*% lambda)   # the total speed b_i in each interval
  if (any(b <= 0)) return(-Inf)

  out <- sum(log(b) - b * x)
  if (!is.finite(out)) return(-Inf)
  out
}
