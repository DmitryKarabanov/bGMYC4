#' Trace plots of a single-tree bGMYC MCMC
#'
#' 2x2 panel: `py`, `pc`, threshold, log-posterior.
#' @param x A `"singlebgmyc"` object.
#' @param burnin,thinning Applied to the ALREADY thinned retained samples;
#'   leave at defaults (0, 1). A warning is issued otherwise.
#' @param ... Passed to [graphics::plot()].
#' @return No return value; called for side effects.
#' @export
plot.singlebgmyc <- function(x, burnin = 0, thinning = 1, ...) {
  if (burnin != 0 || thinning != 1)
    warning("'x$par' already contains only retained (post burn-in, thinned) samples; ",
            "'burnin'/'thinning' here re-filter the retained samples, not the original chain.")
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  result <- x
  par(mfrow = c(2, 2))
  ylabels <- c("p.div", "p.coal", "threshold", "logposterior")
  for (i in 1:4) {
    # Burnin/Thinning-Aware Index Extraction
    idx <- seq(burnin + 1, nrow(result$par), by = thinning)
    vals <- result$par[idx, i]
    # Check for NA/Inf
    finite_mask <- is.finite(vals)
    if (any(finite_mask)) {
      plot(vals, xlab = "generations", ylab = ylabels[i], ...)
      # Highlight invalid points in gray
      if (!all(finite_mask)) {
        points(which(!finite_mask), vals[!finite_mask], col = "gray", pch = 20)
      }
    } else {
      plot(0, 0, type = "n", xlab = "generations", ylab = ylabels[i],
           main = "No finite values", xaxt = "n", yaxt = "n")
      text(0, 0, labels = "Chain unstable / NA or Inf", col = "red", cex = 1.2)
    }
  }
}
