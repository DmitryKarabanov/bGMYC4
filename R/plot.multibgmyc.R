#' Trace plots of pooled multi-tree bGMYC MCMC
#'
#' @param x A `"multibgmyc"` object.
#' @param plot Logical; if `FALSE`, returns the pooled parameter matrix instead
#'   of plotting.
#' @param ... Passed to [graphics::plot()].
#' @return If `plot = TRUE`: side effects only; else a pooled matrix (`py`,
#'   `pc`, `t`, log-posterior).
#' @export
plot.multibgmyc <- function(x, plot = TRUE, ...) {
  result <- x
  # Accumulate in a list instead of rbind() in a loop
  par_list <- lapply(result, function(res) res$par)
  parmat <- do.call(rbind, par_list)

  if (plot) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    par(mfrow = c(2, 2))
    ylabels <- c("p.div", "p.coal", "threshold", "logposterior")
    for (i in 1:4) {
      vals <- parmat[, i]
      # Protection from NA/Inf
      if (any(is.finite(vals))) {
        plot(vals, xlab = "generations", ylab = ylabels[i])
      } else {
        plot(0, 0, type = "n", xlab = "generations", ylab = ylabels[i],
             main = "No finite values", xaxt = "n", yaxt = "n")
        text(0, 0, labels = "Chains unstable", col = "red")
      }
    }
  } else {
    return(parmat)
  }
}
