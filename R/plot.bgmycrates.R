#' @export
plot.bgmycrates <- function(x, ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  checkratesoutput <- x
  par(mfrow = c(2, 2))
  plot(log(checkratesoutput[, 5] / checkratesoutput[, 6], base = 10), ylim = c(-1, 3),
       xlab = "mcmc samples", ylab = "log(coalescence rate/Yule rate)")
  plot((checkratesoutput[, 5] / checkratesoutput[, 6]), ylim = c(-10, 100),
       xlab = "mcmc samples", ylab = "coalescence rate/Yule rate")
  hist(log(checkratesoutput[, 5] / checkratesoutput[, 6], base = 10), breaks = 20,
       xlab = "log(coalescence rate/Yule rate)", main = NULL)
  hist((checkratesoutput[, 5] / checkratesoutput[, 6]), breaks = 20,
       xlab = "coalescence rate/Yule rate", main = NULL)
}
