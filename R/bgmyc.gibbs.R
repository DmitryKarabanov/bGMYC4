bgmyc.gibbs <- function(data, m, burnin=1, thinning=1, py1, py2, pc1, pc2, t1, t2, scale=c(20, 10, 5.00), start=c(1.0, 0.5, 50.0), likelihood, prior) {
  NNodes <- data$tree$Nnode
  p <- length(start)
  vth <- array(0, dim = c(m, p+1))
  f0 <- likelihood(start, data) + prior(start, py1, py2, pc1, pc2, t1, t2)
  arate <- array(0, dim = c(1, p))
  th0 <- start
  th1 <- th0
  mover <- function(index, initial) {
    if (index == 1) return(rgamma(1, shape=scale[1], rate=(scale[1]/initial[1])))
    if (index == 2) return(rgamma(1, shape=scale[2], rate=(scale[2]/initial[2])))
    if (index == 3) return(round(initial[3] + rnorm(1) * scale[3]))
  }
  for (i in 1:m) {
    th1 <- th0
    for (j in 1:p) {
      th1[j] <- mover(j, th0)
      if (j < 3) {
        f1 <- likelihood(th1, data) + prior(th1, py1, py2, pc1, pc2, t1, t2)
        acc <- exp(f1 - f0) * (dgamma(th0[j], shape=scale[j], rate=scale[j]/th1[j]) /
                                 dgamma(th1[j], shape=scale[j], rate=scale[j]/th0[j]))
      } else {
        if (th1[3] >= 2 && th1[3] <= NNodes) {
          f1 <- likelihood(th1, data) + prior(th1, py1, py2, pc1, pc2, t1, t2)
        } else {
          f1 <- -Inf
        }
        acc <- exp(f1 - f0)
      }
      if (!is.finite(acc)) acc <- 0
      if (acc > 1) acc <- 1
      u <- runif(1) < acc
      if (u) {
        th0[j] <- th1[j]
        f0 <- f1
      }
      arate[j] <- arate[j] + u
    }
    vth[i, 1:p] <- th0
    vth[i, p+1] <- f0
    if (i %% round(m/10) == 0 || i == m) {
      message((i/m)*100, "%")
    }
  }
  arate <- arate / m
  keep_idx <- seq(burnin + thinning, m, by = thinning)
  if (length(keep_idx) == 0) keep_idx <- 1:m
  stuff <- list(
    par = vth[keep_idx, , drop = FALSE],
    accept = arate,
    tree = data$tree,
    mrca = data$mrca.nodes
  )
  message("acceptance rates py pc th: ", paste(round(stuff$accept, 3), collapse = " "))
  class(stuff) <- "singlebgmyc"
  return(stuff)
}
