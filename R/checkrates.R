#' Branching rates and model-fit diagnostics
#'
#' Recomputes the Yule and coalescent branching rates for every retained MCMC
#' sample. GMYC assumes the two processes act on very different time scales
#' (roughly an order of magnitude apart); if the rates are similar, the
#' threshold is indistinct and results unreliable.
#'
#' @param result A `"singlebgmyc"` or `"multibgmyc"` object.
#' @return Object of class `"bgmycrates"`: matrix with columns `p.div`, `p.coal`,
#'   `threshold`, `lambda.div`, `lambda.coal`, `lambda.div.mod`. The last is
#'   `lambda.div * t^p.div / t`, i.e. the fold change in speciation rate from
#'   the root to the last speciation event (Reid & Carstens 2012).
#' @examples
#' set.seed(123)
#' data(east10)
#' res <- bgmyc.singlephy(east10[[1]], mcmc = 100, burnin = 10, thinning = 10)
#' rates <- checkrates(res)
#' plot(rates)
#' @export
checkrates <- function(result) {

  bgmyc.lambda <- function(params, data) {
    t_idx  <- as.integer(params[3])
    n      <- data$n[[t_idx]]
    p      <- c(rep(params[2], n), params[1])
    lam_coal <- sum(data$list.s.nod[[t_idx]][1:n, ]) /
      sum(data$list.i.mat[[t_idx]][1:n, ]^p[-(n + 1)] %*% data$internod)
    lam_div <- sum(data$list.s.nod[[t_idx]][n + 1, ]) /
      as.numeric(data$list.i.mat[[t_idx]][n + 1, ]^p[n + 1] %*% data$internod)
    return(c(lam_div, lam_coal))  # [1]=lambda.div, [2]=lambda.coal
  }

  rates_list <- list()
  idx <- 1
  if (inherits(result, "multibgmyc")) {
    for (i in seq_along(result)) {
      message("Processing tree ", i, " of ", length(result))
      dat <- bgmyc.dataprep(result[[i]]$tree)
      for (j in seq_len(nrow(result[[i]]$par))) {
        rates_list[[idx]] <- c(result[[i]]$par[j, 1:3],
                               bgmyc.lambda(result[[i]]$par[j, 1:3], dat))
        idx <- idx + 1
      }
    }
  }
  if (inherits(result, "singlebgmyc")) {
    dat <- bgmyc.dataprep(result$tree)
    for (j in seq_len(nrow(result$par))) {
      rates_list[[idx]] <- c(result$par[j, 1:3],
                             bgmyc.lambda(result$par[j, 1:3], dat))
      idx <- idx + 1
    }
  }

  branchrates <- do.call(rbind, rates_list)
  ratemod <- ((branchrates[, 3]^branchrates[, 1])) / branchrates[, 3]
  branchrates <- cbind(branchrates, branchrates[, 4] * ratemod)
  colnames(branchrates) <- c("p.div", "p.coal", "threshold",
                             "lambda.div", "lambda.coal", "lambda.div.mod")
  class(branchrates) <- "bgmycrates"
  return(branchrates)
}
