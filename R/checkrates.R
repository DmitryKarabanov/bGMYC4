checkrates <- function(result) {
  bgmyc.lambda <- function(params, data) {
    n <- data$n[[params[3]]]
    p <- c(rep(params[2], n), params[1])
    lambda <- sum(data$list.s.nod[[params[3]]][1:n, ]) /
      sum(data$list.i.mat[[params[3]]][1:n, ]^p[-(n + 1)] %% data$internod)
    lambda <- c(lambda, sum(data$list.s.nod[[params[3]]][n + 1, ]) /
                  (data$list.i.mat[[params[3]]][n + 1, ]^p[n + 1] %% data$internod))
    return(lambda[c(2, 1)])
  }
  
  # Accumulate rates in a list instead of rbind() in a loop
  rates_list <- list()
  idx <- 1
  if (inherits(result, "multibgmyc")) {
    for (i in seq_along(result)) {
      message("Processing tree ", i, " of ", length(result))
      dat <- bgmyc.dataprep(result[[i]]$tree)
      n_rows <- nrow(result[[i]]$par)
      for (j in seq_len(n_rows)) {
        rates_list[[idx]] <- c(result[[i]]$par[j, 1:3],
                               bgmyc.lambda(result[[i]]$par[j, 1:3], dat))
        idx <- idx + 1
      }
    }
  }
  if (inherits(result, "singlebgmyc")) {
    dat <- bgmyc.dataprep(result$tree)
    n_rows <- nrow(result$par)
    for (j in seq_len(n_rows)) {
      rates_list[[idx]] <- c(result$par[j, 1:3],
                             bgmyc.lambda(result$par[j, 1:3], dat))
      idx <- idx + 1
    }
  }
  branchrates <- do.call(rbind, rates_list)
  ratemod <- ((branchrates[, 3]^branchrates[, 1])) / branchrates[, 3]
  branchrates <- cbind(branchrates, branchrates[, 4] * ratemod)
  colnames(branchrates) <- c("p.div", "p.coal", "threshold", "lambda.div",
                             "lambda.coal", "lambda.div.mod")
  class(branchrates) <- "bgmycrates"
  return(branchrates)
}
