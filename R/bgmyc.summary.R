#' Shortest discrete credible interval (HPD analogue)
#'
#' Shortest contiguous interval of integer values carrying at least `prob`
#' posterior mass. For multimodal posteriors a true highest-density set may be
#' disjoint; this is the standard shortest-interval analogue used for the
#' threshold (number of species).
#' @param x Numeric/integer vector of posterior samples (e.g. threshold column).
#' @param prob Coverage probability (default 0.95).
#' @return Integer vector `c(lo, hi)`.
#' @export
hpd_discrete <- function(x, prob = 0.95) {
  tb <- table(as.integer(x))
  v  <- as.integer(names(tb))
  p  <- as.numeric(tb) / sum(tb)
  cs <- cumsum(p); n <- length(v)
  best <- c(v[1L], v[n]); bestw <- Inf; j <- 1L
  for (i in seq_len(n)) {
    base <- if (i > 1L) cs[i - 1L] else 0
    while (j <= n && cs[j] - base < prob) j <- j + 1L
    if (j > n) break
    if (v[j] - v[i] < bestw) { bestw <- v[j] - v[i]; best <- c(v[i], v[j]) }
  }
  best
}

#' Posterior summaries for bGMYC analyses
#'
#' Threshold mean/median, 95% shortest credible interval and equal-tailed
#' interval; parameter summaries; ESS (via `mcmcse`); prior-vs-posterior
#' contraction. For `"multibgmyc"`, thresholds are pooled across trees
#' (marginal posterior); per-tree diagnostics are reported separately and
#' Rhat is NOT computed across trees.
#' @param object A `"singlebgmyc"`, `"bgmycchains"` or `"multibgmyc"` object.
#' @param prob Coverage probability (default 0.95).
#' @param priors Optional named vector `py1, py2, pc1, pc2, t1, t2` overriding
#'   the priors stored in the object.
#' @param ... Unused.
#' @return Object of class `"bgmycsummary"` (print method provided).
#' @name bgmyc_summary
NULL

# Contraction of posterior 95% width vs uniform prior 95% width + edge mass.
#' @noRd
.prior_check <- function(post, lo, hi) {
  pw <- 0.95 * (hi - lo)
  qw <- diff(unname(quantile(post, c(0.025, 0.975))))
  edge <- mean(post > hi - 0.025*(hi - lo) | post < lo + 0.025*(hi - lo))
  list(contraction = unname(round(1 - qw/pw, 3)),
       edge_mass = round(edge, 3),
       informative = unname((1 - qw/pw) > 0.2))
}

#' @rdname bgmyc_summary
#' @export
summary.singlebgmyc <- function(object, prob = 0.95, priors = NULL, ...) {
  par <- object$par; th <- par[, 3]
  pri <- if (!is.null(priors)) priors else object$priors
  hpd <- hpd_discrete(th, prob)
  qs  <- unname(quantile(th, c(0.025, 0.975)))
  params <- data.frame(
    parameter = c("py", "pc"),
    mean = colMeans(par[, 1:2, drop = FALSE]),
    q025 = unname(apply(par[, 1:2, drop = FALSE], 2, quantile, probs = 0.025)),
    q975 = unname(apply(par[, 1:2, drop = FALSE], 2, quantile, probs = 0.975)),
    row.names = NULL)
  ess <- if (requireNamespace("mcmcse", quietly = TRUE))
    apply(par, 2, function(z) min(mcmcse::ess(z), length(z))) else rep(NA_real_, 4)
  prior <- if (!is.null(pri)) list(
    py = .prior_check(par[,1], pri[1], pri[2]),
    pc = .prior_check(par[,2], pri[3], pri[4]),
    th = .prior_check(th,      pri[5], pri[6])) else NULL
  out <- list(type = "single tree", n_samples = nrow(par),
              threshold = c(mean = mean(th), median = median(th),
                            hpd_lo = hpd[1], hpd_hi = hpd[2],
                            q025 = qs[1], q975 = qs[2]),
              params = params, accept = object$accept, ESS = ess, prior_check = prior)
  class(out) <- "bgmycsummary"; out
}

#' @rdname bgmyc_summary
#' @export
summary.bgmycchains <- function(object, prob = 0.95, priors = NULL, ...) {
  pooled <- bgmyc.pool(object)
  s <- summary(pooled, prob = prob, priors = priors)
  conv <- bgmyc.convergence(object)
  s$type <- paste0(object$nchains, " chains, one tree")
  s$Rhat <- conv$rhat; s$ESS <- conv$ess_total
  s
}

#' @rdname bgmyc_summary
#' @export
summary.multibgmyc <- function(object, prob = 0.95, ...) {
  nt <- length(object)
  th_all <- unlist(lapply(object, function(r) r$par[, 3]))
  per_tree <- data.frame(
    tree = seq_len(nt),
    n_samples = vapply(object, function(r) nrow(r$par), integer(1)),
    th_mean = vapply(object, function(r) mean(r$par[,3]), numeric(1)),
    hpd_lo = vapply(object, function(r) as.numeric(hpd_discrete(r$par[,3], prob)[1]), numeric(1)),
    hpd_hi = vapply(object, function(r) as.numeric(hpd_discrete(r$par[,3], prob)[2]), numeric(1)),
    acc_py = vapply(object, function(r) r$accept[1], numeric(1)),
    acc_pc = vapply(object, function(r) r$accept[2], numeric(1)),
    acc_th = vapply(object, function(r) r$accept[3], numeric(1)))
  if (requireNamespace("mcmcse", quietly = TRUE)) {
    ess1 <- function(x) min(mcmcse::ess(x), length(x))   # ESS <= n по определению
    per_tree$ESS_py <- vapply(object, function(r) ess1(r$par[,1]), numeric(1))
    per_tree$ESS_pc <- vapply(object, function(r) ess1(r$par[,2]), numeric(1))
    per_tree$ESS_th <- vapply(object, function(r) ess1(r$par[,3]), numeric(1))
  }
  hpd <- hpd_discrete(th_all, prob)
  qs  <- unname(quantile(th_all, c(0.025, 0.975)))
  out <- list(type = paste0(nt, " trees (pooled marginal)"), n_samples = length(th_all),
              threshold = c(mean = mean(th_all), median = median(th_all),
                            hpd_lo = hpd[1], hpd_hi = hpd[2],
                            q025 = qs[1], q975 = qs[2]),
              per_tree = per_tree,
              note = "Rhat is NOT computed across trees (each tree defines a distinct conditional posterior).")
  class(out) <- "bgmycsummary"; out
}

#' @rdname bgmyc_summary
#' @param x A `"bgmycsummary"` object returned by `summary()`.
#' @export
print.bgmycsummary <- function(x, ...) {
  cat("bGMYC posterior summary:", x$type, "|", x$n_samples, "samples\n")
  cat("Threshold (species): mean =", round(x$threshold["mean"], 2),
      "| 95% HPD [", x$threshold["hpd_lo"], ",", x$threshold["hpd_hi"], "]",
      "| equal-tailed [", x$threshold["q025"], ",", x$threshold["q975"], "]\n")
  if (!is.null(x$params)) print(x$params, row.names = FALSE)
  if (!is.null(x$Rhat)) { cat("Split-Rhat:\n"); print(round(x$Rhat, 4)) }
  if (!is.null(x$ESS))  { cat("ESS:\n"); print(round(x$ESS, 1)) }
  if (!is.null(x$prior_check)) {
    cat("Prior contraction (1 - post/prior 95% width):\n")
    print(sapply(x$prior_check, function(z) z$contraction))
  }
  if (!is.null(x$per_tree)) { cat("Per-tree diagnostics:\n"); print(x$per_tree, row.names = FALSE) }
  if (!is.null(x$note)) cat(x$note, "\n")
  invisible(x)
}
