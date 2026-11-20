#' MCMC sampler for the Bayesian GMYC model
#'
#' Metropolis-within-Gibbs sampler used by [bgmyc.singlephy()],
#' [bgmyc.multiphylo()] and [bgmyc.chains()]; rarely needed directly.
#'
#' @param data Output of [bgmyc.dataprep()].
#' @param m Total number of MCMC iterations.
#' @inheritParams bgmyc.singlephy
#'
#' @details `py` and `pc` use mean-preserving gamma proposals with the correct
#'   Hastings correction; `t` uses a symmetric rounded-normal random walk
#'   (Hastings ratio 1), rejected outside `[2, Nnode]`. Retained samples are
#'   `seq(burnin + thinning, m, by = thinning)`.
#'
#' @return An object of class `"singlebgmyc"` (see [bgmyc.singlephy()]).
#' @examples
#' data(east10)
#' dat <- bgmyc.dataprep(east10[[1]])
#' out <- bgmyc.gibbs(dat, m = 100, burnin = 10, thinning = 10,
#'                    py1 = 0, py2 = 2, pc1 = 0, pc2 = 2, t1 = 2, t2 = 100,
#'                    start = c(1, 0.5, 50),
#'                    likelihood = bgmyc.lik, prior = bgmyc.prior)
#' @export
bgmyc.gibbs <- function(data, m, burnin = 1, thinning = 1,
                        py1, py2, pc1, pc2, t1, t2,
                        scale = c(20, 10, 5), start = c(1, 0.5, 50),
                        likelihood, prior) {
  NNodes <- data$tree$Nnode
  p <- length(start)

  ## Input check: it does not change the results in correct runs.
  if (m < 1 || burnin < 0 || burnin >= m || thinning < 1)
    stop("Need mcmc >= 1, 0 <= burnin < mcmc, thinning >= 1.")
  if (length(scale) != 3 || any(scale <= 0))
    stop("'scale' must be three positive numbers: proposal widths for py, pc and t.")
  if (!(py1 < py2 && pc1 < pc2 && t1 < t2))
    stop("Prior bounds must satisfy py1 < py2, pc1 < pc2 and t1 < t2.")

  if (any(start[1:2] <= 0))
    stop("start[1] and start[2] must be > 0 (gamma proposals).")
  if (start[3] < 2 || start[3] > NNodes)
    stop("start[3] must be in [2, Nnode].")
  if (start[1] < py1 || start[1] > py2 ||
      start[2] < pc1 || start[2] > pc2 ||
      start[3] < t1  || start[3] > t2)
    stop("start values must lie inside the priors.")
  if (py2 > 2)
    warning("'py2' > 2: Reid & Carstens (2012) consider such Yule rate-change values unrealistic; ",
            "they recommend U(0, 1.2) (about 2-fold speciation rate increase) for empirical data.")
  if (pc2 > 2)
    warning("'pc2' > 2: coalescent rate-change values above 2 are not meaningfully interpretable demographically.")
  if (t1 < 2)
    warning("'t1' < 2: the threshold (number of species) cannot be below 2; the effective prior is truncated to [", max(2, t1), ", ", min(t2, NNodes), "].")
  if (t2 > NNodes)
    warning("'t2' > Nnode (", NNodes, "): thresholds above Nnode are unreachable; the effective prior is truncated to [", max(2, t1), ", ", NNodes, "].")
  if (m < 100)
    warning("'mcmc' < 100: too few iterations for reliable MCMC inference; Reid & Carstens (2012) recommend >= 10 000 per tree.")
  if (burnin / m > 0.9)
    warning(sprintf("'burnin' is %.0f%% of 'mcmc': almost the entire chain will be discarded.", 100 * burnin / m))

  vth   <- matrix(0, nrow = m, ncol = p + 1)
  th0   <- start
  f0    <- likelihood(th0, data) + prior(th0, py1, py2, pc1, pc2, t1, t2)
  if (!is.finite(f0))
    stop("Initial state has zero posterior density; check start, priors and tree.")
  arate <- numeric(p)

  mover <- function(j, th) {
    if (j == 1) return(rgamma(1, shape = scale[1], rate = scale[1]/th[1]))
    if (j == 2) return(rgamma(1, shape = scale[2], rate = scale[2]/th[2]))
    round(th[3] + rnorm(1) * scale[3])   # симметричное предложение, Hastings = 1
  }

  step <- max(1L, round(m/10))           # protection against %% 0 in short circuits
  for (i in 1:m) {
    for (j in 1:p) {
      th1 <- th0
      th1[j] <- mover(j, th0)

      if (j < 3) {
        f1 <- likelihood(th1, data) + prior(th1, py1, py2, pc1, pc2, t1, t2)
        acc <- exp(f1 - f0) *
          (dgamma(th0[j], shape = scale[j], rate = scale[j]/th1[j]) /
             dgamma(th1[j], shape = scale[j], rate = scale[j]/th0[j]))
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
      if (u) { th0 <- th1; f0 <- f1 }
      arate[j] <- arate[j] + u
    }
    vth[i, 1:p] <- th0
    vth[i, p+1] <- f0                    # Always log-posterior(th0)
    if (i %% step == 0 || i == m) message(round(100*i/m), "%")
  }

  keep <- seq(burnin + thinning, m, by = thinning)
  if (length(keep) == 0)
    stop("No samples retained: check mcmc/burnin/thinning.")

  out <- list(par = vth[keep, , drop = FALSE], accept = arate/m,
              tree = data$tree, mrca = data$mrca.nodes,
              priors = c(py1 = py1, py2 = py2, pc1 = pc1, pc2 = pc2, t1 = t1, t2 = t2))
  class(out) <- "singlebgmyc"
  message("acceptance rates (py, pc, th): ",
          paste(round(out$accept, 3), collapse = " "))
  out
}
