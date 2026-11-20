#' Bayesian GMYC analysis of multiple trees (phylogenetic uncertainty)
#'
#' Runs an independent MCMC chain for each tree of a `"multiPhylo"` object
#' (e.g. a post-burn-in sample from BEAST). Pooling samples across trees yields
#' marginal posterior probabilities of species limits that integrate over
#' phylogenetic uncertainty (Reid & Carstens 2012).
#'
#' @param multiphylo An `ape` `"multiPhylo"` object; all trees must share the same
#'   tip labels and meet the requirements of [bgmyc.singlephy()].
#' @param cores Number of parallel PSOCK workers. Default 1 (sequential,
#'   CRAN-safe). Workers require the package to be installed (not just
#'   `load_all()`-ed).
#' @param seed Optional integer seed. If given, reproducible per-tree L'Ecuyer-CMRG
#'   streams are used and the caller's RNG state is restored on exit; if `NULL`,
#'   results are not reproducible.
#' @inheritParams bgmyc.singlephy
#'
#' @details Recommended protocol of Reid & Carstens (2012): about 10,000
#'   iterations per tree, burn-in 1,000, thinning 100 (about 90 samples per tree;
#'   100 trees -> 9,000 pooled samples). Precision comes from pooling across
#'   trees; Gelman-Rubin Rhat is NOT computed across trees (each tree defines a
#'   distinct conditional posterior). For replicate chains of ONE tree use
#'   [bgmyc.chains()] + [bgmyc.convergence()].
#'
#' @return Object of class `"multibgmyc"`: a list of `"singlebgmyc"` objects.
#' @references Reid & Carstens (2012) \doi{10.1186/1471-2148-12-196}.
#' @examples
#' data(east10)
#' res <- bgmyc.multiphylo(east10[1:2], mcmc = 100, burnin = 10, thinning = 10, seed = 42)
#' pm <- spec.probmat(res)
#' @export
bgmyc.multiphylo <- function(multiphylo, mcmc, burnin, thinning,
                             py1 = 0, py2 = 2, pc1 = 0, pc2 = 2,
                             t1 = 2, t2 = 51,
                             scale = c(20, 10, 5), start = c(1, 0.5, 50),
                             sampler = bgmyc.gibbs, likelihood = bgmyc.lik,
                             prior = bgmyc.prior, cores = 1L, seed = NULL) {
  # Force all arguments so run_one's closure captures plain values, not lazy
  # promises. On parallel workers the caller's frame (e.g. a script-local
  # `params`) is not available, so an unforced promise like mcmc = params$mcmc
  # would fail with "object 'params' not found".
  force(multiphylo); force(mcmc); force(burnin); force(thinning)
  force(py1); force(py2); force(pc1); force(pc2); force(t1); force(t2)
  force(scale); force(start); force(sampler); force(likelihood); force(prior)
  force(cores); force(seed)
  
  ntre <- length(multiphylo)
  if (ntre < 1L) stop("multiphylo is empty.")
  labs <- multiphylo[[1]]$tip.label
  if (!all(vapply(multiphylo, function(t) setequal(t$tip.label, labs), logical(1))))
    stop("All trees must share the same tip label set (prune or rename tips before analysis).")
  message("bGMYC4 multi-tree: ", ntre, " trees, cores = ", cores,
          if (is.null(seed)) " (seed = NULL: NOT reproducible)" else paste0(", seed = ", seed))
  
  # one independent thread of random numbers for each tree
  seeds <- NULL
  if (!is.null(seed)) {
    oldkind <- RNGkind()                       # remember the user generator
    oldseed <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({                                  # restoring the user generator
      RNGkind(oldkind[1L], oldkind[2L], oldkind[3L])
      if (!is.null(oldseed)) assign(".Random.seed", oldseed, envir = .GlobalEnv)
    }, add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    s <- .Random.seed
    seeds <- vector("list", ntre)
    for (i in seq_len(ntre)) { s <- parallel::nextRNGStream(s); seeds[[i]] <- s }
  }
  
  run_one <- function(i) {
    if (!is.null(seeds)) {
      RNGkind("L'Ecuyer-CMRG")
      assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
    }
    dat <- bgmyc.dataprep(multiphylo[[i]])
    sampler(dat, m = mcmc, burnin = burnin, thinning = thinning,
            py1 = py1, py2 = py2, pc1 = pc1, pc2 = pc2, t1 = t1, t2 = t2,
            scale = scale, start = start,
            likelihood = likelihood, prior = prior)
  }
  
  if (cores <= 1L || ntre == 1L) {
    out <- lapply(seq_len(ntre), run_one)      
  } else {
    cl <- parallel::makeCluster(min(cores, ntre))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { library(bGMYC4); library(ape); NULL })
    out <- parallel::parLapply(cl, seq_len(ntre), run_one)
  }
  class(out) <- "multibgmyc"
  out
}
