#' Multiple independent MCMC chains for ONE tree
#'
#' Runs `nchains` independent chains on the same tree (optionally from
#' overdispersed starts drawn from the priors) for Gelman-Rubin diagnostics.
#'
#' @param phylo An `ape` `"phylo"` object (see [bgmyc.singlephy()]).
#' @param nchains Number of chains (>= 2).
#' @param seed Optional integer; per-chain L'Ecuyer-CMRG streams, caller RNG restored.
#' @param cores Number of parallel PSOCK workers (default 1).
#' @param dispersed Logical; draw overdispersed starting values from the priors.
#' @inheritParams bgmyc.singlephy
#' @return Object of class `"bgmycchains"`: `chains` (list of `"singlebgmyc"`),
#'   `nchains`, `dispersed`.
#' @seealso [bgmyc.convergence()], [bgmyc.pool()]
#' @examples
#' data(east10)
#' ch <- bgmyc.chains(east10[[1]], mcmc = 100, burnin = 10, thinning = 10,
#'                    nchains = 2, seed = 123)
#' bgmyc.convergence(ch)
#' @export
bgmyc.chains <- function(phylo, mcmc, burnin, thinning,
                         py1 = 0, py2 = 2, pc1 = 0, pc2 = 2,
                         t1 = 2, t2 = NULL,
                         scale = c(20, 10, 5), start = NULL,
                         sampler = bgmyc.gibbs, likelihood = bgmyc.lik,
                         prior = bgmyc.prior,
                         nchains = 4L, seed = NULL, cores = 1L,
                         dispersed = TRUE) {
  # Force all arguments so the worker closure captures plain values, not lazy
  # promises. On parallel workers the caller's frame (e.g. a script-local
  # `params`) is not available, so an unforced promise like mcmc = params$mcmc
  # would fail with "object 'params' not found".
  force(phylo); force(mcmc); force(burnin); force(thinning)
  force(py1); force(py2); force(pc1); force(pc2); force(t1); force(t2)
  force(scale); force(start); force(sampler); force(likelihood); force(prior)
  force(nchains); force(seed); force(cores); force(dispersed)

  nchains <- as.integer(nchains)
  if (nchains < 2L) stop("nchains must be >= 2 (Gelman-Rubin needs >= 2 chains).")
  if (is.null(seed))
    message("NOTE: seed = NULL -> chains are not reproducible.")

  dat <- bgmyc.dataprep(phylo)
  NNodes <- dat$tree$Nnode
  if (is.null(t2))    t2 <- NNodes
  if (is.null(start)) start <- c(max(1e-3, (py1 + py2) / 2),
                                 max(1e-3, (pc1 + pc2) / 2),
                                 max(2, (max(2, t1) + min(t2, NNodes)) / 2))
  th_lo <- max(2, t1); th_hi <- min(t2, NNodes)
  if (th_lo > th_hi) stop("Empty threshold range after clamping to [2, Nnode].")

  ## per-chain RNG streams (same scheme as bgmyc.multiphylo)
  seeds <- NULL
  if (!is.null(seed)) {
    oldkind <- RNGkind()
    oldseed <- if (exists(".Random.seed", envir = .GlobalEnv))
                 get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      RNGkind(oldkind[1L], oldkind[2L], oldkind[3L])
      if (!is.null(oldseed)) assign(".Random.seed", oldseed, envir = .GlobalEnv)
    }, add = TRUE)
    RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
    s <- .Random.seed
    seeds <- vector("list", nchains)
    for (i in seq_len(nchains)) { s <- parallel::nextRNGStream(s); seeds[[i]] <- s }
  }

  run_chain <- function(i) {
    if (!is.null(seeds)) {
      RNGkind("L'Ecuyer-CMRG")
      assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
    }
    st <- start
    if (dispersed) {                      # overdispersed starts from the priors
      st[1] <- max(runif(1, py1, py2), 1e-3)
      st[2] <- max(runif(1, pc1, pc2), 1e-3)
      st[3] <- th_lo + sample.int(th_hi - th_lo + 1L, 1L) - 1L
    }
    sampler(dat, m = mcmc, burnin = burnin, thinning = thinning,
            py1 = py1, py2 = py2, pc1 = pc1, pc2 = pc2, t1 = t1, t2 = t2,
            scale = scale, start = st, likelihood = likelihood, prior = prior)
  }

  message("bGMYC4 chains: ", nchains, " x ", mcmc, " steps, cores = ", cores,
          if (dispersed) ", dispersed starts" else ", common start")
  if (cores <= 1L) {
    runs <- lapply(seq_len(nchains), run_chain)
  } else {
    cl <- parallel::makeCluster(min(cores, nchains))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { library(bGMYC4); library(ape); NULL })
    runs <- parallel::parLapply(cl, seq_len(nchains), run_chain)
  }
  structure(list(chains = runs, nchains = nchains, dispersed = dispersed),
            class = "bgmycchains")
}

# split-Rhat: n x m matrix (retained samples x chains)
#' @noRd
.split_rhat <- function(x) {
  n <- nrow(x); m <- ncol(x)
  if (n < 4L || m < 2L) return(NA_real_)
  half <- floor(n / 2)
  sp <- matrix(NA_real_, half, 2L * m)
  for (j in seq_len(m)) {
    sp[, 2L*j - 1L] <- x[1:half, j]
    sp[, 2L*j]      <- x[(n - half + 1L):n, j]
  }
  nn <- half; mm <- 2L * m
  cm <- colMeans(sp)
  B  <- nn * sum((cm - mean(cm))^2) / (mm - 1L)
  W  <- mean(apply(sp, 2L, var))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt(((nn - 1L)/nn * W + B/nn) / W)
}

#' Convergence diagnostics for replicate chains of one tree
#'
#' Split-Rhat (target < 1.05) and per-chain ESS (target > 200; requires
#' `mcmcse`). Rhat is computed ONLY across chains of the same tree, never
#' across trees.
#' @param chains A `"bgmycchains"` object from [bgmyc.chains()].
#' @return List: `summary` (data frame with Rhat, ESS, acceptance, flags),
#'   `per_chain_accept`, `per_chain_ess`, `rhat`, `ess_total`.
#' @export
bgmyc.convergence <- function(chains) {
  if (!inherits(chains, "bgmycchains")) stop("Expected an object from bgmyc.chains().")
  runs <- chains$chains; m <- length(runs)
  n <- min(vapply(runs, function(r) nrow(r$par), integer(1)))
  mats <- lapply(runs, function(r) r$par[seq_len(n), , drop = FALSE])
  pnames <- c("py", "pc", "th", "logpost")

  rhat <- vapply(1:4, function(k) {
    X <- do.call(cbind, lapply(mats, function(z) z[, k]))
    .split_rhat(X)
  }, numeric(1))
  names(rhat) <- pnames

  ess <- matrix(NA_real_, m, 4L,
                dimnames = list(paste0("chain", seq_len(m)), pnames))
  if (requireNamespace("mcmcse", quietly = TRUE)) {
    for (j in seq_len(m)) ess[j, ] <- apply(mats[[j]], 2L, mcmcse::ess)
  } else message("Install 'mcmcse' for ESS: install.packages('mcmcse')")

  acc <- do.call(rbind, lapply(runs, function(r) as.numeric(r$accept)))
  colnames(acc) <- c("py", "pc", "th")
  ess_total <- colSums(ess); ess_min <- apply(ess, 2L, min)

  summary <- data.frame(
    parameter     = pnames,
    Rhat          = round(rhat, 4),
    ESS_total     = round(ess_total, 1),
    ESS_min_chain = round(ess_min, 1),
    accept_mean   = c(round(colMeans(acc), 3), NA),
    flag = ifelse(!is.na(rhat) & rhat > 1.05, "Rhat>1.05: run longer / tune scale",
           ifelse(!is.na(ess_total) & ess_total < 200, "ESS<200: run longer", "OK")),
    stringsAsFactors = FALSE
  )
  list(summary = summary, per_chain_accept = acc, per_chain_ess = ess,
       rhat = rhat, ess_total = ess_total)
}

#' Pool replicate chains of one tree into a single sample
#'
#' Concatenates retained samples of chains targeting the SAME conditional
#' posterior into one `"singlebgmyc"` object.
#' @param chains A `"bgmycchains"` object.
#' @return A `"singlebgmyc"` object (see [bgmyc.singlephy()]).
#' @export
bgmyc.pool <- function(chains) {
  if (!inherits(chains, "bgmycchains")) stop("Expected an object from bgmyc.chains().")
  runs <- chains$chains
  out <- list(
    par    = do.call(rbind, lapply(runs, function(r) r$par)),
    accept = colMeans(do.call(rbind, lapply(runs, function(r) as.numeric(r$accept)))),
    tree   = runs[[1]]$tree,
    mrca   = runs[[1]]$mrca,          # <- запятая обязательна
    priors = runs[[1]]$priors         # <- новое поле
  )
  class(out) <- "singlebgmyc"
  out
}
