#' Matrix of pairwise posterior probabilities of conspecificity
#'
#' Computes the marginal posterior probability that each pair of tips belongs
#' to the same species, pooling over all retained MCMC samples (and, for
#' `"multibgmyc"`, over all trees). Samples are aggregated by unique threshold
#' instead of looping over every retained sample, which is much faster on
#' large outputs.
#'
#' @param res Output of [bgmyc.singlephy()] or [bgmyc.multiphylo()].
#' @return A square matrix of class `"bgmycprobmat"` with tip-label dimnames;
#'   the diagonal is always 1.
#' @export
spec.probmat <- function(res) {
  if (inherits(res, "singlebgmyc")) {
    tr     <- res$tree
    numtip <- length(tr$tip.label)
    desc   <- .tip_descendants(tr)
    freq   <- table(factor(res$par[, 3], levels = seq_len(tr$Nnode)))
    keep   <- as.integer(names(freq[freq > 0]))
    probmat <- matrix(0, numtip, numtip,
                      dimnames = list(tr$tip.label, tr$tip.label))
    tot <- sum(freq)                       
    for (t in keep) {
      cid <- .threshold_clusters(res, t, desc, numtip)
      w   <- freq[as.character(t)] / tot
      for (u in unique(cid)) {
        idx <- which(cid == u)
        probmat[idx, idx] <- probmat[idx, idx] + w
      }
    }
  } else if (inherits(res, "multibgmyc")) {
    tr     <- res[[1]]$tree
    numtip <- length(tr$tip.label)
    probmat <- matrix(0, numtip, numtip,
                      dimnames = list(tr$tip.label, tr$tip.label))
    samp_vec <- vapply(res, function(r) nrow(r$par), integer(1))
    if (length(unique(samp_vec)) > 1L)
      warning("Trees contributed unequal numbers of samples; ",
              "pooled probabilities then weight trees by sample count. ",
              "Use identical mcmc/burnin/thinning for equal-weight model averaging (Reid & Carstens 2012).")
    for (q in seq_along(res)) {
      sq <- spec.probmat(res[[q]])
      attr(sq, "class") <- NULL
      sq <- sq[tr$tip.label, tr$tip.label, drop = FALSE]
      probmat <- probmat + sq * samp_vec[q]
    }
    probmat <- probmat / sum(samp_vec)
  } else {
    stop("spec.probmat() expects a 'singlebgmyc' or 'multibgmyc' object.")
  }
  class(probmat) <- "bgmycprobmat"
  probmat
}
