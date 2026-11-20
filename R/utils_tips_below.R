# Internal helpers for fast post-processing (spec.probmat, bgmyc.spec).
# The species clustering is a deterministic function of (tree, threshold):
# tip descendants are computed once per tree, then clusters are built
# once per unique threshold and weighted by its frequency.

#' @noRd
.tip_descendants <- function(tr) {
  numtip <- length(tr$tip.label)
  nnode  <- tr$Nnode
  child  <- split(tr$edge[, 2], tr$edge[, 1])
  desc   <- vector("list", nnode)
  bt     <- -ape::branching.times(tr)      # more negative = older
  ord    <- order(bt, decreasing = TRUE)   # youngest nodes first
  for (v in ord) {
    cd    <- child[[as.character(numtip + v)]]
    tips  <- cd[cd <= numtip]
    nodes <- cd[cd >  numtip]
    desc[[v]] <- c(tips, unlist(desc[nodes - numtip]))
  }
  desc
}

#' @noRd
.threshold_clusters <- function(single, t, desc, numtip) {
  mrc <- single$mrca[[t]]
  cid <- integer(numtip)
  for (i in seq_along(mrc)) cid[desc[[mrc[i]]]] <- i
  singles <- which(cid == 0L)
  if (length(singles)) cid[singles] <- seq_along(singles) + length(mrc)
  cid
}
