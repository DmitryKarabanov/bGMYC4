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
  get_desc <- function(v) {
    if (!is.null(desc[[v]])) return(desc[[v]])
    cd <- child[[as.character(numtip + v)]]
    if (is.null(cd)) cd <- integer(0)
    res <- c(cd[cd <= numtip],
             unlist(lapply(cd[cd > numtip] - numtip, get_desc)))
    desc[[v]] <<- res
    res
  }
  for (v in seq_len(nnode)) get_desc(v)
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
