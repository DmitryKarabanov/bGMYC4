#' Fix rounding-level non-ultrametricity
#'
#' Rounds branch lengths to `digits` decimals; stops if the tree remains
#' non-ultrametric (re-export from BEAST or check the clock model instead).
#' Preprocessing only — the GMYC likelihood is unaffected.
#' @param tr A `"phylo"` object.
#' @param digits Number of decimals for rounding (default 8).
#' @return The (possibly rounded) tree.
#' @export
fix_ultrametric <- function(tr, digits = 8) {
  if (!inherits(tr, "phylo")) stop("tr must be a 'phylo' object.")
  tr <- ape::reorder.phylo(tr, "cladewise")
  n <- length(tr$tip.label)
  
  ## 1. Round to clean BEAST export noise.
  tr$edge.length <- round(tr$edge.length, digits)
  
  ## 2. Root-to-tip depths and path edge counts (cladewise: parents before children).
  nn <- n + tr$Nnode
  depth <- numeric(nn)
  nedge <- integer(nn)
  for (k in seq_len(nrow(tr$edge))) {
    a <- tr$edge[k, 1]; b <- tr$edge[k, 2]
    depth[b] <- depth[a] + tr$edge.length[k]
    nedge[b] <- nedge[a] + 1L
  }
  tipdepth <- depth[seq_len(n)]
  spread <- max(tipdepth) - min(tipdepth)
  budget <- 10^-digits * max(nedge[seq_len(n)])
  
  ## 3. Reject genuinely non-clock trees.
  if (spread > budget)
    stop("Tree is genuinely non-clock: root-to-tip spread ", format(spread),
         " exceeds the rounding budget ", format(budget),
         ". Re-export from BEAST or check the clock model.")
  
  ## 4. Align all tips to the deepest one by extending ONLY terminal branches.
  target <- max(tipdepth)
  term <- which(tr$edge[, 2] <= n)
  for (k in term) {
    b <- tr$edge[k, 2]
    tr$edge.length[k] <- tr$edge.length[k] + (target - depth[b])
  }
  tr
}
