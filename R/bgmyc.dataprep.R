#' Prepare an ultrametric tree for GMYC analysis
#'
#' Decomposes the tree into waiting times between branching events and, for
#' every possible threshold, precomputes lineage-count (`list.i.mat`) and
#' event-count (`list.s.nod`) matrices for the Yule and coalescent processes.
#'
#' @param tr An `ape` `"phylo"` object: ultrametric, fully bifurcating, no
#'   zero-length terminal branches.
#' @return A list:
#'   \item{mrca.nodes}{MRCA node indices (node number minus number of tips) per threshold.}
#'   \item{nod.types}{Node types per threshold: 0 = Yule, 1 = coalescent, 2 = species MRCA.}
#'   \item{n}{Number of coalescent processes (multi-tip species) per threshold.}
#'   \item{list.s.nod}{Event-count matrices per threshold.}
#'   \item{list.i.mat}{Lineage-count matrices per threshold (coalescent rows: n(n-1)).}
#'   \item{internod}{Waiting times between successive branching events (last: to present).}
#'   \item{tree}{The input tree.}
#' @export
bgmyc.dataprep <- function(tr) {
  if (!ape::is.ultrametric(tr)) {
    stop("Your input tree is not ultrametric. This method requires that trees be ultrametric. ",
         "If this is rounding noise from BEAST, try fix_ultrametric(tr) first.")
  }
  if (!ape::is.binary(tr)) {
    stop("Your input tree is not fully bifurcating, please resolve with zero branch lengths.")
  }
  if (0 %in% tr$edge.length[which(tr$edge[,2] <= length(tr$tip.label))]) {
    stop("Your tree contains tip branches with zero length. This will wreak havoc with the GMYC model.")
  }

  numtip <- length(tr$tip.label)
  if (numtip < 3) stop("GMYC requires at least 3 tips.")
  numnod <- tr$Nnode

  # Calculating branching times directly
  bt <- -ape::branching.times(tr)
  bt[bt > -1e-06] <- -1e-06
  names(bt) <- NULL
  sb <- sort(bt)
  nthresh <- numnod
  internod <- sb[seq(2, numnod)] - sb[seq_len(numnod - 1)]
  internod[numnod] <- 0 - sb[numnod]

  # Internal descendants of every node (ancestors are not needed by the model)
  get_nested <- function(x) {
    desc <- as.integer(tr$edge[tr$edge[, 1] == x, 2])
    if (length(desc) == 0) return(integer(0))
    res <- integer(0)
    for (d in desc) {
      if (d > numtip) res <- c(res, d, get_nested(d))
    }
    res
  }

  nested <- vector("list", numnod)
  for (i in seq_len(numnod)) {
    nested[[i]] <- get_nested(i + numtip)
  }

  # Parent/child branching times. The root has no parent: -Inf = older than any threshold.
  ancs <- cbind(tr$edge[match((1:numnod + numtip), tr$edge[, 2]), 1], (1:numnod + numtip))
  bt.ancs <- cbind(bt[ancs[, 1] - numtip], bt[ancs[, 2] - numtip])
  bt.ancs[1, 1] <- -Inf

  # Preallocation of basic data structures
  mrca.nodes <- vector("list", nthresh)
  nod.types  <- vector("list", nthresh)
  n_vec      <- vector("list", nthresh)
  list.i.mat <- vector("list", nthresh)
  list.s.nod <- vector("list", nthresh)

  ord_bt <- order(bt)

  for (j in seq(2, nthresh)) {
    threshy <- sb[j]
    tmp <- (bt.ancs[, 1] < threshy) & (bt.ancs[, 2] >= threshy)
    nod.type <- tmp + (bt >= threshy)

    mrca_idx <- which(nod.type == 2)
    mrca.nodes[[j]] <- mrca_idx
    if (nod.type[1] == 1) nod.type[1] <- 2
    nod.types[[j]] <- nod.type

    n_val <- length(mrca_idx)
    n_vec[[j]] <- n_val

    # Preallocation of fixed-size matrices
    list.s.nod[[j]] <- matrix(0, nrow = n_val + 1, ncol = numnod)
    list.i.mat[[j]] <- matrix(0, nrow = n_val + 1, ncol = numnod)

    nod_j <- nod.type[ord_bt]

    for (i in seq_len(n_val)) {
      m_idx <- mrca_idx[i]
      list.s.nod[[j]][i, m_idx] <- 2

      if (!is.null(nested[[m_idx]])) {
        nested_idx <- nested[[m_idx]] - numtip
        list.s.nod[[j]][i, nested_idx] <- 1
      }
      list.s.nod[[j]][i, ] <- list.s.nod[[j]][i, ord_bt]

      # Vectorized filling
      mask <- ifelse(list.s.nod[[j]][i, ] == 2, 2,
                     ifelse(list.s.nod[[j]][i, ] == 1, 1, 0))
      list.i.mat[[j]][i, ] <- cumsum(mask)
    }

    list.s.nod[[j]][list.s.nod[[j]] == 2] <- 1
    list.i.mat[[j]] <- list.i.mat[[j]] * (list.i.mat[[j]] - 1)

    last_row <- n_val + 1
    is_zero <- (nod_j == 0)
    is_two  <- (nod_j == 2)
    list.s.nod[[j]][last_row, ] <- is_zero
    list.i.mat[[j]][last_row, ] <- cumsum(ifelse(is_zero, 1, ifelse(is_two, -1, 0))) + 1
  }

  # Direct return of a list
  list(
    mrca.nodes = mrca.nodes,
    nod.types  = nod.types,
    n          = n_vec,
    list.s.nod = list.s.nod,
    list.i.mat = list.i.mat,
    internod   = internod,
    tree       = tr
  )
}
