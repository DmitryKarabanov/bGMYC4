bgmyc.dataprep <- function(tr) {
  if (!is.ultrametric(tr)) {
    stop("Your input tree is not ultrametric. This method requires that trees be ultrametric.")
  }
  if (!is.binary(tr)) {
    stop("Your input tree is not fully bifurcating, please resolve with zero branch lengths.")
  }
  if (0 %in% tr$edge.length[which(tr$edge[,2] <= length(tr$tip.label))]) {
    stop("Your tree contains tip branches with zero length. This will wreak havoc with the GMYC model.")
  }

  numtip <- length(tr$tip.label)
  numnod <- tr$Nnode
  numall <- numnod + numtip

  # Calculating branching times directly
  bt <- -branching.times(tr)
  bt[bt > -1e-06] <- -1e-06
  names(bt) <- NULL
  sb <- sort(bt)
  nthresh <- numnod
  internod <- sb[2:numnod] - sb[1:(numnod - 1)]
  internod[numnod] <- 0 - sb[numnod]

  # Tree traversal functions
  get_nesting <- function(x) {
    if (x < numtip + 2) return(integer(0))
    anc <- as.integer(tr$edge[tr$edge[, 2] == x, 1])
    if (length(anc) == 0) return(integer(0))
    if (anc >= numtip + 2) c(anc, get_nesting(anc)) else anc
  }

  get_nested <- function(x) {
    desc <- as.integer(tr$edge[tr$edge[, 1] == x, 2])
    if (length(desc) == 0) return(integer(0))
    res <- integer(0)
    for (d in desc) {
      if (d > numtip) res <- c(res, d, get_nested(d))
    }
    res
  }

  # Preallocating lists
  nesting <- vector("list", numnod)
  nested <- vector("list", numnod)
  for (i in seq_len(numnod)) {
    idx <- i + numtip
    nesting[[i]] <- get_nesting(idx)
    nested[[i]]  <- get_nested(idx)
  }

  ancs <- cbind(tr$edge[pmatch((1:numnod + numtip), tr$edge[, 2]), 1], (1:numnod + numtip))
  bt.ancs <- cbind(bt[ancs[, 1] - numtip], bt[ancs[, 2] - numtip])

  # Preallocation of basic data structures
  mrca.nodes <- vector("list", nthresh)
  nod.types  <- vector("list", nthresh)
  n_vec      <- vector("list", nthresh)
  list.i.mat <- vector("list", nthresh)
  list.s.nod <- vector("list", nthresh)
  nod_list   <- vector("list", nthresh)

  ord_bt <- order(bt)

  for (j in 2:nthresh) {
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
    nod_list[[j]] <- nod_j

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
