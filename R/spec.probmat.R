spec.probmat <- function(res) {
  # Optimized internal function for getting list of clusters
  spec.list <- function(res, thresh) {
    tr       <- res$tree
    spec     <- tr$tip.label
    numtip   <- length(tr$tip.label)
    max.mrca <- res$mrca[[thresh]] + numtip
    numspec  <- length(max.mrca)
    
    # Recursive collection of tip indexes
    nest.tip <- function(nod, tr) {
      child <- tr$edge[tr$edge[, 1] == nod, 2]
      tips  <- integer(0)
      for (ch in child) {
        if (ch <= numtip) {
          tips <- c(tips, ch)
        } else {
          tips <- c(tips, nest.tip(ch, tr))
        }
      }
      tips
    }
    
    # Generating a list of clusters
    res_list <- vector("list", length(max.mrca))
    for (i in seq_along(max.mrca)) {
      tip.name <- tr$tip.label[nest.tip(max.mrca[i], tr)]
      res_list[[i]] <- data.frame(spec = i, name = tip.name, stringsAsFactors = FALSE)
    }
    
    assigned_names <- unlist(lapply(res_list, `[[`, 2))
    missing <- spec[is.na(match(spec, assigned_names))]
    
    if (length(missing) > 0) {
      
      for (s in missing) {
        numspec <- numspec + 1
        res_list[[length(res_list) + 1]] <- data.frame(spec = numspec, name = s, stringsAsFactors = FALSE)
      }
    }
    
    res2 <- do.call(rbind, res_list)
    colnames(res2) <- c("GMYC_spec", "sample_name")
    res2$GMYC_spec  <- as.factor(res2$GMYC_spec)
    res2$sample_name <- as.factor(res2$sample_name)
    return(res2)
  }
  
  if (inherits(res, "singlebgmyc")) {
    numtip    <- length(res$tree$tip.label)
    probmat   <- matrix(0, nrow = numtip, ncol = numtip)
    rownames(probmat) <- colnames(probmat) <- res$tree$tip.label
    n_samples <- nrow(res$par)
    add_val   <- 1 / n_samples
    
    for (j in seq_len(n_samples)) {
      assignlists <- spec.list(res, res$par[j, 3])
      
      clusters <- split(as.character(assignlists$sample_name), assignlists$GMYC_spec)
      for (cl in clusters) {
        probmat[cl, cl] <- probmat[cl, cl] + add_val
      }
    }
  }
  
  if (inherits(res, "multibgmyc")) {
    numtip    <- length(res[[1]]$tree$tip.label)
    probmat   <- matrix(0, nrow = numtip, ncol = numtip)
    rownames(probmat) <- colnames(probmat) <- res[[1]]$tree$tip.label
    ntrees    <- length(res)
    totalsamp <- ntrees * nrow(res[[1]]$par)
    
    for (q in seq_len(ntrees)) {
      for (j in seq_len(nrow(res[[q]]$par))) {
        assignlists <- spec.list(res[[q]], res[[q]]$par[j, 3])
        clusters <- split(as.character(assignlists$sample_name), assignlists$GMYC_spec)
        for (cl in clusters) {
          probmat[cl, cl] <- probmat[cl, cl] + 1
        }
      }
    }
    probmat <- probmat / totalsamp
  }
  
  class(probmat) <- "bgmycprobmat"
  return(probmat)
}
