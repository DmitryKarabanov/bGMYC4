#' Posterior distribution of species clusters
#'
#' Processes the output of [bgmyc.singlephy()] or [bgmyc.multiphylo()] and
#' tabulates every distinct species partition with its posterior probability.
#' Optionally returns the posterior collection of species-by-locality matrices
#' for downstream community analyses (Picante, Vegan). Cluster strings are
#' generated once per unique threshold and weighted by its frequency.
#'
#' @param res A `"singlebgmyc"` or `"multibgmyc"` object.
#' @param filename Optional path; the cluster table is also written to this file.
#' @param cmatrix Optional tips-by-localities matrix (rownames = tip labels).
#' @return A list with `specprobs` and, optionally, `matrices`.
#' @export
bgmyc.spec <- function(res, filename = NULL, cmatrix = NULL) {
    ## validate cmatrix BEFORE agg_tree() uses it
    if (!is.null(cmatrix)) {
      tips <- if (inherits(res, "singlebgmyc")) res$tree$tip.label
      else if (inherits(res, "multibgmyc")) res[[1]]$tree$tip.label
      else stop("bgmyc.spec() expects a 'singlebgmyc' or 'multibgmyc' object.")
      if (is.null(rownames(cmatrix)) || !all(tips %in% rownames(cmatrix)))
        stop("cmatrix must have rownames covering all tree tip labels.")
    }
  agg_tree <- function(single) {
    tr     <- single$tree
    numtip <- length(tr$tip.label)
    desc   <- .tip_descendants(tr)
    freq   <- table(factor(single$par[, 3], levels = seq_len(tr$Nnode)))
    keep   <- as.integer(names(freq[freq > 0]))
    strings <- character(0)
    counts  <- integer(0)
    cl   <- list()      # per unique threshold: list of sorted tip-label vectors
    mats <- list()      # per unique threshold: species-by-locality matrix
    for (t in keep) {
      cid <- .threshold_clusters(single, t, desc, numtip)
      uu  <- sort(unique(cid))
      st  <- character(length(uu))
      clt <- vector("list", length(uu))
      for (k in seq_along(uu)) {
        labs <- sort(tr$tip.label[which(cid == uu[k])])
        clt[[k]] <- labs
        st[k] <- if (length(labs) > 1L) paste(labs, collapse = "\t") else labs
      }
      strings <- c(strings, st)
      counts  <- c(counts, rep(as.integer(freq[as.character(t)]), length(uu)))
      cl[[as.character(t)]] <- clt
      if (!is.null(cmatrix)) {
        rows <- vapply(clt, function(l) l[1], character(1))
        m <- matrix(NA_real_, nrow = length(uu), ncol = ncol(cmatrix),
                    dimnames = list(rows, colnames(cmatrix)))
        for (k in seq_along(uu)) {
          m[k, ] <- colSums(cmatrix[clt[[k]], , drop = FALSE])
        }
        mats[[as.character(t)]] <- m
      }
    }
    list(strings = strings, counts = counts, n_samples = sum(freq),
         clusters = cl, matrices = mats, thresholds = keep)
  }

  if (inherits(res, "singlebgmyc")) {
    a       <- agg_tree(res)
    per_t   <- a
    strings <- a$strings
    counts  <- a$counts
    tot     <- a$n_samples
  } else if (inherits(res, "multibgmyc")) {
    samp_vec <- vapply(res, function(r) nrow(r$par), integer(1))
    if (length(unique(samp_vec)) > 1L)
      warning("Trees contributed unequal numbers of samples; ",
              "pooled probabilities then weight trees by sample count. ",
              "Use identical mcmc/burnin/thinning for equal-weight model averaging (Reid & Carstens 2012).")
    parts   <- lapply(res, agg_tree)
    strings <- unlist(lapply(parts, `[[`, "strings"))
    counts  <- unlist(lapply(parts, `[[`, "counts"))
    tot     <- sum(vapply(parts, `[[`, integer(1), "n_samples"))
    per_t   <- parts
  } else {
    stop("bgmyc.spec() expects a 'singlebgmyc' or 'multibgmyc' object.")
  }

  # specprobs: identical table to the per-sample implementation
  tab <- tapply(counts, strings, sum)
  ord <- order(as.numeric(tab), decreasing = TRUE)
  specprobs <- data.frame(count       = as.integer(tab[ord]),
                          probability = as.numeric(tab[ord]) / tot,
                          cluster     = names(tab)[ord],
                          stringsAsFactors = FALSE)
  if (!is.null(filename))
    write.table(specprobs, filename, row.names = FALSE, col.names = FALSE)

  output <- list(specprobs = specprobs)
  if (!is.null(cmatrix)) {
    ## one matrix per retained sample, taken from the per-threshold cache
    if (inherits(res, "singlebgmyc")) {
      thr <- res$par[, 3]
      output$matrices <- lapply(thr, function(t) per_t$matrices[[as.character(t)]])
    } else {
      output$matrices <- unlist(
        lapply(seq_along(res), function(q)
          lapply(res[[q]]$par[, 3], function(t)
            per_t[[q]]$matrices[[as.character(t)]])),
        recursive = FALSE)
    }
  }
  output
}
