bgmyc.multiphylo <- function(multiphylo, mcmc, burnin, thinning, py1=0, py2=2, pc1=0, pc2=2, t1=2, t2=51, scale=c(20, 10, 5.00), start=c(1.0, 0.5, 50.0), sampler=bgmyc.gibbs, likelihood=bgmyc.lik, prior=bgmyc.prior) {
  ntre <- length(multiphylo)
  message("You are running a multi tree analysis on ", ntre, " trees.")
  message("These trees each contain ", length(multiphylo[[1]]$tip.label), " tips.")
  message("The Yule process rate change parameter has a uniform prior ranging from ", py1, " to ", py2, ".")
  message("The coalescent process rate change parameter has a uniform prior ranging from ", pc1, " to ", pc2, ".")
  message("The threshold parameter, which is equal to the number of species, has a uniform prior ranging from ", t1, " to ", t2, ". The upper bound of this prior should not be more than the number of tips in your trees.")
  message("The MCMC will start with the Yule parameter set to ", start[1], ".")
  message("The MCMC will start with the coalescent parameter set to ", start[2], ".")
  message("The MCMC will start with the threshold parameter set to ", start[3], ". If this number is greater than the number of tips in your tree, an error will result.")
  message("Given your settings for mcmc, burnin and thinning, your analysis will result in ", ((mcmc-burnin)/thinning)*ntre, " samples being retained.")
  
  # Robust worker detection for CRAN/Check environments
  n_workers <- parallel::detectCores(logical = FALSE)
  if (is.na(n_workers) || n_workers < 1) n_workers <- 2L
  # Never spawn more workers than trees
  n_workers <- min(n_workers, ntre)
  message(sprintf("[INFO] Starting parallel analysis (using %d cores)...", n_workers))
  
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # Load packages into each worker.
  # IMPORTANT: Works only after devtools::install(), not load_all()!
  parallel::clusterEvalQ(cl, {
    library(bGMYC4)
    library(ape)
    NULL
  })
  
  # Closure captures all arguments from parent environment
  process_tree <- function(idx) {
    data <- bgmyc.dataprep(multiphylo[[idx]])
    sampler(data, m = mcmc, burnin = burnin, thinning = thinning,
            py1 = py1, py2 = py2, pc1 = pc1, pc2 = pc2, t1 = t1, t2 = t2,
            scale = scale, start = start, likelihood = likelihood, prior = prior)
  }
  
  # Clean parallel run without progress bars (they block PSOCK sockets on Windows)
  outputlist <- parallel::parLapply(cl, seq_len(ntre), process_tree)
  class(outputlist) <- "multibgmyc"
  message("[OK] All trees processed. Assembling results...")
  return(outputlist)
}
