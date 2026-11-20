#' Point estimate of species limits from a conspecificity matrix
#'
#' Greedy single-linkage clustering: a sequence joins a cluster if its
#' probability of conspecificity with ANY member exceeds `ppcutoff`; clusters
#' with overlapping memberships are merged.
#'
#' @param probmat A `"bgmycprobmat"` matrix from [spec.probmat()].
#' @param ppcutoff Posterior probability threshold for lumping tips into
#'   species. [spec.probmat()] returns a marginal probability pooled over
#'   thresholds and trees, so it is not comparable to a clade posterior
#'   probability. Because clustering is single linkage, a species is held
#'   together by its weakest link: raising the cutoff removes bridges between
#'   intraspecific lineages and fragments real species, whereas pairs from
#'   different species usually have probabilities near zero. Empirical data
#'   typically show a stability plateau at low values (0.001 to 0.05); 0.05 is
#'   the value used in the original bGMYC documentation. Choose the cutoff from
#'   a sweep over several values and report the plateau in your methods.
#' @return A list of character vectors, each giving the members of one species.
#'   Because of phylogenetic uncertainty, not every pair within a returned
#'   species necessarily exceeds `ppcutoff`; members are connected through
#'   intermediate individuals (single linkage). A warning flags such chains.
#' @examples
#' set.seed(123)
#' data(east10)
#' res <- bgmyc.singlephy(east10[[1]], mcmc = 100, burnin = 10, thinning = 10)
#' bgmyc.point(spec.probmat(res), ppcutoff = 0.5)
#' @export
bgmyc.point <- function(probmat, ppcutoff){
  seqs <- rownames(probmat)
  if (length(seqs) < 2L) return(list(seqs))   # guard: single-taxon matrix
  pointest <- list(seqs[1])
  nspec <- 1
  for(i in 2:length(seqs)){
    flip <- 0
    for (j in seq_along(pointest)) {
      if(any(probmat[i, pointest[[j]]] > ppcutoff)){
        pointest[[j]] <- c(pointest[[j]], seqs[i])
        flip <- flip + 1
      }
    }
    if(flip == 0){
      nspec <- nspec + 1
      pointest[[nspec]] <- seqs[i]
    }
    if(flip > 1){
      warning("Not all members of each cluster have conspecificity probability
              greater than specified with every other member in the cluster.")
    }
  }
  while(any(duplicated(unlist(pointest)))){
    message("Merging clusters due to overlapping memberships.")
    dupl <- duplicated(unlist(pointest))
    dupes <- unlist(pointest)[dupl]
    lvec <- c()
    for (i in seq_along(pointest)) { lvec <- c(lvec, dupes[1] %in% pointest[[i]]) }
    newspec <- unlist(pointest[lvec])
    newspec <- newspec[!duplicated(newspec)]
    pointest <- pointest[!lvec]
    pointest[[length(pointest) + 1]] <- newspec
  }
  return(pointest)
}
