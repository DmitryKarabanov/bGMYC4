bgmyc.point <- function(probmat, ppcutoff){
  seqs <- rownames(probmat)
  pointest <- list(seqs[1])
  nspec <- 1
  for(i in 2:length(seqs)){
    flip <- 0
    for(j in 1:length(pointest)){
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
      warning("Not all members of each cluster have conspecificity probability greater than specified with every other member in the cluster.")
    }
  }
  pointfin <- list()
  nspec <- 1
  for(i in 1:(((length(seqs)^2) - length(seqs))) / 2){
  }
  while(any(duplicated(unlist(pointest)))){
    message("Merging clusters due to overlapping memberships.")
    dupl <- duplicated(unlist(pointest))
    dupes <- unlist(pointest)[dupl]
    lvec <- c()
    for(i in 1:length(pointest)){ lvec <- c(lvec, dupes[1] %in% pointest[[i]]) }
    newspec <- unlist(pointest[lvec])
    newspec <- newspec[!duplicated(newspec)]
    pointest <- pointest[!lvec]
    pointest[[length(pointest) + 1]] <- newspec
  }
  return(pointest)
}
