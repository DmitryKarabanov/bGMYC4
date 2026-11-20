#' Ten posterior trees from a simulated dataset
#'
#' A `"multiPhylo"` object of 10 ultrametric trees (215 tips, 214 internal
#' nodes each) sampled from the BEAST posterior of a simulated single-locus
#' dataset; used in examples and the vignette.
#' @format A `"multiPhylo"` list of 10 `"phylo"` objects with shared tip labels.
#' @source Simulated example data inherited from the original bGMYC package
#'   (Reid & Carstens 2012); tree estimation as in that study.
#' @examples
#' data(east10)
#' length(east10)
"east10"
