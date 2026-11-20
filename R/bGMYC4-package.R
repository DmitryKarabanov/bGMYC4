#' bGMYC4: Bayesian General Mixed Yule-Coalescent Model for Species Delimitation
#'
#' Bayesian implementation of the GMYC model (Pons et al. 2006; Reid & Carstens
#' 2012). MCMC sampling of the threshold between Yule (speciation) and coalescent
#' (within-species) branching on ultrametric trees. Pooling samples over a
#' posterior sample of trees propagates phylogenetic uncertainty into marginal
#' posterior probabilities of species limits and pairwise conspecificity.
#'
#' @section Workflow:
#' 1. [bgmyc.singlephy()] (one tree) or [bgmyc.multiphylo()] (many trees);
#' 2. diagnostics: `plot()`, `summary()`, [checkrates()], [bgmyc.chains()] +
#'    [bgmyc.convergence()];
#' 3. [spec.probmat()] -> `plot()` / [plot.interactive.probmat()];
#' 4. point estimate [bgmyc.point()]; cluster posterior [bgmyc.spec()].
#'
#' @references
#' Reid NM, Carstens BC (2012) Phylogenetic estimation error can decrease the
#' accuracy of species delimitation: a Bayesian implementation of the general
#' mixed Yule-coalescent model. *BMC Evolutionary Biology* 12:196.
#' \doi{10.1186/1471-2148-12-196}
#'
#' Pons J, Barraclough TG, Gomez-Zurita J, et al. (2006) Sequence-based species
#' delimitation for the DNA taxonomy of undescribed insects.
#' *Systematic Biology* 55(4):595--609. \doi{10.1080/10635150600852011}
#'
#' @keywords internal
"_PACKAGE"
#' @importFrom ape is.ultrametric is.binary branching.times reorder.phylo plot.phylo
#' @importFrom stats dunif dgamma rgamma rnorm runif quantile var median
#' @importFrom graphics par image plot hist points text abline
#' @importFrom grDevices pdf png dev.off colorRampPalette heat.colors
#' @importFrom utils write.table
#' @importFrom parallel makeCluster parLapply clusterEvalQ stopCluster nextRNGStream
NULL
