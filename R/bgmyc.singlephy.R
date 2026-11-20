#' Bayesian GMYC analysis of a single ultrametric tree
#'
#' Runs an MCMC sampler for the general mixed Yule-coalescent model on one
#' ultrametric, fully bifurcating tree with unique haplotypes only.
#'
#' @param phylo An `ape` `"phylo"` object. Must be ultrametric, fully resolved,
#'   and without zero-length terminal branches. Prune identical haplotypes
#'   BEFORE tree estimation (Reid & Carstens 2012); see [fix_ultrametric()] for
#'   rounding-level non-ultrametricity.
#' @param mcmc Total number of MCMC iterations (not the number of retained samples).
#' @param burnin Number of initial iterations discarded as burn-in.
#' @param thinning Interval at which samples are retained; the number of retained
#'   samples is `floor((mcmc - burnin) / thinning)`.
#' @param py1,py2 Bounds of the uniform prior on the Yule rate change parameter
#'   `py`. The speciation rate with `n` lineages is proportional to `n^py`
#'   (Pons et al. 2006). With 50 species, `py = 2` allows a 50-fold rate increase
#'   from root to tips, blurring the Yule/coalescent transition; `py2 <= 1.5`
#'   (about 7-fold) is more realistic, and Reid & Carstens (2012) also used
#'   U(0, 1.2). This is the most influential of the three priors.
#' @param pc1,pc2 Bounds of the uniform prior on the coalescent rate change
#'   parameter `pc`; the within-species rate is proportional to `(n(n-1))^pc`.
#'   In principle `pc < 1` mimics population decline and `pc > 1` growth, but
#'   because identical haplotypes must be excluded, a strict demographic
#'   interpretation is not possible.
#' @param t1,t2 Bounds of the uniform prior on the threshold parameter `t`,
#'   which equals the number of delimited species (clusters + singleton tips).
#'   Use `2 <= t1 <= t2 <= Nnode(tree)`; larger values are unreachable and merely
#'   truncate the effective prior.
#' @param scale Numeric(3): proposal widths for `py`, `pc` (mean-preserving gamma
#'   proposals; larger = more conservative) and `t` (rounded normal random walk;
#'   larger = more liberal). Tune via acceptance rates (target roughly 0.2--0.4).
#' @param start Numeric(3): starting values `c(py, pc, t)`; must lie inside the
#'   priors, with `start[3]` in `[2, Nnode(tree)]`.
#' @param sampler MCMC sampler function; default [bgmyc.gibbs()].
#' @param likelihood Log-likelihood function; default [bgmyc.lik()] (profile ML
#'   estimates of branching rates, Nee 2001).
#' @param prior Log-prior function; default [bgmyc.prior()] (uniform). Replace
#'   `dunif` with `dgamma`/`dlnorm` for other prior families.
#'
#' @details Branching rates lambda are NOT sampled: at every step they are fixed
#'   at their maximum-likelihood estimates conditional on the current threshold
#'   and rate-change parameters, exactly as in Reid & Carstens (2012); rates and
#'   rate-change parameters are shared across all coalescent processes. The
#'   resulting posterior is therefore an empirical-Bayes (profile) posterior.
#'   Check stationarity with `plot(output)` and `summary(output)`; acceptance
#'   rates are computed over all iterations, including burn-in.
#'
#' @return An object of class `"singlebgmyc"`:
#'   \item{par}{Matrix of retained samples: `py`, `pc`, `t`, log-posterior.}
#'   \item{accept}{Acceptance rates for the three parameters.}
#'   \item{tree}{The input tree.}
#'   \item{mrca}{MRCA node indices per threshold (from [bgmyc.dataprep()]).}
#'   \item{priors}{Named vector `py1, py2, pc1, pc2, t1, t2`.}
#'
#' @references Reid & Carstens (2012) \doi{10.1186/1471-2148-12-196};
#'   Pons et al. (2006) \doi{10.1080/10635150600852011}.
#' @seealso [bgmyc.multiphylo()], [bgmyc.chains()], [spec.probmat()], [checkrates()]
#' @examples
#' data(east10)
#' res <- bgmyc.singlephy(east10[[1]], mcmc = 100, burnin = 10, thinning = 10)
#' plot(res)
#' @export
bgmyc.singlephy <-
function(phylo, mcmc, burnin, thinning, py1=0, py2=2, pc1=0, pc2=2, t1=2, t2=51, scale=c(20, 10, 5.00), start=c(1.0, 0.5, 50.0), sampler=bgmyc.gibbs, likelihood=bgmyc.lik, prior=bgmyc.prior){
	
  message("You are running bGMYC on a single phylogenetic tree.\n")
  message("This tree contains ", length(phylo$tip.label), " tips.\n")
  message("The Yule process rate change parameter has a uniform prior ranging from ", py1, " to ", py2, ".\n")
  message("The coalescent process rate change parameter has a uniform prior ranging from ", pc1, " to ", pc2, ".\n")
  message("The threshold parameter, which is equal to the number of species, has a uniform prior ranging from ", t1, " to ", t2, ". The upper bound of this prior should not be more than the number of tips in your trees.\n")
  message("The MCMC will start with the Yule parameter set to ", start[1], ".\n")
  message("The MCMC will start with the coalescent parameter set to ", start[2], ".\n")
  message("The MCMC will start with the threshold parameter set to ", start[3], ". If this number is greater than the number of tips in your tree, an error will result.\n")
  message("Given your settings for mcmc, burnin and thinning, your analysis will result in ", (mcmc-burnin)/thinning, " samples being retained.\n")
  
	bgmyc.dataprep(phylo)->data
	sampler(data, m=mcmc, burnin, thinning, py1, py2, pc1, pc2, t1, t2, scale=scale, start=start, likelihood=likelihood, prior=prior)->output
	return(output)
	
	}
