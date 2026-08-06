\name{bGMYC}
\alias{bGMYC}
\docType{package}
\title{A Bayesian MCMC implementation of the general mixed Yule-coalescent model for species delimitation.}
\description{
This package implements a Bayesian version of the general mixed Yule-coalescent
model for species delimitation introduced by Pons et al. 2006. It uses
Markov Chain Monte Carlo simulation to evaluate the posterior distribution
of the model and uses a flexible method of specifying prior distributions on the model parameters. 
The primary strength of this method over the pre-existing implementation is that it allows the use of multiple trees sampled from an MCMC analysis, 
eliminating the conditioning of results on a single point estimate of the phylogenetic tree, which often is associated with high uncertainty.
}
\details{
\tabular{ll}{
Package: \tab bGMYC4\cr
Type: \tab Package\cr
Version: \tab 4.1.0\cr
Date: \tab 2026-05-01\cr
License: \tab GPL (>= 3)\cr
LazyData: \tab true\cr
}
}
\author{
Noah M. Reid
Maintainer: Dmitry Karabanov <dk@ibiw.ru>
}
\references{
Pons et al. 2006
Reid and Carstens 2012
}
\examples{
# For user data, ape tree format is required. 
# Easiest to input a nexus tree file using ape function "read.nexus".
# In this case, load test data, an ape multiPhylo object.
data(east10)

# Two primary functions for running analysis: 
# bgmyc.singlephy and bgmyc.multiphylo.
# Run bGMYC on a single tree (phylo object).
result.single <- bgmyc.singlephy(phylo = east10[[1]], mcmc = 100, burnin = 1, thinning = 10)

# Visualize MCMC output.
plot(result.single)

# Run bGMYC on multiple trees (multiphylo object).
# Wrapped in \donttest{} to bypass CRAN check parallel worker limits.
\donttest{
result.multi <- bgmyc.multiphylo(multiphylo = east10[1:2], mcmc = 100, burnin = 1, thinning = 10)

# For multiple trees, visualize MCMC output OR pool samples 
# from across trees into a single matrix using plot=FALSE.
plot(result.multi, plot = FALSE) -> parameter.matrix

# Get the posterior probabilities of each possible species cluster.
bgmyc.spec(res = result.multi) -> result.spec

# Produce a sequence by sequence matrix of probabilities of conspecificity.
spec.probmat(result.multi) -> result.probmat

# Plot the matrix of conspecificity probabilities.
# Provide a single tree to order the matrix.
plot(result.probmat, east10[[1]])
}
}
\keyword{package}
