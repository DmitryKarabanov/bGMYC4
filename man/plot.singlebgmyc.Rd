\name{plot.singlebgmyc}
\alias{plot.singlebgmyc}
\title{
Plot MCMC output.
}
\description{
Plots trace plots of MCMC parameters for a single tree analysis.
}
\usage{
\method{plot}{singlebgmyc}(x, burnin = 0, thinning = 1, ...)
}
\arguments{
\item{x}{
object of class "singlebgmyc", output from bgmyc.singlephy.}
\item{burnin}{
the number of samples to discard as burn-in.}
\item{thinning}{
the interval at which samples are retained from the Markov Chain.}
\item{...}{
additional arguments passed to plot.}
}
\value{
No return value, called for side effects (produces a 2x2 panel of trace plots).
}
\author{
Noah M. Reid
Maintainer: Dmitry Karabanov <dk@ibiw.ru>
}
