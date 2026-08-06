\name{bgmyc.spec}
\alias{bgmyc.spec}
\title{
bGMYC output processing.
}
\description{
Processes output from bgmyc.singlephy or bgmyc.multiphylo.
}
\usage{
bgmyc.spec(res, filename = NULL, cmatrix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
output from bgmyc.singlephy or bgmyc.multiphylo}
  \item{filename}{
if given, a file that contains all species sampled in the MCMC and their posterior probabilities}
  \item{cmatrix}{
this option may not work at the moment. a matrix of tree leaves by sampling localities. if given, returns a matrix of species by samping localities suitable for community ecological analysis}
}
\details{
this function outputs a list with one or two components. the first is a data frame containing all species sampled in the MCMC and their posterior probabilities (also optionally output to a text file if "filename" is given) and the second is an optional collection of matrices if "cmatrix" is specified. cmatrix should be a tree leaves by sampling locality matrix and if given, a collection of species by sampling locality matrices are returned that represent the posterior distribution of species occurrences.
}
\value{
A list with one or two components:
\item{specprobs}{a data frame containing all species sampled in the MCMC and their posterior probabilities}
\item{cmatrix}{(optional) a collection of species by sampling locality matrices if cmatrix argument is specified}
}
\author{
Noah M. Reid
Maintainer: Dmitry Karabanov <dk@ibiw.ru>
}
