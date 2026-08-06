\name{checkrates}
\alias{checkrates}
\title{
bGMYC output processing.
}
\description{
Processes output from bgmyc.singlephy or bgmyc.multiphylo.
}
\usage{
checkrates(result)
}
\arguments{
  \item{result}{
output from bgmyc.singlephy or bgmyc.multiphylo}
}
\details{
This function takes the results of a bgmyc.singlephy or bgmyc.multiphylo run and outputs a matrix of all parameters, 
including the branching rates for the Yule and coalescent processes (this is the only place to get those last two). 
The branching rates are may be useful in evaluating the reliability of the model output (although I haven't tested this extensively). 
The GMYC is based on the idea that coalescent and speciation events happen at drastically different (perhaps on around an order of magnitude different) rates 
and thus can be differentiated in a phylogenetic tree. If the rates are similar, then the threshold point in the model might be relatively indistinct 
and the results unreliable. Comparison of the rates might help identify this.
}
\value{
An object of class \code{"bgmycrates"}: a matrix with columns:
\item{p.div}{Yule rate change parameter}
\item{p.coal}{coalescent rate change parameter}
\item{threshold}{threshold parameter (number of species)}
\item{lambda.div}{Yule branching rate}
\item{lambda.coal}{coalescent branching rate}
\item{lambda.div.mod}{modified Yule branching rate}
}
\author{
Noah M. Reid
Maintainer: Dmitry Karabanov <dk@ibiw.ru>
}
