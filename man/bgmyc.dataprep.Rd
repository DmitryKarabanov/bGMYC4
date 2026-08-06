\name{bgmyc.dataprep}
\alias{bgmyc.dataprep}
\title{Prepares a tree for analysis.}
\description{
Takes an ape "phylo" object and turns it into input for downstream functions.
}
\usage{
bgmyc.dataprep(tr)
}

\arguments{
  \item{tr}{
an ape "phylo" object. 
}
}
\value{
A list with the following components:
\item{mrca.nodes}{list of MRCA node indices for each threshold}
\item{nod.types}{list of node type vectors for each threshold}
\item{n}{list of the number of coalescent entities for each threshold}
\item{list.s.nod}{list of matrices indicating node states}
\item{list.i.mat}{list of interval matrices used in likelihood calculation}
\item{internod}{vector of internode distances}
\item{tree}{the input phylogenetic tree}
}
\author{
Noah M. Reid
Maintainer: Dmitry Karabanov <dk@ibiw.ru>
}
