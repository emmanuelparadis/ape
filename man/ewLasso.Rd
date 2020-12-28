\name{ewLasso}
\alias{ewLasso}
\title{
  Incomplete distances and edge weights of unrooted topology
}
\description{
  This function implements a method for checking whether an incomplete
  set of distances satisfy certain conditions that might make it
  uniquely determine the edge weights of a given topology, T. It prints
  information about whether the graph with vertex set the set of leaves,
  denoted by X, and edge set the set of non-missing distance pairs,
  denoted by L, is connected or strongly non-bipartite. It then also
  checks whether L is a triplet cover for T.
}
\usage{
ewLasso(X, phy)
}
\arguments{
  \item{X}{a distance matrix.}
  \item{phy}{an unrooted tree of class \code{"phylo"}.}
}
\details{
  Missing values must be represented by either \code{NA} or a negative value.

  This implements a method for checking whether an incomplete set of
  distances satisfies certain conditions that might make it uniquely
  determine the edge weights of a given topology, T. It prints
  information about whether the graph, G, with vertex set the set of
  leaves, denoted by X, and edge set the set of non-missing distance
  pairs, denoted by L, is connected or strongly non-bipartite. It also
  checks whether L is a triplet cover for T. If G is not connected, then
  T does not need to be the only topology satisfying the input
  incomplete distances. If G is not strongly non-bipartite then the
  edge-weights of the edges of T are not the unique ones for which the
  input distance is satisfied. If L is a triplet cover, then the input
  distance matrix uniquely determines the edge weights of T. See Dress
  et al. (2012) for details.
}
\value{
  NULL, the results are printed in the console.
}
\references{
  Dress, A. W. M., Huber, K. T., and Steel, M. (2012) `Lassoing' a
  phylogentic tree I: basic properties, shellings and covers.
  \emph{Journal of Mathematical Biology}, \bold{65(1)}, 77--105.
}
\author{Andrei Popescu}
\keyword{multivariate}
