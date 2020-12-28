\name{is.binary}
\alias{is.binary}
\alias{is.binary.phylo}
\alias{is.binary.multiPhylo}
\alias{is.binary.tree}
\title{Test for Binary Tree}
\description{
  This function tests whether a phylogenetic tree is binary.
}
\usage{
is.binary(phy)
\method{is.binary}{phylo}(phy)
\method{is.binary}{multiPhylo}(phy)
\method{is.binary}{tree}(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
}
\details{
  The test differs whether the tree is rooted or not. An urooted tree is
  considered binary if all its nodes are of degree three (i.e., three
  edges connect to each node). A rooted tree is considered binary if all
  nodes (including the root node) have exactly two descendant nodes, so
  that they are of degree three expect the root which is of degree 2.

  \code{is.binary.tree} is deprecated and will be removed soon:
  currently it calls \code{is.binary}.
}
\value{
  a logical vector.
}
\seealso{
\code{\link{is.rooted}}, \code{\link{is.ultrametric}}, \code{\link{multi2di}}
}
\author{Emmanuel Paradis}
\examples{
is.binary(rtree(10))
is.binary(rtree(10, rooted = FALSE))
is.binary(stree(10))
x <- setNames(rmtree(10, 10), LETTERS[1:10])
is.binary(x)
}
\keyword{logic}
