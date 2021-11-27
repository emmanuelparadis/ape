\name{degree}
\alias{degree}
\alias{degree.phylo}
\alias{degree.evonet}
\title{Vertex Degrees in Trees and Networks}
\description{
  \code{degree} is a generic function to calculate the degree of all
  nodes in a tree or in a network.
}
\usage{
degree(x, ...)
\method{degree}{phylo}(x, details = FALSE, ...)
\method{degree}{evonet}(x, details = FALSE, ...)
}
\arguments{
  \item{x}{an object (tree, network, \dots).}
  \item{details}{whether to return the degree of each node in the tree,
    or a summary table (the default).}
  \item{\dots}{arguments passed to methods.}
}
\details{
  The degree of a node (or vertex) in a network is defined by the number
  of branches (or edges) that connect to this node. In a phylogenetic
  tree, the tips (or terminal nodes) are of degree one, and the
  (internal) nodes are of degree two or more.

  There are currently two methods for the classes \code{"phylo"} and
  \code{"evonet"}. The default of these functions is to return a summary
  table with the degrees observed in the tree or network in the first
  column, and the number of nodes in the second column. If \code{details
  = TRUE}, a vector giving the degree of each node (as numbered in the
  \code{edge} matrix) is returned.

  The validity of the object is not checked, so \code{degree} can be
  used to check problems with badly conformed trees.
}
\value{
  a data frame if \code{details = FALSE}, or a vector of integers
  otherwise.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{checkValidPhylo}}}
\examples{
data(bird.orders)
degree(bird.orders)
degree(bird.orders, details = TRUE)

data(bird.families)
degree(bird.families)

degree(rtree(10)) # 10, 1, 8
degree(rtree(10, rooted = FALSE)) # 10, 0, 8
degree(stree(10)) # 10 + 1 node of degree 10
}
\keyword{manip}
