\name{node.depth}
\alias{node.depth}
\alias{node.depth.edgelength}
\alias{node.height}
\title{Depth and Heights of Nodes and Tips}
\description{
  These functions return the depths or heights of nodes and tips.
}
\usage{
node.depth(phy, method = 1)
node.depth.edgelength(phy)
node.height(phy, clado.style = FALSE)
}
\arguments{
  \item{phy}{an object of class "phylo".}
  \item{method}{an integer value (1 or 2); 1: the node depths are
    proportional to the number of tips descending from each node, 2:
    they are evenly spaced.}
  \item{clado.style}{a logical value; if \code{TRUE}, the node heights
    are calculated for a cladogram.}
}
\details{
  \code{node.depth} computes the depth of a node depending on the value
  of \code{method} (see the option \code{node.depth} in
  \code{\link{plot.phylo}}). The value of 1 is given to the tips.

  \code{node.depth.edgelength} does the same but using branch lengths.

  \code{node.height} computes the heights of nodes and tips as plotted
  by a phylogram or a cladogram.
}
\value{
  A numeric vector indexed with the node numbers of the matrix `edge' of
  \code{phy}.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{plot.phylo}}}
\keyword{manip}
