\name{branching.times}
\alias{branching.times}
\title{Branching Times of a Phylogenetic Tree}
\usage{
branching.times(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
}
\description{
  This function computes the branching times of a phylogenetic tree,
  that is the distance from each node to the tips, under the assumption that
  the tree is ultrametric. Note that the function does not check that the
  tree is effectively ultrametric, so if it is not, the returned result
  may not be meaningful.
}
\value{
  a numeric vector with the branching times. If the phylogeny \code{phy}
  has an element \code{node.label}, this is used as names for the
  returned vector; otherwise the numbers (of mode character) of the
  matrix \code{edge} of \code{phy} are used as names.
}
\author{Emmanuel Paradis}
\seealso{
\code{\link{is.ultrametric}}
}
\keyword{manip}
