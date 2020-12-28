\name{compute.brtime}
\alias{compute.brtime}
\title{Compute and Set Branching Times}
\description{
  This function computes the branch lengths of a tree giving its
  branching times (aka node ages or heights).
}
\usage{
compute.brtime(phy, method = "coalescent", force.positive = NULL)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{method}{either \code{"coalescent"} (the default), or a numeric
    vector giving the branching times.}
  \item{force.positive}{a logical value (see details).}
}
\details{
  By default, a set of random branching times is generated from a simple
  coalescent, and the option \code{force.positive} is set to \code{TRUE}
  so that no branch length is negative.

  If a numeric vector is passed to \code{method}, it is taken as the
  branching times of the nodes with respect to their numbers (i.e., the
  first element of \code{method} is the branching time of the node
  numbered \eqn{n + 1} [= the root], the second element of the node
  numbered \eqn{n + 2}, and so on), so \code{force.positive} is set to
  \code{FALSE}. This may result in negative branch lengths. To avoid
  this, one should use \code{force.positive = TRUE} in which case the
  branching times are eventually reordered.
}
\value{
  An object of class \code{"phylo"} with branch lengths and ultrametric.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{compute.brlen}}, \code{\link{branching.times}}
}
\examples{
tr <- rtree(10)
layout(matrix(1:4, 2))
plot(compute.brtime(tr)); axisPhylo()
plot(compute.brtime(tr, force.positive = FALSE)); axisPhylo()
plot(compute.brtime(tr, 1:9)); axisPhylo() # a bit nonsense
plot(compute.brtime(tr, 1:9, TRUE)); axisPhylo()
layout(1)
}
\keyword{manip}
