\name{is.ultrametric}
\alias{is.ultrametric}
\alias{is.ultrametric.phylo}
\alias{is.ultrametric.multiPhylo}
\title{Test if a Tree is Ultrametric}
\description{
  This function tests whether a tree is ultrametric using the distances
  from each tip to the root.
}
\usage{
is.ultrametric(phy, ...)
\method{is.ultrametric}{phylo}(phy, tol = .Machine$double.eps^0.5, option = 1, ...)
\method{is.ultrametric}{multiPhylo}(phy, tol = .Machine$double.eps^0.5, option = 1, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{tol}{a numeric >= 0, variation below this value are considered
    non-significant.}
  \item{option}{an integer (1 or 2; see details).}
  \item{\dots}{arguments passed among methods.}
}
\details{
  The test is based on the distances from each tip to the root and a
  criterion: if \code{option = 1}, the criterion is the scaled range
  ((max - min/max)), if \code{option = 2}, the variance is used (this
  was the method used until ape 3.5). The default criterion is invariant
  to linear changes of the branch lengths.
}
\value{
  a logical vector.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{is.binary}}, \code{\link[base]{.Machine}}
}
\examples{
is.ultrametric(rtree(10))
is.ultrametric(rcoal(10))
}
\keyword{utilities}
