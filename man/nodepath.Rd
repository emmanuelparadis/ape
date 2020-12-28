\name{nodepath}
\alias{nodepath}
\title{Find Paths of Nodes}
\description{
  This function finds paths of nodes in a tree. The nodes can be
  internal and/or terminal (i.e., tips).
}
\usage{
nodepath(phy, from = NULL, to = NULL)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{from, to}{integers giving node or tip numbers.}
}
\details{
  By default, this function returns all the paths from the root to each
  tip of the tree. If both arguments \code{from} and \code{to} are
  specified, the shortest path of nodes linking them is returned.
}
\value{
  a list of vectors of integers (by default), or a single vector of
  integers.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{getMRCA}}}
\examples{
tr <- rtree(2)
nodepath(tr)
nodepath(tr, 1, 2)
}
\keyword{manip}
