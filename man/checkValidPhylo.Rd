\name{checkValidPhylo}
\alias{checkValidPhylo}
\title{Check the Structure of a "phylo" Object}
\description{
  This function takes as single argument an object (phy), checks its
  elements, and prints a diagnostic. All problems are printed with a
  label: FATAL (will likely cause an error or a crash) or MODERATE (may
  cause some problems).

  This function is mainly intended for developers creating
  \code{"phylo"} objects from scratch.
}
\usage{
checkValidPhylo(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
}
\value{
  NULL.
}
\author{Emmanuel Paradis}
\examples{
tr <- rtree(3)
checkValidPhylo(tr)
tr$edge[1] <- 0
checkValidPhylo(tr)
}
\keyword{manip}
