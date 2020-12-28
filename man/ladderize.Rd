\name{ladderize}
\alias{ladderize}
\title{Ladderize a Tree}
\usage{
ladderize(phy, right = TRUE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{right}{a logical specifying whether the smallest clade is on the
    right-hand side (when the tree is plotted upwards), or the opposite
    (if \code{FALSE}).}
}
\description{
  This function reorganizes the internal structure of the tree to get
  the ladderized effect when plotted.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{reorder.phylo}}
}
\examples{
tr <- rcoal(50)
layout(matrix(1:4, 2, 2))
plot(tr, main = "normal")
plot(ladderize(tr), main = "right-ladderized")
plot(ladderize(tr, FALSE), main = "left-ladderized")
layout(matrix(1, 1))
}
\keyword{manip}
