\name{bird.orders}
\alias{bird.orders}
\title{Phylogeny of the Orders of Birds From Sibley and Ahlquist}
\description{
  This data set describes the phylogenetic relationships of the orders
  of birds as reported by Sibley and Ahlquist (1990). The branch lengths
  were calculated from the values of \eqn{\Delta T_{50}H}{Delta T50H} as
  found in Sibley and Ahlquist (1990, Fig. 353).
}
\usage{
data(bird.orders)
}
\format{
  The data are stored as an object of class \code{"phylo"} which
  structure is described in the help page of the function
  \code{\link{read.tree}}.
}
\source{
  Sibley, C. G. and Ahlquist, J. E. (1990) Phylogeny and classification
  of birds: a study in molecular evolution. New Haven: Yale University Press.
}
\seealso{
  \code{\link{read.tree}}, \code{\link{bird.families}}
}
\examples{
data(bird.orders)
plot(bird.orders)
}
\keyword{datasets}
