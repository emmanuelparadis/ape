\name{bird.families}
\alias{bird.families}
\title{Phylogeny of the Families of Birds From Sibley and Ahlquist}
\description{
  This data set describes the phylogenetic relationships of the families
  of birds as reported by Sibley and Ahlquist (1990). Sibley and
  Ahlquist inferred this phylogeny from an extensive number of DNA/DNA
  hybridization experiments. The ``tapestry'' reported by these two
  authors (more than 1000 species out of the ca. 9000 extant bird
  species) generated a lot of debates.

  The present tree is based on the relationships among families. A few
  families were not included in the figures in Sibley and Ahlquist, and
  thus are not included here as well. The branch lengths were calculated
  from the values of \eqn{\Delta T_{50}H}{Delta T50H} as found in Sibley
  and Ahlquist (1990, figs. 354, 355, 356, and 369).
}
\usage{
data(bird.families)
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
  \code{\link{read.tree}}, \code{\link{bird.orders}}
}
\examples{
data(bird.families)
op <- par()
par(cex = 0.3)
plot(bird.families)
par(op)
}
\keyword{datasets}
