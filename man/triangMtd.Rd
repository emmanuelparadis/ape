\name{triangMtd}
\alias{triangMtd}
\alias{triangMtds}
\title{Tree Reconstruction Based on the Triangles Method}
\usage{
triangMtd(X)
triangMtds(X)
}
\arguments{
  \item{X}{a distance matrix}.
}
\description{
  Fast distance-based construction method. Should only be used when
  distance measures are fairly reliable.
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  \url{http://archive.numdam.org/ARCHIVE/RO/RO_2001__35_2/RO_2001__35_2_283_0/RO_2001__35_2_283_0.pdf}
}
\author{Andrei Popescu}
\seealso{
  \code{\link{nj}}, \code{\link{bionj}}, \code{\link{fastme}},
  \code{\link{njs}}, \code{\link{mvrs}}, \code{\link{SDM}}
}
\examples{
data(woodmouse)
tr <- triangMtd(dist.dna(woodmouse))
plot(tr)
}
\keyword{models}
