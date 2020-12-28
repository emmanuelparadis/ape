\name{BIONJ}
\alias{bionj}
\title{
  Tree Estimation Based on an Improved Version of the NJ Algorithm
}
\description{
  This function performs the BIONJ algorithm of Gascuel (1997).
}
\usage{
bionj(X)
}
\arguments{
  \item{X}{a distance matrix; may be an object of class \code{"dist"}.}
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  Gascuel, O. (1997) BIONJ: an improved version of the NJ algorithm
  based on a simple model of sequence data.
  \emph{Molecular Biology and Evolution}, \bold{14:}, 685--695.
}
\author{
  original C code by Hoa Sien Cuong and Olivier Gascuel; adapted and
  ported to \R by Vincent Lefort \email{vincent.lefort@lirmm.fr}
}
\seealso{
  \code{\link{nj}}, \code{\link{fastme}}, \code{\link{mvr}},
  \code{\link{bionjs}}, \code{\link{SDM}}, \code{\link{dist.dna}}
}
\examples{
### From Saitou and Nei (1987, Table 1):
x <- c(7, 8, 11, 13, 16, 13, 17, 5, 8, 10, 13,
       10, 14, 5, 7, 10, 7, 11, 8, 11, 8, 12,
       5, 6, 10, 9, 13, 8)
M <- matrix(0, 8, 8)
M[lower.tri(M)] <- x
M <- t(M)
M[lower.tri(M)] <- x
dimnames(M) <- list(1:8, 1:8)
tr <- bionj(M)
plot(tr, "u")
### a less theoretical example
data(woodmouse)
trw <- bionj(dist.dna(woodmouse))
plot(trw)
}
\keyword{models}

