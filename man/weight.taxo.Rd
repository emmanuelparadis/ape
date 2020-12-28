\name{weight.taxo}
\alias{weight.taxo}
\alias{weight.taxo2}
\title{Define Similarity Matrix}
\usage{
  weight.taxo(x)
  weight.taxo2(x, y)
}
\arguments{
  \item{x, y}{a vector or a factor.}
}
\description{
  \code{weight.taxo} computes a matrix whose entries [i, j] are set to 1
  if x[i] == x[j], 0 otherwise.

  \code{weight.taxo2} computes a matrix whose entries [i, j] are set to 1
  if x[i] == x[j] AND y[i] != y[j], 0 otherwise.

  The diagonal [i, i] is always set to 0.

  The returned matrix can be used as a weight matrix in
  \code{\link{Moran.I}}. \code{x} and \code{y} may be vectors of
  factors.

  See further details in \code{vignette("MoranI")}.
}
\value{
  a square numeric matrix.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{Moran.I}}, \code{\link{correlogram.formula}}
}
\keyword{manip}
