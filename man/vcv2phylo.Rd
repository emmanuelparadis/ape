\name{vcv2phylo}
\alias{vcv2phylo}
\title{Variance-Covariance Matrix to Tree}
\description{
  This function transforms a variance-covariance matrix into a
  phylogenetic tree.
}
\usage{
vcv2phylo(mat, tolerance = 1e-7)
}
\arguments{
  \item{mat}{a square symmetric (positive-definite) matrix.}
  \item{tolerance}{the numeric tolerance used to compare the branch
    lengths.}
}
\details{
  The function tests if the matrix is symmetric and positive-definite
  (i.e., all its eigenvalues positive within the specified tolerance).
}
\value{
  an object of class \code{"phylo"}.
}
\author{Simon Blomberg}
\seealso{
  \code{\link{vcv}}, \code{\link{corPhyl}}
}
\examples{
tr <- rtree(10)
V <- vcv(tr) # VCV matrix assuming Brownian motion
z <- vcv2phylo(V)
identical(tr, z) # FALSE
all.equal(tr, z) # TRUE
}
\keyword{manip}
\keyword{multivariate}
