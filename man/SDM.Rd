\name{SDM}
\alias{SDM}
\title{Construction of Consensus Distance Matrix With SDM}
\description{
  This function implements the SDM method of Criscuolo et al. (2006) for
  a set of n distance matrices.
}
\usage{
SDM(...)
}
\arguments{
  \item{\dots}{2n elements (with n > 1), the first n elements are the
    distance matrices: these can be (symmetric) matrices, objects of
    class \code{"dist"}, or a mix of both. The next n elements are the
    sequence length from which the matrices have been estimated (can be
    seen as a degree of confidence in matrices).}
}
\details{
  Reconstructs a consensus distance matrix from a set of input distance
  matrices on overlapping sets of taxa. Potentially missing values in
  the supermatrix are represented by \code{NA}. An error is returned if
  the input distance matrices can not resolve to a consensus matrix.
}
\value{
  a 2-element list containing a distance matrix labelled by the union of
  the set of taxa of the input distance matrices, and a variance matrix
  associated to the returned distance matrix.
}
\references{
  Criscuolo, A., Berry, V., Douzery, E. J. P. , and Gascuel, O. (2006)
  SDM: A fast distance-based approach for (super)tree building in
  phylogenomics. \emph{Systematic Biology}, \bold{55}, 740--755.
}
\author{Andrei Popescu}
\seealso{
  \code{\link{bionj}}, \code{\link{fastme}}, \code{\link{njs}},
  \code{\link{mvrs}}, \code{\link{triangMtd}}
}
\keyword{models}
