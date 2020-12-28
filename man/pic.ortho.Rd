\name{pic.ortho}
\alias{pic.ortho}
\title{Phylogenetically Independent Orthonormal Contrasts}
\description{
  This function computes the orthonormal contrasts using the method
  described by Felsenstein (2008). Only a single trait can be analyzed;
  there can be several observations per species.
}
\usage{
pic.ortho(x, phy, var.contrasts = FALSE, intra = FALSE)
}
\arguments{
  \item{x}{a numeric vector or a list of numeric vectors.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{var.contrasts}{logical, indicates whether the expected
    variances of the contrasts should be returned (default to
    \code{FALSE}).}
  \item{intra}{logical, whether to return the intraspecific contrasts.}
}
\details{
  The data \code{x} can be in two forms: a vector if there is a single
  observation for each species, or a list whose elements are vectors
  containing the individual observations for each species. These vectors
  may be of different lengths.

  If \code{x} has names, its values are matched to the tip labels of
  \code{phy}, otherwise its values are taken to be in the same order
  than the tip labels of \code{phy}.
}
\value{
  either a vector of contrasts, or a two-column matrix with the
  contrasts in the first column and their expected variances in the
  second column (if \code{var.contrasts = TRUE}). If the tree has node
  labels, these are used as labels of the returned object.

  If \code{intra = TRUE}, the attribute \code{"intra"}, a list of
  vectors with the intraspecific contrasts or \code{NULL} for the
  species with a one observation, is attached to the returned object.
}
\references{
  Felsenstein, J. (2008) Comparative methods with sampling error and
  within-species variation: Contrasts revisited and revised.
  \emph{American Naturalist}, \bold{171}, 713--725.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{pic}}, \code{\link{varCompPhylip}}
}
\examples{
tr <- rcoal(30)
### a single observation per species:
x <- rTraitCont(tr)
pic.ortho(x, tr)
pic.ortho(x, tr, TRUE)
### different number of observations per species:
x <- lapply(sample(1:5, 30, TRUE), rnorm)
pic.ortho(x, tr, intra = TRUE)
}
\keyword{regression}
