\name{pcoa}
\alias{pcoa}
\alias{biplot.pcoa}
\title{ Principal Coordinate Analysis }
\description{
Function \code{\link{pcoa}} computes principal coordinate decomposition (also called classical scaling) of a distance matrix D (Gower 1966). It implements two correction methods for negative eigenvalues.
}
\usage{
pcoa(D, correction="none", rn=NULL)

\method{biplot}{pcoa}(x, Y=NULL, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1, rn=NULL, main=NULL, ...)
}

\arguments{
  \item{D}{ A distance matrix of class \code{dist} or \code{matrix}. }
  \item{correction}{ Correction methods for negative eigenvalues (details below): \code{"lingoes"} and \code{"cailliez"}. Default value: \code{"none"}.  }
  \item{rn}{ An optional vector of row names, of length n, for the n objects. }
  \item{x}{ Output object from \code{\link{pcoa}}. }
  \item{Y}{ Any rectangular data table containing explanatory variables to be projected onto the ordination plot. That table may contain, for example, the community composition data used to compute D, or any transformation of these data; see examples. }
  \item{plot.axes}{ The two PCoA axes to plot. }
  \item{dir.axis1}{ = -1 to revert axis 1 for the projection of points and variables. Default value: +1. }
  \item{dir.axis2}{ = -1 to revert axis 2 for the projection of points
    and variables. Default value: +1. }
  \item{main}{An optional title.}
  \item{...}{ Other graphical arguments passed to function. }
}

\details{
This function implements two methods for correcting for negative values in principal coordinate analysis (PCoA). Negative eigenvalues can be produced in PCoA when decomposing distance matrices produced by coefficients that are not Euclidean (Gower and Legendre 1986, Legendre and Legendre 1998).

In \code{pcoa}, when negative eigenvalues are present in the decomposition results, the distance matrix D can be modified using either the Lingoes or the Cailliez procedure to produce results without negative eigenvalues.

In the Lingoes (1971) procedure, a constant c1, equal to twice absolute value of the largest negative value of the original principal coordinate analysis, is added to each original squared distance in the distance matrix, except the diagonal values. A newe principal coordinate analysis, performed on the modified distances, has at most (n-2) positive eigenvalues, at least 2 null eigenvalues, and no negative eigenvalue.

In the Cailliez (1983) procedure, a constant c2 is added to the original distances in the distance matrix, except the diagonal values. The calculation of c2 is described in Legendre and Legendre (1998). A new principal coordinate analysis, performed on the modified distances, has at most (n-2) positive eigenvalues, at least 2 null eigenvalues, and no negative eigenvalue.

In all cases, only the eigenvectors corresponding to positive eigenvalues are shown in the output list. The eigenvectors are scaled to the square root of the corresponding eigenvalues. Gower (1966) has shown that eigenvectors scaled in that way preserve the original distance (in the D matrix) among the objects. These eigenvectors can be used to plot ordination graphs of the objects.

We recommend not to use PCoA to produce ordinations from the chord, chi-square, abundance profile, or Hellinger distances. It is easier to first transform the community composition data using the following transformations, available in the \code{decostand} function of the \code{vegan} package, and then carry out a principal component analysis (PCA) on the transformed data:

\describe{
\item{ }{Chord transformation: decostand(spiders,"normalize") }
\item{ }{Transformation to relative abundance profiles: decostand(spiders,"total") }
\item{ }{Hellinger transformation: decostand(spiders,"hellinger") }
\item{ }{Chi-square transformation: decostand(spiders,"chi.square") }
}

The ordination results will be identical and the calculations shorter. This two-step ordination method, called transformation-based PCA (tb-PCA), was described by Legendre and Gallagher (2001).

The \code{biplot.pcoa} function produces plots for any pair of principal coordinates. The original variables can be projected onto the ordination plot.
}

\value{

  \item{correction }{The values of parameter \code{correction} and variable 'correct' in the function. }
  \item{note }{A note describing the type of correction done, if any. }
  \item{values }{The eigenvalues and related information: }
  \item{Eigenvalues}{All eigenvalues (positive, null, negative). }
  \item{Relative_eig}{Relative eigenvalues. }
  \item{Corr_eig}{Corrected eigenvalues (Lingoes correction); Legendre and Legendre (1998, p. 438, eq. 9.27). }
  \item{Rel_corr_eig}{Relative eigenvalues after Lingoes or Cailliez correction. }
  \item{Broken_stick}{Expected fractions of variance under the broken stick model. }
  \item{Cumul_eig}{Cumulative relative eigenvalues. }
  \item{Cum_corr_eig}{Cumulative corrected relative eigenvalues. }
  \item{Cumul_br_stick}{Cumulative broken stick fractions. }
  \item{vectors}{The principal coordinates with positive eigenvalues. }
  \item{trace}{The trace of the distance matrix. This is also the sum of all eigenvalues, positive and negative. }
  \item{vectors.cor }{The principal coordinates with positive eigenvalues from the distance matrix corrected using the method specified by parameter \code{correction}. }
  \item{trace.cor }{The trace of the corrected distance matrix. This is also the sum of its eigenvalues. }
}

\references{
Cailliez, F. (1983) The analytical solution of the additive constant problem. \emph{Psychometrika}, \bold{48}, 305--308.

Gower, J. C. (1966) Some distance properties of latent root and vector methods used in multivariate analysis. \emph{Biometrika}, \bold{53}, 325--338.

Gower, J. C. and Legendre, P. (1986) Metric and Euclidean properties of
dissimilarity coefficients. \emph{Journal of Classification}, \bold{3}, 5--48.

Legendre, P. and Gallagher, E. D. (2001) Ecologically meaningful transformations for ordination of species data. \emph{Oecologia}, \bold{129}, 271--280.

Legendre, P. and Legendre, L. (1998) \emph{Numerical Ecology, 2nd English
  edition.} Amsterdam: Elsevier Science BV.

Lingoes, J. C. (1971) Some boundary conditions for a monotone analysis of symmetric matrices. \emph{Psychometrika}, \bold{36}, 195--203.
}

\author{ Pierre Legendre, Universite de Montreal }

\examples{
# Oribatid mite data from Borcard and Legendre (1994)

\dontrun{
if (require(vegan)) {
data(mite)  # Community composition data, 70 peat cores, 35 species

# Select rows 1:30. Species 35 is absent from these rows. Transform to log
mite.log <- log(mite[1:30,-35]+1)  # Equivalent: log1p(mite[1:30,-35])

# Principal coordinate analysis and simple ordination plot
mite.D <- vegdist(mite.log, "bray")
res <- pcoa(mite.D)
res$values
biplot(res)

# Project unstandardized and standardized species on the PCoA ordination plot
mite.log.st = apply(mite.log, 2, scale, center=TRUE, scale=TRUE)

par(mfrow=c(1,2))
biplot(res, mite.log)
biplot(res, mite.log.st)

# Reverse the ordination axes in the  plot
par(mfrow=c(1,2))
biplot(res, mite.log, dir.axis1=-1, dir.axis2=-1)
biplot(res, mite.log.st, dir.axis1=-1, dir.axis2=-1)
}
}
}

\keyword{ multivariate }
