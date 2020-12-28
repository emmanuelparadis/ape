\name{varCompPhylip}
\alias{varCompPhylip}
\title{Variance Components with Orthonormal Contrasts}
\description{
  This function calls Phylip's contrast program and returns the
  phylogenetic and phenotypic variance-covariance components for one or
  several traits. There can be several observations per species.
}
\usage{
varCompPhylip(x, phy, exec = NULL)
}
\arguments{
  \item{x}{a numeric vector, a matrix (or data frame), or a list.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{exec}{a character string giving the name of the executable
    contrast program (see details).}
}
\details{
  The data \code{x} can be in several forms: (i) a numeric vector if
  there is single trait and one observation per species; (ii) a
  matrix or data frame if there are several traits (as columns) and a
  single observation of each trait for each species; (iii) a list of
  vectors if there is a single trait and several observations per
  species; (iv) a list of matrices or data frames: same than (ii) but
  with several traits and the rows are individuals.

  If \code{x} has names, its values are matched to the tip labels of
  \code{phy}, otherwise its values are taken to be in the same order
  than the tip labels of \code{phy}.

  Phylip (version 3.68 or higher) must be accessible on your computer. If
  you have a Unix-like operating system, the executable name is assumed
  to be \code{"phylip contrast"} (as in Debian); otherwise it is set
  to \code{"contrast"}. If this doesn't suit your system, use the
  option \code{exec} accordingly. If the executable is not in the path, you
  may need to specify it, e.g., \code{exec = "C:/Program Files/Phylip/contrast"}.
}
\value{
  a list with elements \code{varA} and \code{varE} with the phylogenetic
  (additive) and phenotypic (environmental) variance-covariance
  matrices. If a single trait is analyzed, these contains its variances.
}
\references{
  Felsenstein, J. (2004) Phylip (Phylogeny Inference Package) version
  3.68. Department of Genetics, University of Washington, Seattle, USA.
  \url{http://evolution.genetics.washington.edu/phylip/phylip.html}.

  Felsenstein, J. (2008) Comparative methods with sampling error and
  within-species variation: Contrasts revisited and revised.
  \emph{American Naturalist}, \bold{171}, 713--725.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{pic}}, \code{\link{pic.ortho}}, \code{\link{compar.lynch}}
}
\examples{
\dontrun{
tr <- rcoal(30)
### Five traits, one observation per species:
x <- replicate(5, rTraitCont(tr, sigma = 1))
varCompPhylip(x, tr) # varE is small
x <- replicate(5, rnorm(30))
varCompPhylip(x, tr) # varE is large
### Five traits, ten observations per species:
x <- replicate(30, replicate(5, rnorm(10)), simplify = FALSE)
varCompPhylip(x, tr)
}}
\keyword{regression}
