\name{seg.sites}
\alias{seg.sites}
\title{
  Find Segregating Sites in DNA Sequences
}
\usage{
seg.sites(x, strict = FALSE, trailingGapsAsN = TRUE)
}
\arguments{
  \item{x}{a matrix or a list which contains the DNA sequences.}
  \item{strict}{a logical value; if \code{TRUE}, ambiguities and gaps in
    the sequences are not interpreted in the usual way.}
  \item{trailingGapsAsN}{a logical value; if \code{TRUE} (the default),
    the leading and trailing alignment gaps are considered as unknown
    bases (i.e., N).}
}
\description{
  This function gives the indices of segregating (polymorphic) sites in
  a sample of DNA sequences.
}
\details{
  If the sequences are in a list, they must all be of the same length.

  If \code{strict = FALSE} (the default), the following rule is used to
  determine if a site is polymorphic or not in the presence of ambiguous
  bases: `A' and `R' are not interpreted as different, `A' and `Y' are
  interpreted as different, and `N' and any other base (ambiguous or
  not) are interpreted as not different. If \code{strict = TRUE}, all
  letters are considered different.

  Alignment gaps are considered different from all letters except for
  the leading and trailing gaps if \code{trailingGapsAsN = TRUE} (which
  is the default).
}
\value{
  A numeric (integer) vector giving the indices of the segregating
  sites.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{base.freq}}, \code{\link[pegas]{theta.s}},
  \code{\link[pegas]{nuc.div}}
}
\examples{
data(woodmouse)
y <- seg.sites(woodmouse)
y
length(y)
}
\keyword{univar}
