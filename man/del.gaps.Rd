\name{del.gaps}
\alias{del.gaps}
\alias{del.colgapsonly}
\alias{del.rowgapsonly}
\title{Delete Alignment Gaps in DNA Sequences}
\description{
  These functions remove gaps (\code{"-"}) in a sample of DNA sequences.
}
\usage{
del.gaps(x)
del.colgapsonly(x, threshold = 1, freq.only = FALSE)
del.rowgapsonly(x, threshold = 1, freq.only = FALSE)
}
\arguments{
  \item{x}{a matrix, a list, or a vector containing the DNA
    sequences; only matrices for \code{del.colgapsonly} and for
    \code{del.rowgapsonly}.}
  \item{threshold}{the largest gap proportion to delete the column or row.}
  \item{freq.only}{if \code{TRUE}, returns only the numbers of gaps for
    each column or row.}
}
\details{
  \code{del.gaps} remove all gaps, so the returned sequences may not
  have all the same lengths and are therefore returned in a list.

  \code{del.colgapsonly} removes the columns with a proportion at least
  \code{threshold} of gaps. Thus by default, only the columns with gaps
  only are removed (useful when a small matrix is extracted from a large
  alignment). \code{del.rowgapsonly} does the same for the rows.

  The sequences can be either in \code{"DNAbin"} or in another format,
  but the returned object is always of class \code{"DNAbin"}.
}
\value{
  \code{del.gaps} returns a vector (if there is only one input sequence)
  or a list of class \code{"DNAbin"}; \code{del.colgapsonly} and
  \code{del.rowgapsonly} return a matrix of class \code{"DNAbin"} or a
  numeric vector (with names for the second function) if \code{freq.only
  = TRUE}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{base.freq}}, \code{\link{seg.sites}},
  \code{\link{image.DNAbin}}, \code{\link{checkAlignment}}
}
\keyword{univar}
