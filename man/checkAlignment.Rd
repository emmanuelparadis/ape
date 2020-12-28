\name{checkAlignment}
\alias{checkAlignment}
\title{Check DNA Alignments}
\description{
  This function performs a series of diagnostics on a DNA alignement.
}
\usage{
checkAlignment(x, check.gaps = TRUE, plot = TRUE, what = 1:4)
}
\arguments{
\item{x}{an object of class \code{"DNAbin"}.}
\item{check.gaps}{a logical value specifying whether to check the
  distribution of alignment gaps.}
\item{plot}{a logical value specifying whether to do the plots.}
\item{what}{an integer value giving the plot to be done. By default,
  four plots are done on the same figure.}
}
\details{
  This function prints on the console a series of diagnostics on the
  set a aligned DNA sequences. If alignment gaps are present, their
  width distribution is analysed, as well as the width of contiguous
  base segments. The pattern of nucleotide diversity on each site is
  also analysed, and a relevant table is printed.

  If \code{plot = TRUE}, four plots are done: an image of the
  alignement, the distribution of gap widths (if present), the Shannon
  index of nucleotide diversity along the sequence, and the number of
  observed bases along the sequence.

  If the sequences contain many gaps, it might be better to set
  \code{check.gaps = FALSE} to skip the analysis of contiguous
  segments.
}
\value{NULL}
\author{Emmanuel Paradis}
\seealso{
\code{\link{alview}}, \code{\link{image.DNAbin}}, \code{\link{all.equal.DNAbin}}
}
\examples{
data(woodmouse)
checkAlignment(woodmouse)
layout(1)
}