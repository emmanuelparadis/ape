\name{base.freq}
\alias{base.freq}
\alias{GC.content}
\alias{Ftab}
\title{Base frequencies from DNA Sequences}
\description{
  \code{base.freq} computes the frequencies (absolute or relative) of
  the four DNA bases (adenine, cytosine, guanine, and thymidine) from a
  sample of sequences.

  \code{GC.content} computes the proportion of G+C (using the previous
  function). All missing or unknown sites are ignored.

  \code{Ftab} computes the contingency table with the absolute
  frequencies of the DNA bases from a pair of sequences.
}
\usage{
base.freq(x, freq = FALSE, all = FALSE)
GC.content(x)
Ftab(x, y = NULL)
}
\arguments{
  \item{x}{a vector, a matrix, or a list which contains the DNA
    sequences.}
  \item{y}{a vector with a single DNA sequence.}
  \item{freq}{a logical specifying whether to return the proportions
    (the default) or the absolute frequencies (counts).}
  \item{all}{a logical; by default only the counts of A, C, G, and T are
    returned. If \code{all = TRUE}, all counts of bases, ambiguous codes,
    missing data, and alignment gaps are returned.}
}
\details{
  The base frequencies are computed over all sequences in the
  sample.

  For \code{Ftab}, if the argument \code{y} is given then both \code{x}
  and \code{y} are coerced as vectors and must be of equal length. If
  \code{y} is not given, \code{x} must be a matrix or a list and only
  the two first sequences are used.
}
\value{
  A numeric vector with names \code{c("a", "c", "g", "t")} (and possibly
  \code{"r", "m", ...}, a single numeric value, or a four by four matrix
  with similar dimnames.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{seg.sites}}, \code{\link[pegas]{nuc.div}},
  \code{\link{DNAbin}}
}
\examples{
data(woodmouse)
base.freq(woodmouse)
base.freq(woodmouse, TRUE)
base.freq(woodmouse, TRUE, TRUE)
GC.content(woodmouse)
Ftab(woodmouse)
Ftab(woodmouse[1, ], woodmouse[2, ]) # same than above
Ftab(woodmouse[14:15, ]) # between the last two
}
\keyword{univar}
\keyword{manip}
