\name{trans}
\alias{trans}
\alias{complement}
\title{Translation from DNA to Amino Acid Sequences}
\description{
  \code{trans} translates DNA sequences into amino acids.
  \code{complement} returns the (reverse) complement sequences.
}
\usage{
trans(x, code = 1, codonstart = 1)
complement(x)
}
\arguments{
\item{x}{an object of class \code{"DNAbin"} (vector, matrix or list).}
\item{code}{an integer value giving the genetic code to be
  used. Currently only the genetic codes  1 to 6 are supported.}
\item{codonstart}{an integer giving where to start the translation. This
  should be 1, 2, or 3, but larger values are accepted and have for
  effect to start the translation further towards the 3'-end of the sequence.}
}
\details{
  With \code{trans}, if the sequence length is not a multiple of three,
  a warning message is printed. Alignment gaps are simply ignored (i.e.,
  \code{AG-} returns \code{X} with no special warning or message). Base
  ambiguities are taken into account where relevant: for instance,
  \code{GGN}, \code{GGA}, \code{GGR}, etc, all return \code{G}.

  See the link given in the References for details about the taxonomic
  coverage and alternative codons of each code.
}
\value{
  an object of class \code{"AAbin"} or \code{"DNAbin"}, respectively.
}
\note{
  These functions are equivalent to \code{translate} and \code{comp} in
  the package \pkg{seqinr} with the difference that there is no need to
  convert the sequences into character strings.
}
\references{
\url{https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes}
}
\author{Emmanuel Paradis}
\seealso{
\code{\link{AAbin}}, \code{\link{checkAlignment}}, \code{\link{alview}}
}
\examples{
data(woodmouse)
X <- trans(woodmouse) # not correct
X2 <- trans(woodmouse, 2) # using the correct code
identical(X, X2)
alview(X[1:2, 1:60]) # some 'Stop' codons (*)
alview(X2[, 1:60])
X2
}