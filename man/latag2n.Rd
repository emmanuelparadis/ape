\name{latag2n}
\alias{latag2n}
\title{Leading and Trailing Alignment Gaps to N}
\description{
  Substitutes leading and trailing alignment gaps in aligned sequences
  into \code{N} (i.e., A, C, G, or T). The gaps in the middle of the
  sequences are left unchanged.
}
\usage{
latag2n(x)
}
\arguments{
  \item{x}{an object of class \code{"DNAbin"} with the aligned sequences.}
}
\details{
  This function is called by others in \pkg{ape} and in \pkg{pegas}. It
  is documented here in case it needs to be called by other packages.
}
\value{an object of class \code{"DNAbin"}.}
\author{Emmanuel Paradis}
\seealso{\code{\link{DNAbin}}}
\examples{
x <- as.DNAbin(matrix(c("-", "A", "G", "-", "T", "C"), 2, 3))
y <- latag2n(x)
alview(x)
alview(y)
}
\keyword{manip}
