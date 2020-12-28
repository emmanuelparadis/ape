\name{DNAbin2indel}
\alias{DNAbin2indel}
\title{Recode Blocks of Indels}
\description{
  This function scans a set of aligned DNA sequences and returns a
  matrix with information of the localisations and lengths on alignment
  gaps.
}
\usage{
DNAbin2indel(x)
}
\arguments{
  \item{x}{an object of class \code{"DNAbin"}.}
}
\details{
  The output matrix has the same dimensions than the input one with,
  either a numeric value where an alignment gap starts giving the length
  of the gap, or zero. The rownames are kept.
}
\value{
  a numeric matrix.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{DNAbin}}, \code{\link{as.DNAbin}}, \code{\link{del.gaps}},
  \code{\link{seg.sites}}, \code{\link{image.DNAbin}},
  \code{\link{checkAlignment}}
}
\keyword{manip}
