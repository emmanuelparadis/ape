\name{DNAbin2indel}
\alias{DNAbin2indel}
\title{Recode Blocks of Indels}
\description{
  This function scans a set of aligned DNA sequences and returns a
  matrix with information of the localisations and lengths on alignment
  gaps.

  The example below gives an efficient function returning the results in
  a more convenient form.
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
  of the gap, or zero. The rownames are kept in the output.
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
\examples{
## The output of the function below is a matrix with 3 columns where
##   each row  gives the position of each indel in 'x':
## column 1 (or 2): index of the row (or column) of 'x' where the indel
##   starts
## column 3: length of the indel (i.e., number of "-")

indels <- function(x)
{
    tmp <- DNAbin2indel(x)
    stmp <- which(tmp > 0, TRUE)
    cbind(unname(stmp), l = tmp[stmp])
}
}
\keyword{manip}
