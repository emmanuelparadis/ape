\name{where}
\alias{where}
\title{Find Patterns in DNA Sequences}
\description{
  This function finds patterns in a single or a set of DNA sequences.
}
\usage{
where(x, pattern)
}
\arguments{
  \item{x}{an object of class \code{"DNAbin"}.}
  \item{pattern}{a character string to be searched in \code{x}.}
}
\details{
  If \code{x} is a vector, the function returns a single vector giving
  the position(s) where the pattern was found. If \code{x} is a matrix
  or a list, it returns a list with the positions of the pattern for
  each sequence.

  Patterns may be overlapping. For instance, if \code{pattern = "tata"}
  and the sequence starts with `tatata', then the output will be c(1, 3).
}
\value{
  a vector of integers or a list of such vectors.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{DNAbin}}, \code{\link{image.DNAbin}}
}
\examples{
data(woodmouse)
where(woodmouse, "tata")
}
\keyword{manip}
