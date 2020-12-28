\name{rDNAbin}
\alias{rDNAbin}
\title{Random DNA Sequences}
\description{
  This function generates random sets of DNA sequences.
}
\usage{
rDNAbin(n, nrow, ncol, base.freq = rep(0.25, 4), prefix = "Ind_")
}
\arguments{
\item{n}{a vector of integers giving the lengths of the sequences. Can
  be missing in which case \code{nrow} and \code{ncol} must be given.}
\item{nrow, ncol}{two single integer values giving the number of
  sequences and the number of sites, respectively (ignored if \code{n}
  is given).}
\item{base.freq}{the base frequencies.}
\item{prefix}{the prefix used to give labels to the sequences; by
  default these are Ind_1, \dots Ind_n (or Ind_nrow).}
}
\details{
  If \code{n} is used, this function generates a list with sequence lengths given by the values in \code{n}. If \code{n} is missing, a matrix is
  generated.

  The purpose of this function is to generate a set of sequences of a
  specific size. To simulate sequences on a phylogenetic tree, see
  \code{\link[phangorn]{simSeq}} in \pkg{phangorn} (very efficient), and
  the package \pkg{phylosim} (more for pedagogy).
}
\value{
  an object of class \code{"DNAbin"}.
}
\note{
  It is not recommended to use this function to generate objects larger
  than two billion bases (2 Gb).
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{DNAbin}}
}
\examples{
rDNAbin(1:10)
rDNAbin(rep(10, 10))
rDNAbin(nrow = 10, ncol = 10)
}
\keyword{datagen}
