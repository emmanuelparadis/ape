\name{dnds}
\alias{dnds}
\title{dN/dS Ratio}
\description{
  This function computes the pairwise ratios dN/dS for a set of aligned
  DNA sequences using Li's (1993) method.
}
\usage{
dnds(x, code = 1, codonstart = 1, quiet = FALSE)
}
\arguments{
\item{x}{an object of class \code{"DNAbin"} (matrix or list) with the
  aligned sequences.}
\item{code}{an integer value giving the genetic code to be
  used. Currently only 1 (standard code) and 2 (vertebrate mitochondrial
  code) are supported.}
\item{codonstart}{an integer giving where to start the translation. This
  should be 1, 2, or 3, but larger values are accepted and have for
  effect to start the translation further within the sequence.}
\item{quiet}{single logical value: whether to indicate progress of
  calculations.}
}
\value{
  an object of class \code{"dist"}.
}
\references{
  Li, W.-H. (1993) Unbiased estimation of the rates of synonymous and
  nonsynonymous substitution. \emph{Journal of Molecular Evolution},
  \bold{36}, 96--99.
}
\author{Emmanuel Paradis}
\seealso{
\code{\link{AAbin}}, \code{\link{trans}}, \code{\link{alview}}
}
\examples{
data(woodmouse)
res <- dnds(woodmouse, quiet = TRUE) # NOT correct
res2 <- dnds(woodmouse, code = 2, quiet = TRUE) # using the correct code
identical(res, res2)
}