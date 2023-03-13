\name{dnds}
\alias{dnds}
\title{dN/dS Ratio}
\description{
  This function computes the pairwise ratios dN/dS for a set of aligned
  DNA sequences using Li's (1993) method.
}
\usage{
dnds(x, code = 1, codonstart = 1, quiet = FALSE,
     details = FALSE, return.categories = FALSE)
}
\arguments{
\item{x}{an object of class \code{"DNAbin"} (matrix or list) with the
  aligned sequences.}
\item{code}{an integer value giving the genetic code to be
  used. Currently, the codes 1 to 6 are supported.}
\item{codonstart}{an integer giving where to start the translation. This
  should be 1, 2, or 3, but larger values are accepted and have for
  effect to start the translation further within the sequence.}
\item{quiet}{single logical value: whether to indicate progress of
  calculations.}
\item{details}{single logical value (see details).}
\item{return.categories}{a logical value: if \code{TRUE}, a matrix of
  the same size than \code{x} is returned giving the degeneracy category
  of each base in the original alignment.}
}
\details{
  Since \pkg{ape} 5.6, the degeneracy of each codon is calculated
  directly from the genetic code using the function
  \code{\link{trans}}. A consequence is that ambiguous bases are ignored
  (see \code{\link{solveAmbiguousBases}}).

  If \code{details = TRUE}, a table is printed for each pair of
  sequences giving the numbers of transitions and transversions for each
  category of degeneracy (nondegenerate, twofold, and fourfold). This is
  helpful when non-meaningful values are returned (e.g., NaN, Inf,
  negative values).
}
\value{
  an object of class \code{"dist"}, or a numeric matrix if
  \code{return.categories = TRUE}.

}
\references{
  Li, W.-H. (1993) Unbiased estimation of the rates of synonymous and
  nonsynonymous substitution. \emph{Journal of Molecular Evolution},
  \bold{36}, 96--99.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{AAbin}}, \code{\link{trans}}, \code{\link{alview}},
  \code{\link{solveAmbiguousBases}}
}
\examples{
data(woodmouse)
res <- dnds(woodmouse, quiet = TRUE) # NOT correct
res2 <- dnds(woodmouse, code = 2, quiet = TRUE) # using the correct code
identical(res, res2) # FALSE...
cor(res, res2) # ... but very close
## There a few N's in the woodmouse data, but this does not affect
## greatly the results:
res3 <- dnds(solveAmbiguousBases(woodmouse), code = 2, quiet = TRUE)
cor(res, res3)

## a simple example showing the usefulness of 'details = TRUE'
X <- as.DNAbin(matrix(c("C", "A", "G", "G", "T", "T"), 2, 3))
alview(X)
dnds(X, quiet = TRUE) # NaN
dnds(X, details = TRUE) # only a TV at a nondegenerate site
}