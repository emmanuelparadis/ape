\name{solveAmbiguousBases}
\alias{solveAmbiguousBases}
\title{Solve Ambiguous Bases in DNA Sequences}
\description{
  Replaces ambiguous bases in DNA sequences (R, Y, W, \dots) by A, G, C,
  or T.
}
\usage{
solveAmbiguousBases(x, method = "columnwise", random = TRUE)
}
\arguments{
  \item{x}{a matrix of class \code{"DNAbin"}; a list is accepted and is
    converted into a matrix.}
  \item{method}{the method used (no other choice than the default for
    the moment; see details).}
  \item{random}{a logical value (see details).}
}
\details{
  The replacements of ambiguous bases are done columwise. First, the
  base frequencies are counted: if no ambiguous base is found in the
  column, nothing is done. By default (i.e., if \code{random = TRUE}),
  the replacements are done by random sampling using the frequencies of
  the observed compatible, non-ambiguous bases. For instance, if the
  ambiguous base is Y, it is replaced by either C or T using their
  observed frequencies as probabilities. If \code{random = FALSE}, the
  greatest of these frequencies is used. If there are no compatible
  bases in the column, equal probabilities are used. For instance, if
  the ambiguous base is R, and only C and T are observed, then it is
  replaced by either A or G with equal probabilities.

  Alignment gaps are not changed; see the function \code{\link{latag2n}}
  to change the leading and trailing gaps.
}
\value{
  a matrix of class \code{"DNAbin"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{base.freq}}, \code{\link{latag2n}}, \code{\link{dnds}}
}
\examples{
X <- as.DNAbin(matrix(c("A", "G", "G", "R"), ncol = 1))
alview(solveAmbiguousBases(X)) # R replaced by either A or G
alview(solveAmbiguousBases(X, random = FALSE)) # R always replaced by G
}
\keyword{manip}
