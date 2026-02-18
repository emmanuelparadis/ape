\name{howmanytrees}
\alias{howmanytrees}
\alias{LargeNumber}
\alias{print.LargeNumber}
\title{Calculate Numbers of Phylogenetic Trees}
\description{
  \code{howmanytrees} calculates the number of possible phylogenetic
  trees for a given number of tips.

  \code{LargeNumber} is a utility function to compute (approximately)
  large numbers from the power \eqn{a^b}.
}
\usage{
howmanytrees(n, rooted = TRUE, binary = TRUE,
             labeled = TRUE, detail = FALSE)
LargeNumber(a, b)
\method{print}{LargeNumber}(x, latex = FALSE, digits = 1, ...)
}
\arguments{
  \item{n}{a positive numeric integer giving the number of tips.}
  \item{rooted}{a logical indicating whether the trees are rooted
    (default is \code{TRUE}).}
  \item{binary}{a logical indicating whether the trees are bifurcating
    (default is \code{TRUE}).}
  \item{labeled}{a logical indicating whether the trees have tips
    labeled (default is \code{TRUE}).}
  \item{detail}{a logical indicating whether the eventual intermediate
    calculations should be returned (default is \code{FALSE}). This
    applies only for the multifurcating trees, and the bifurcating,
    rooted, unlabeled trees (aka tree shapes).}
  \item{a, b}{two numbers.}
  \item{x}{an object of class \code{"LargeNumber"}.}
  \item{latex}{a logical value specifying whether to print the number in
    LaTeX code in addition to return it.}
  \item{digits}{the number of digits printed for the real part of the
    large number (unused if \code{latex = FALSE}).}
  \item{\dots}{(unused).}
}
\details{
  In the cases of labeled binary trees, the calculation is done directly
  and a single numeric value is returned (or an object of class
  \code{"LargeNumber"}).

  For multifurcating trees, and bifurcating, rooted, unlabeled trees,
  the calculation is done iteratively for 1 to \code{n} tips. Thus the
  user can print all the intermediate values if \code{detail = TRUE}, or
  only a single value if \code{detail = FALSE} (the default).

  For multifurcating trees, if \code{detail = TRUE}, a matrix is
  returned with the number of tips as rows (named from \code{1} to
  \code{n}), and the number of nodes as columns (named from \code{1} to
  \code{n - 1}). For bifurcating, rooted, unlabeled trees, a vector is
  returned with names equal to the number of tips (from \code{1} to
  \code{n}).

  The number of unlabeled trees (aka tree shapes) can be computed only
  for the rooted binary cases.

  Note that if an infinite value (\code{Inf}) is returned this does not
  mean that there is an infinite number of trees (this cannot be if the
  number of tips is finite), but that the calculation is beyond the
  limits of the computer. Only for the cases of rooted, binary, labeled
  topologies an approximate number is returned in the form a
  \code{"LargeNumber"} object.
}
\value{
  a single numeric value, an object of class \code{"LargeNumber"}, or in
  the case where \code{detail = TRUE} is used, a named vector or
  matrix.
}
\references{
  Felsenstein, J. (2004) \emph{Inferring Phylogenies}. Sunderland:
  Sinauer Associates.
}
\author{Emmanuel Paradis}
\seealso{
  The vignette "RandomTopologies" gives more details and many examples
  of these calculations.
}
\examples{
### Table 3.1 in Felsenstein 2004:
for (i in c(1:20, 30, 40, 50))
  cat(paste(i, howmanytrees(i), sep = "\t"), sep ="\n")
### Table 3.6:
howmanytrees(8, binary = FALSE, detail = TRUE)
}
\keyword{arith}
\keyword{math}
