\name{c.phylo}
\alias{c.phylo}
\alias{c.multiPhylo}
\alias{.compressTipLabel}
\alias{.uncompressTipLabel}
\title{Building Lists of Trees}
\description{
  These functions help to build lists of trees of class \code{"multiPhylo"}.
}
\usage{
\method{c}{phylo}(..., recursive = TRUE)
\method{c}{multiPhylo}(..., recursive = TRUE)
.compressTipLabel(x, ref = NULL)
.uncompressTipLabel(x)
}
\arguments{
  \item{\dots}{one or several objects of class \code{"phylo"} and/or
    \code{"multiPhylo"}.}
  \item{recursive}{see details.}
  \item{x}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{ref}{an optional vector of mode character to constrain the order
    of the tips. By default, the order from the first tree is used.}
}
\details{
  These \code{c} methods check all the arguments, and return by default
  a list of single trees unless some objects are not trees or lists of
  trees, in which case \code{recursive} is switched to FALSE and a
  warning message is given. If \code{recursive = FALSE}, the objects are
  simply concatenated into a list. Before \pkg{ape} 4.0, \code{recursive}
  was always set to FALSE.

  \code{.compressTipLabel} transforms an object of class
  \code{"multiPhylo"} by checking that all trees have the same tip
  labels and renumbering the tips in the \code{edge} matrix so that the
  tip numbers are also the same taking the first tree as the reference
  (duplicated labels are not allowed). The returned object has a unique
  vector of tip labels (\code{attr(x, "TipLabel")}).

  \code{.uncompressTipLabel} does the reverse operation.
}
\value{
  An object of class \code{"multiPhylo"}.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{summary.phylo}}, \code{\link{multiphylo}}}
\examples{
x <- c(rtree(4), rtree(2))
x
y <- c(rtree(4), rtree(4))
z <- c(x, y)
z
print(z, TRUE)
try(.compressTipLabel(x)) # error
a <- .compressTipLabel(y)
.uncompressTipLabel(a) # back to y
## eventually compare str(a) and str(y)
}
\keyword{manip}
