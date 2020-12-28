\name{print.phylo}
\alias{print.phylo}
\alias{print.multiPhylo}
\alias{str.multiPhylo}
\title{Compact Display of a Phylogeny}
\usage{
\method{print}{phylo}(x, printlen = 6 ,...)
\method{print}{multiPhylo}(x, details = FALSE ,...)
\method{str}{multiPhylo}(object, ...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{object}{an object of class \code{"multiPhylo"}.}
  \item{printlen}{the number of labels to print (6 by default).}
  \item{details}{a logical indicating whether to print information on
    all trees.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  These functions prints a compact summary of a phylogeny, or a list of
  phylogenies, on the console.
}
\value{
  NULL.
}
\author{Ben Bolker and Emmanuel Paradis}
\seealso{
  \code{\link{read.tree}}, \code{\link{summary.phylo}},
  \code{\link[base]{print}} for the generic \R function
}
\examples{
x <- rtree(10)
print(x)
print(x, printlen = 10)
x <- rmtree(2, 10)
print(x)
print(x, TRUE)
str(x)
}
\keyword{manip}
