\name{comparePhylo}
\alias{comparePhylo}
\alias{print.comparePhylo}
\title{Compare Two "phylo" Objects}
\description{
  This function compares two phylogenetic trees, rooted or unrooted, and
  returns a detailed report of this comparison.
}
\usage{
comparePhylo(x, y, plot = FALSE, force.rooted = FALSE,
             use.edge.length = FALSE)
\method{print}{comparePhylo}(x, ...)
}
\arguments{
  \item{x, y}{two objects of class \code{"phylo"}.}
  \item{plot}{a logical value. If \code{TRUE}, the two trees are plotted
    on the same device and their similarities are shown.}
  \item{force.rooted}{a logical value. If \code{TRUE}, the trees are
    considered rooted even if \code{is.rooted} returns \code{FALSE}.}
  \item{use.edge.length}{a logical value passed to
    \code{\link{plot.phylo}} (see below).}
  \item{\dots}{unused.}
}
\details{
  In all cases, the numbers of tips and of nodes and the tip labels are
  compared.

  If both trees are rooted, or if \code{force.rooted = TRUE}, the clade
  compositions of each tree are compared. If both trees are also
  ultrametric, their branching times are compared.

  If both trees are unrooted and have the same number of nodes, the
  bipartitions (aka splits) are compared.

  If \code{plot = TRUE}, the edge lengths are not used by default
  because in some situations with unrooted trees, some splits might not
  be visible if the corresponding internal edge length is very short. To
  use edge lengths, set \code{use.edge.length = TRUE}.
}
\value{
  an object of class \code{"comparePhylo"} which is a list with messages
  from the comparison and, optionally, tables comparing branching times.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{all.equal.phylo}}}
\examples{
## two unrooted trees but force comparison as rooted:
a <- read.tree(text = "(a,b,(c,d));")
b <- read.tree(text = "(a,c,(b,d));")
comparePhylo(a, b, plot = TRUE, force.rooted = TRUE)
## two random unrooted trees:
c <- rtree(5, rooted = FALSE)
d <- rtree(5, rooted = FALSE)
comparePhylo(c, d, plot = TRUE)
}
\keyword{manip}
