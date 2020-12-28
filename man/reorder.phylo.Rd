\name{reorder.phylo}
\alias{reorder.phylo}
\alias{reorder.multiPhylo}
\alias{cladewise}
\alias{postorder}
\title{Internal Reordering of Trees}
\description{
  \code{reorder} changes the internal structure of a phylogeny stored as
  an object of class \code{"phylo"}. The tree returned is the same than
  the one input, but the ordering of the edges could be different.

  \code{cladewise} and \code{postorder} are convenience functions to
  return only the indices of the reordered edge matrices (see examples).
}
\usage{
\method{reorder}{phylo}(x, order = "cladewise", index.only = FALSE, ...)
\method{reorder}{multiPhylo}(x, order = "cladewise", ...)
cladewise(x)
postorder(x)
}
\arguments{
  \item{x}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{order}{a character string: either \code{"cladewise"} (the
    default), \code{"postorder"}, \code{"pruningwise"}, or any
    unambiguous abbreviation of these.}
  \item{index.only}{should the function return only the ordered indices
    of the rows of the edge matrix?}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  Because in a tree coded as an object of class \code{"phylo"} each
  branch is represented by a row in the element `edge', there is an
  arbitrary choice for the ordering of these rows. \code{reorder} allows
  to reorder these rows according to three rules: in the
  \code{"cladewise"} order each clade is formed by a series of
  contiguous rows. In the \code{"postorder"} order, the rows are
  arranged so that computations following pruning-like algorithm the
  tree (or postorder tree traversal) can be done by descending along
  these rows (conversely, a preorder tree traversal can be performed by
  moving from the last to the first row). The \code{"pruningwise"} order
  is an alternative ``pruning'' order which is actually a bottom-up
  traversal order (Valiente 2002). (This third choice might be removed
  in the future as it merely duplicates the second one which is more
  efficient.) The possible multichotomies and branch lengths are preserved.

  Note that for a given order, there are several possible orderings of
  the rows of `edge'.
}
\value{
  an object of class \code{"phylo"} (with the attribute \code{"order"}
  set accordingly), or a numeric vector if \code{index.only = TRUE}; if
  \code{x} is of class \code{"multiPhylo"}, then an object of the same
  class.
}
\references{
  Valiente, G. (2002) \emph{Algorithms on Trees and Graphs.} New York:
  Springer.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.tree}} to read tree files in Newick format,
  \code{\link[stats]{reorder}} for the generic function
}
\examples{
data(bird.families)
tr <- reorder(bird.families, "postorder")
all.equal(bird.families, tr) # uses all.equal.phylo actually
all.equal.list(bird.families, tr) # bypasses the generic

## get the number of descendants for each tip or node:
nr_desc <-  function(x) {
    res <- numeric(max(x$edge))
    res[1:Ntip(x)] <- 1L
    for (i in postorder(x)) {
       tmp <- x$edge[i,1]
       res[tmp] <- res[tmp] + res[x$edge[i, 2]]
   }
   res
}
## apply it to a random tree:
tree <- rtree(10)
plot(tree, show.tip.label = FALSE)
tiplabels()
nodelabels()
nr_desc(tree)
}
\keyword{manip}
