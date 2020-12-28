\name{subtreeplot}
\alias{subtreeplot}
\title{Zoom on a Portion of a Phylogeny by Successive Clicks}
\description{
  This function plots simultaneously a whole phylogenetic tree
  (supposedly large) and a portion of it determined by clicking on the nodes of the phylogeny. On exit, returns the last subtree visualized.
}
\usage{
subtreeplot(x, wait=FALSE, ...)
}

\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{wait}{a logical indicating whether the node beeing processed should be printed (useful for big phylogenies).}
  \item{\dots}{further arguments passed to \code{plot.phylo}.}
}
\details{
  This function aims at easily exploring very large trees. The main argument is
  a phylogenetic tree, and the second one is a logical indicating whether a waiting message should be printed while the calculation is being processed. 

  The whole tree is plotted on the left-hand side in half of the device. The
  subtree is plotted on the right-hand side in the other half. The user clicks on the nodes in the complete tree and the subtree corresponding to this node is ploted in the right-hand side. There is no limit for the number of clicks that can be done. On exit, the subtree on the right hand side is returned. 

  To use a subtree as the new tree in which to zoom, the user has to use the function many times. This can however be done in a single command line (see example 2).
}
\author{Damien de Vienne \email{damien.de-vienne@u-psud.fr}}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{drop.tip}}, \code{\link{subtrees}}
}
\examples{
\dontrun{
#example 1: simple
tree1 <- rtree(50)
tree2 <- subtreeplot(tree1, wait = TRUE) # on exit, tree2 will be a subtree of tree1

#example 2: more than one zoom
tree1 <- rtree(60)
tree2 <- subtreeplot(subtreeplot(subtreeplot(tree1))) # allow three succssive zooms
}
}
\keyword{hplot}
