\name{subtrees}
\alias{subtrees}
\title{All subtrees of a Phylogenetic Tree}
\usage{
subtrees(tree, wait=FALSE)
}
\arguments{
  \item{tree}{an object of class \code{"phylo"}.}
  \item{wait}{a logical indicating whether the node beeing processed should be printed (useful for big phylogenies).}
}
\description{
  This function returns a list of all the subtrees of a phylogenetic tree.
}
\author{Damien de Vienne \email{damien.de-vienne@u-psud.fr}}
\seealso{
  \code{\link{zoom}}, \code{\link{subtreeplot}} for functions extracting particular subtrees.
}
\value{
  \code{subtrees} returns a list of trees of class \code{"phylo"} and
  returns invisibly for each subtree a list with the following
  components:

  \item{tip.label}{}
  \item{node.label}{}
  \item{Ntip}{}
  \item{Nnode}{}
}
\examples{
### Random tree with 12 leaves
phy<-rtree(12)
par(mfrow=c(4,3))
plot(phy, sub="Complete tree")

### Extract the subtrees
l<-subtrees(phy)

### plot all the subtrees
for (i in 1:11) plot(l[[i]], sub=paste("Node", l[[i]]$node.label[1]))
par(mfrow=c(1,1))
}
\keyword{manip}
