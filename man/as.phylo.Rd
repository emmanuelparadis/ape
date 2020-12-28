\name{as.phylo}
\alias{as.phylo}
\alias{as.phylo.hclust}
\alias{as.phylo.phylog}
\alias{as.hclust.phylo}
\alias{old2new.phylo}
\alias{new2old.phylo}
\alias{as.network.phylo}
\alias{as.igraph}
\alias{as.igraph.phylo}
\title{Conversion Among Tree and Network Objects}
\description{
  \code{as.phylo} is a generic function which converts an object into a
  tree of class \code{"phylo"}. There are currently two methods for
  objects of class \code{"hclust"} and of class \code{"phylog"}
  (implemented in the package ade4).

  \code{as.hclust.phylo} is a method of the generic
  \code{\link[stats]{as.hclust}} which converts an object of class
  \code{"phylo"} into one of class \code{"hclust"}. This can used to
  convert an object of class \code{"phylo"} into one of class
  \code{"dendrogram"} (see examples).

  \code{as.network} and \code{as.igraph} convert trees of class
  \code{"phylo"} into these respective classes defined in the packages
  of the same names (where the generics are defined).

  \code{old2new.phylo} and \code{new2old.phylo} are utility functions
  for converting between the old and new coding of the class
  \code{"phylo"}.
}
\usage{
as.phylo(x, ...)
\method{as.phylo}{hclust}(x, ...)
\method{as.phylo}{phylog}(x, ...)
\method{as.hclust}{phylo}(x, ...)
old2new.phylo(phy)
new2old.phylo(phy)
\method{as.network}{phylo}(x, directed = is.rooted(x), ...)
\method{as.igraph}{phylo}(x, directed = is.rooted(x), use.labels = TRUE, ...)
}
\arguments{
  \item{x}{an object to be converted into another class.}
  \item{directed}{a logical value: should the network be directed? By
    default, this depends on whether the tree is rooted or not.}
  \item{use.labels}{a logical specifying whether to use labels to build
    the network of class \code{"igraph"}. If \code{TRUE} and the tree
    has no node labels, then some default labels are created first. If
    \code{FALSE}, the network is built with integers.}
  \item{\dots}{further arguments to be passed to or from other methods.}
  \item{phy}{an object of class \code{"phylo"}.}
}
\value{
  An object of class \code{"hclust"}, \code{"phylo"}, \code{"network"},
  or \code{"igraph"}.
}
\note{
  In an object of class \code{"hclust"}, the \code{height} gives the
  distance between the two sets that are being agglomerated. So these
  distances are divided by two when setting the branch lengths of a
  phylogenetic tree.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link[stats]{hclust}}, \code{\link[stats]{as.hclust}},
  \code{\link[stats]{dendrogram}}, \code{\link[ade4]{phylog}},
  \code{\link{as.phylo.formula}}
}
\examples{
data(bird.orders)
hc <- as.hclust(bird.orders)
tr <- as.phylo(hc)
all.equal(bird.orders, tr) # TRUE

### shows the three plots for tree objects:
dend <- as.dendrogram(hc)
layout(matrix(c(1:3, 3), 2, 2))
plot(bird.orders, font = 1)
plot(hc)
par(mar = c(8, 0, 0, 0)) # leave space for the labels
plot(dend)

### how to get identical plots with
### plot.phylo and plot.dendrogram:
layout(matrix(1:2, 2, 1))
plot(bird.orders, font = 1, no.margin = TRUE, label.offset = 0.4)
par(mar = c(0, 0, 0, 8))
plot(dend, horiz = TRUE)
layout(1)

\dontrun{
### convert into networks:
if (require(network)) {
    x <- as.network(rtree(10))
    print(x)
    plot(x, vertex.cex = 1:4)
    plot(x, displaylabels = TRUE)
}
tr <- rtree(5)
if (require(igraph)) {
    print((x <- as.igraph(tr)))
    plot(x)
    print(as.igraph(tr, TRUE, FALSE))
    print(as.igraph(tr, FALSE, FALSE))
}
}
}
\keyword{manip}
