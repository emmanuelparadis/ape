\name{plot.phylo.extra}
\alias{plot.phylo.extra}
\alias{plotBreakLongEdges}
\alias{drawSupportOnEdges}
\title{Extra Fuctions to Plot and Annotate Phylogenies}
\description{
  These are extra functions to plot and annotate phylogenies, mostly
  calling basic graphical functions in \pkg{ape}.
}
\usage{
plotBreakLongEdges(phy, n = 1, ...)
drawSupportOnEdges(value, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{n}{the numner of long branches to be broken.}
  \item{value}{the values to be printed on the internal branches of the tree.}
  \item{\dots}{further arguments to be passed to \code{plot.phylo} or to
    \code{edgelabels}.}
}
\details{
  \code{drawSupportOnEdges} assumes the tree is unrooted, so the vector
  \code{value} should have as many values than the number of internal
  branches (= number of nodes - 1). If there is one additional value, it
  is assumed that it relates to the root node and is dropped (see examples).
}
\value{NULL}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{edgelabels}},
  \code{\link{boot.phylo}}, \code{\link{plotTreeTime}}
}
\examples{
tr <- rtree(10)
tr$edge.length[c(1, 18)] <- 100
op <- par(mfcol = 1:2)
plot(tr); axisPhylo()
plotBreakLongEdges(tr, 2); axisPhylo()

## from ?boot.phylo:
f <- function(x) nj(dist.dna(x))
data(woodmouse)
tw <- f(woodmouse) # NJ tree with K80 distance
set.seed(1)
## bootstrap with 100 replications:
(bp <- boot.phylo(tw, woodmouse, f, quiet = TRUE))
## the first value relates to the root node and is always 100
## it is ignored below:
plot(tw, "u")
drawSupportOnEdges(bp)
## more readable but the tree is really unrooted:
plot(tw)
drawSupportOnEdges(bp)
par(op)
}
\keyword{hplot}
