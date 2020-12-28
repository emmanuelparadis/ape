\name{identify.phylo}
\alias{identify.phylo}
\title{Graphical Identification of Nodes and Tips}
\usage{
\method{identify}{phylo}(x, nodes = TRUE, tips = FALSE,
                  labels = FALSE, quiet = FALSE, ...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{nodes}{a logical specifying whether to identify the node.}
  \item{tips}{a logical specifying whether to return the tip
    information.}
  \item{labels}{a logical specifying whether to return the labels; by
    default only the numbers are returned.}
  \item{quiet}{a logical controlling whether to print a message inviting
    the user to click on the tree.}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\description{
  This function allows to identify a clade on a plotted tree by clicking
  on the plot with the mouse. The tree, specified in the argument
  \code{x}, must be plotted beforehand.
}
\details{
  By default, the clade is identified by its number as found in the
  `edge' matrix of the tree. If \code{tips = TRUE}, the tips descending
  from the identified node are returned, possibly together with the
  node. If \code{labels = TRUE}, the labels are returned (if the tree
  has no node labels, then the node numbered is returned).

  The node is identified by the shortest distance where the click
  occurs. If the click occurs close to a tip, the function returns its
  information.
}
\note{
  This function does not add anything on the plot, but it can be wrapped
  with, e.g., \code{\link{nodelabels}} (see example), or its results can
  be sent to, e.g., \code{\link{drop.tip}}.
}
\value{
  A list with one or two vectors named \code{"tips"} and/or
  \code{"nodes"} with the identification of the tips and/or of the
  nodes.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{nodelabels}},
  \code{\link[graphics]{identify}} for the generic function
}
\examples{
\dontrun{
tr <- rtree(20)
f <- function(col) {
    o <- identify(tr)
    nodelabels(node=o$nodes, pch = 19, col = col)
}
plot(tr)
f("red") # click close to a node
f("green")
}
}
\keyword{aplot}
