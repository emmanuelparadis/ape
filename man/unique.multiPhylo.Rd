\name{unique.multiPhylo}
\alias{unique.multiPhylo}
\title{Revomes Duplicate Trees}
\description{
  This function scans a list of trees, and returns a list with the
  duplicate trees removed. By default the labelled topologies are
  compared.
}
\usage{
\method{unique}{multiPhylo}(x, incomparables = FALSE,
        use.edge.length = FALSE,
        use.tip.label = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{"multiPhylo"}.}
  \item{incomparables}{unused (for compatibility with the generic).}
  \item{use.edge.length}{a logical specifying whether to consider the edge
    lengths in the comparisons; the default is \code{FALSE}.}
  \item{use.tip.label}{a logical specifying whether to consider the tip
    labels in the comparisons; the default is \code{TRUE}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  an object of class \code{"multiPhylo"} with an attribute
  \code{"old.index"} indicating which trees of the original list are
  similar (the tree of smaller index is taken as reference).
}
\author{Emmanuel Paradis}
\seealso{
  \code{all.equal.phylo}, \code{\link[base]{unique}} for the generic R
  function, \code{read.tree}, \code{read.nexus}
}
\examples{
TR <- rmtree(50, 4)
length(unique(TR)) # not always 15...
howmanytrees(4)
}
\keyword{manip}
