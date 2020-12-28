\name{which.edge}
\alias{which.edge}
\title{Identifies Edges of a Tree}
\description{
  This function identifies the edges that belong to a group (possibly
  non-monophyletic) specified as a set of tips.
}
\usage{
which.edge(phy, group)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{group}{a vector of mode numeric or character specifying the tips
    for which the edges are to be identified.}
}
\details{
  The group of tips specified in `group' may be non-monophyletic
  (paraphyletic or polyphyletic), in which case all edges from the tips
  to their most recent common ancestor are identified.

  The identification is made with the indices of the rows of the matrix
  `edge' of the tree.
}
\value{
  a numeric vector.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{bind.tree}}, \code{\link{drop.tip}}, \code{\link{root}}
}
\keyword{manip}
