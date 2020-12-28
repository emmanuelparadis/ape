\name{stree}
\alias{stree}
\title{Generates Systematic Regular Trees}
\usage{
stree(n, type = "star", tip.label = NULL)
}
\arguments{
  \item{n}{an integer giving the number of tips in the tree.}
  \item{type}{a character string specifying the type of tree to
    generate; four choices are possible: \code{"star"},
    \code{"balanced"}, \code{"left"}, \code{"right"}, or any unambiguous
    abbreviation of these.}
  \item{tip.label}{a character vector giving the tip labels; if not
    specified, the tips "t1", "t2", ..., are given.}
}
\description{
  This function generates trees with regular shapes.
}
\details{
  The types of trees generated are:

  \itemize{
    \item{``star''}{a star (or comb) tree with a single internal node.}
    \item{``balanced''}{a fully balanced dichotomous rooted tree;
      \code{n} must be a power of 2 (2, 4, 8, \dots).}
    \item{``left''}{a fully unbalanced rooted tree where the largest
      clade is on the left-hand side when the tree is plotted upwards.}
    \item{``right''}{same than above but in the other direction.}
  }
}
\value{
  An object of class \code{"phylo"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{compute.brlen}}, \code{\link{rtree}}
}
\examples{
layout(matrix(1:4, 2, 2))
plot(stree(100))
plot(stree(128, "balanced"))
plot(stree(100, "left"))
plot(stree(100, "right"))
}
\keyword{datagen}
