\name{as.matching}
\alias{as.matching}
\alias{matching}
\alias{as.matching.phylo}
\alias{as.phylo.matching}
\title{Conversion Between Phylo and Matching Objects}
\description{
  These functions convert objects between the classes \code{"phylo"} and
  \code{"matching"}.
}
\usage{
as.matching(x, ...)
\method{as.matching}{phylo}(x, labels = TRUE, ...)
\method{as.phylo}{matching}(x, ...)
}
\arguments{
  \item{x}{an object to convert as an object of class \code{"matching"}
    or of class \code{"phylo"}.}
  \item{labels}{a logical specifying whether the tip and node labels
    should be included in the returned matching.}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
  A matching is a representation where each tip and each node are given
  a number, and sibling groups are grouped in a ``matching pair'' (see
  Diaconis and Holmes 1998, for details). This coding system can be used
  only for binary (fully dichotomous) trees.

  Diaconis and Holmes (1998) gave some conventions to insure that a
  given tree has a unique representation as a matching. I have tried to
  follow them in the present functions.
}
\value{
  \code{as.matching} returns an object of class \code{"matching"} with
  the following component:

  \item{matching}{a two-column numeric matrix where the columns
    represent the sibling pairs.}
  \item{tip.label}{(optional) a character vector giving the tip labels
    where the ith element is the label of the tip numbered i in
    \code{matching}.}
  \item{node.label}{(optional) a character vector giving the node
    labels in the same order than in \code{matching} (i.e. the ith
    element is the label of the node numbered i + n in \code{matching},
    with n the number of tips).}

  \code{as.phylo.matching} returns an object of class \code{"phylo"}.
}
\note{
  Branch lengths are not supported in the present version.
}
\author{Emmanuel Paradis}
\references{
  Diaconis, P. W. and Holmes, S. P. (1998) Matchings and phylogenetic
  trees. \emph{Proceedings of the National Academy of Sciences USA},
  \bold{95}, 14600--14602.
}
\seealso{\code{\link{as.phylo}}}
\examples{
data(bird.orders)
m <- as.matching(bird.orders)
str(m)
m
tr <- as.phylo(m)
all.equal(tr, bird.orders, use.edge.length = FALSE)
}
\keyword{manip}
