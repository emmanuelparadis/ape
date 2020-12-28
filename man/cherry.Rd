\name{cherry}
\alias{cherry}
\title{Number of Cherries and Null Models of Trees}
\usage{
cherry(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
}
\description{
  This function calculates the number of cherries (see definition below)
  on a phylogenetic tree, and tests the null hypotheses whether this
  number agrees with those predicted from two null models of trees (the
  Yule model, and the uniform model).
}
\value{
  A NULL value is returned, the results are simply printed.
}
\details{
  A cherry is a pair of adjacent tips on a tree. The tree can be either
  rooted or unrooted, but the present function considers only rooted
  trees. The probability distribution function of the number of cherries
  on a tree depends on the speciation/extinction model that generated
  the tree.

  McKenzie and Steel (2000) derived the probability
  distribution function of the number of cherries for two models: the
  Yule model and the uniform model. Broadly, in the Yule model, each extant
  species is equally likely to split into two daughter-species; in the
  uniform model, a branch is added to tree on any of the already
  existing branches with a uniform probability.

  The probabilities are computed using recursive formulae; however, for
  both models, the probability density function converges to a normal
  law with increasing number of tips in the tree. The function uses
  these normal approximations for a number of tips greater than or equal
  to 20.
}
\references{
  McKenzie, A. and Steel, M. (2000) Distributions of cherries for two
  models of trees. \emph{Mathematical Biosciences}, \bold{164}, 81--92.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{gammaStat}}
}
\keyword{univar}
