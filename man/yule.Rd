\name{yule}
\alias{yule}
\title{Fits the Yule Model to a Phylogenetic Tree}
\usage{
yule(phy, use.root.edge = FALSE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{use.root.edge}{a logical specifying whether to consider the root
    edge in the calculations.}
}
\description{
  This function fits by maximum likelihood a Yule model, i.e., a
  birth-only model to the branching times computed from a phylogenetic
  tree.
}
\details{
  The tree must be fully dichotomous.

  The maximum likelihood estimate of the speciation rate is obtained by
  the ratio of the number of speciation events on the cumulative number
  of species through time; these two quantities are obtained with the
  number of nodes in the tree, and the sum of the branch lengths,
  respectively.

  If there is a `root.edge' element in the phylogenetic tree, and
  \code{use.root.edge = TRUE}, then it is assumed that it has a
  biological meaning and is counted as a branch length, and the root is
  counted as a speciation event; otherwise the number of speciation
  events is the number of nodes - 1.

  The standard-error of lambda is computed with the second derivative of
  the log-likelihood function.
}
\value{
  An object of class "yule" which is a list with the following
  components:
  \item{lambda}{the maximum likelihood estimate of the speciation
    (birth) rate.}
  \item{se}{the standard-error of lambda.}
  \item{loglik}{the log-likelihood at its maximum.}
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{branching.times}}, \code{\link{diversi.gof}},
  \code{\link{diversi.time}}, \code{\link{ltt.plot}},
  \code{\link{birthdeath}}, \code{\link{bd.ext}}, \code{\link{yule.cov}}
}
\keyword{models}
