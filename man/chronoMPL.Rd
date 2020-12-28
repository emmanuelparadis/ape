\name{chronoMPL}
\alias{chronoMPL}
\title{Molecular Dating With Mean Path Lengths}
\usage{
chronoMPL(phy, se = TRUE, test = TRUE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{se}{a logical specifying whether to compute the standard-errors
    of the node ages (\code{TRUE} by default).}
  \item{test}{a logical specifying whether to test the molecular clock
    at each node (\code{TRUE} by default).}
}
\description{
  This function estimates the node ages of a tree using the mean path
  lengths method of Britton et al. (2002). The branch lengths of the
  input tree are interpreted as (mean) numbers of substitutions.
}
\details{
  The mean path lengths (MPL) method estimates the age of a node with
  the mean of the distances from this node to all tips descending from
  it. Under the assumption of a molecular clock, standard-errors of the
  estimates node ages can be computed (Britton et al. 2002).

  The tests performed if \code{test = TRUE} is a comparison of the MPL
  of the two subtrees originating from a node; the null hypothesis is
  that the rate of substitution was the same in both subtrees (Britton
  et al. 2002). The test statistic follows, under the null hypothesis, a
  standard normal distribution. The returned \emph{P}-value is the
  probability of observing a greater absolute value (i.e., a two-sided
  test). No correction for multiple testing is applied: this is left to
  the user.

  Absolute dating can be done by multiplying the edge lengths found by
  calibrating one node age.
}
\note{
  The present version requires a dichotomous tree.
}
\value{
  an object of class \code{"phylo"} with branch lengths as estimated by
  the function. There are, by default, two attributes:

  \item{stderr}{the standard-errors of the node ages.}
  \item{Pval}{the \emph{P}-value of the test of the molecular clock for
    each node.}
}
\references{
  Britton, T., Oxelman, B., Vinnersten, A. and Bremer, K. (2002)
  Phylogenetic dating with confidence intervals using mean path
  lengths. \emph{Molecular Phylogenetics and Evolution}, \bold{24},
  58--65.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{chronopl}}
}
\examples{
tr <- rtree(10)
tr$edge.length <- 5*tr$edge.length
chr <- chronoMPL(tr)
layout(matrix(1:4, 2, 2, byrow = TRUE))
plot(tr)
title("The original tree")
plot(chr)
axisPhylo()
title("The dated MPL tree")
plot(chr)
nodelabels(round(attr(chr, "stderr"), 3))
title("The standard-errors")
plot(tr)
nodelabels(round(attr(chr, "Pval"), 3))
title("The tests")
layout(1)
}
\keyword{models}
