\name{as.bitsplits}
\alias{as.bitsplits}
\alias{as.bitsplits.prop.part}
\alias{print.bitsplits}
\alias{sort.bitsplits}
\alias{bitsplits}
\alias{countBipartitions}
\alias{as.prop.part}
\alias{as.prop.part.bitsplits}
\title{Split Frequencies and Conversion Among Split Classes}
\description{
  \code{bitsplits} returns the bipartitions (aka splits) for a single
  tree or a list of trees.

  \code{countBipartitions} returns the frequencies of the bipartitions
  from a reference tree (phy) observed in a list of trees (X).

  \code{as.bitsplits} and \code{as.prop.part} are generic functions for
  converting between the \code{"bitsplits"} and \code{"prop.part"}
  classes.
}
\usage{
bitsplits(x)
countBipartitions(phy, X)
as.bitsplits(x)
\method{as.bitsplits}{prop.part}(x)
\method{print}{bitsplits}(x, ...)
\method{sort}{bitsplits}(x, decreasing = FALSE, ...)
as.prop.part(x, ...)
\method{as.prop.part}{bitsplits}(x, include.trivial = FALSE, ...)
}
\arguments{
  \item{x}{an object of the appropriate class.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{X}{an object of class \code{"multiPhylo"}.}
  \item{decreasing}{a logical value to sort the bipartitions in
    increasing (the default) or decreasing order of their frequency.}
  \item{include.trivial}{a logical value specifying whether to include
    the trivial split with all tips in the returned object.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  These functions count bipartitions as defined by internal branches, so
  they do not work with rooted trees (see examples). The structure of
  the class \code{"bitsplits"} is described in a separate document
  on ape's web site.
}
\value{
  \code{bitsplits}, \code{as.bitsplits}, and \code{sort} return an object
  of class \code{"bitsplits"}.

  \code{countBipartitions} returns a vector of integers.

  \code{as.prop.part} returns an object of class \code{"prop.part"}.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{prop.part}}, \code{\link{is.compatible}}}
\examples{
tr <- rtree(20)
pp <- prop.part(tr)
as.bitsplits(pp)
## doesn't work for rooted trees...:
countBipartitions(rtree(10), rmtree(100, 10))
## ... but OK with unrooted trees:
countBipartitions(rtree(10, rooted = FALSE), rmtree(100, 10, rooted = FALSE))
}
\keyword{manip}
