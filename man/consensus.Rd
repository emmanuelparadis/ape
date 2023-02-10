\name{consensus}
\alias{consensus}
\title{Concensus Trees}
\usage{
consensus(..., p = 1, check.labels = TRUE, rooted = FALSE)
}
\arguments{
  \item{\dots}{either (i) a single object of class \code{"phylo"}, (ii) a
    series of such objects separated by commas, or (iii) a list
    containing such objects.}
  \item{p}{a numeric value between 0.5 and 1 giving the proportion for a
    clade to be represented in the consensus tree.}
  \item{check.labels}{a logical specifying whether to check the labels
    of each tree. If \code{FALSE} (the default), it is assumed that all
    trees have the same tip labels, and that they are in the same order
    (see details).}
  \item{rooted}{a logical specifying whether the trees should be treated as
    rooted or not.}
}
\description{
  Given a series of trees, this function returns the consensus tree. By
  default, the strict-consensus tree is computed. To get the
  majority-rule consensus tree, use \code{p = 0.5}. Any value between
  0.5 and 1 can be used.
}
\details{
  Using \code{check.labels = FALSE} results in
  considerable decrease in computing times. This requires that all
  trees have the same tip labels, \emph{and} these labels are
  ordered similarly in all trees (in other words, the element
  \code{tip.label} are identical in all trees).

  Until \pkg{ape} 5.6-2, the trees passed to this function were
  implicitly treated as rooted, even when the option \code{rooted =
  FALSE} was used. This is now fixed (see PR65 on GitHub) so that, by
  default, the trees are explicitly treated as unrooted (even if
  \code{\link{is.rooted}} returns \code{TRUE}). Thus, it could
  be that results now differ from previous analyses (setting
  \code{rooted = TRUE} might help to replicate previous results).
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  Felsenstein, J. (2004) \emph{Inferring Phylogenies}. Sunderland:
  Sinauer Associates.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{prop.part}}, \code{\link{dist.topo}}
}
\keyword{manip}
