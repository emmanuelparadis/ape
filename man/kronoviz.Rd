\name{kronoviz}
\alias{kronoviz}
\title{Plot Multiple Chronograms on the Same Scale}
\description{
  The main argument is a list of (rooted) trees which are plotted on the
  same scale.
}
\usage{
kronoviz(x, layout = length(x), horiz = TRUE, ...)
}
\arguments{
  \item{x}{a list of (rooted) trees of class \code{"phylo"}.}
  \item{layout}{an integer giving the number of trees plotted
    simultaneously; by default all.}
  \item{horiz}{a logical specifying whether the trees should be plotted
    rightwards (the default) or upwards.}
  \item{\dots}{further arguments passed to \code{plot.phylo}.}
}
\details{
  The size of the individual plots is proportional to the size of the
  trees.
}
\value{NULL}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}
}
\examples{
TR <- replicate(10, rcoal(sample(11:20, size = 1)), simplify = FALSE)
kronoviz(TR)
kronoviz(TR, horiz = FALSE, type = "c", show.tip.label = FALSE)
}
\keyword{hplot}
