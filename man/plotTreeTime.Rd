\name{plotTreeTime}
\alias{plotTreeTime}
\title{Plot Tree With Time Axis}
\description{
  This function plots a non-ultrametric tree where the tips are not
  contemporary together with their dates on the x-axis.
}
\usage{
plotTreeTime(phy, tip.dates, show.tip.label = FALSE, y.lim = NULL,
             color = TRUE, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{tip.dates}{a vector of the same length than the number of tips
    in \code{phy} (see details).}
  \item{show.tip.label}{a logical value; see \code{\link{plot.phylo}}.}
  \item{y.lim}{by default, one fifth of the plot is left below the tree;
    use this option to change this behaviour.}
  \item{color}{a logical value specifying whether to use colors for the
    lines linking the tips to the time axis. If \code{FALSE}, a grey
    scale is used.}
  \item{\dots}{other arguments to be passed to \code{plot.phylo}.}
}
\details{
  The vector \code{tip.dates} may be numeric or of class
  \dQuote{\link[base]{Date}}. In either case, the time axis is set
  accordingly. The length of this vector must be equal to the number of
  tips of the tree: the dates are matched to the tips numbers. Missing
  values are allowed.
}
\value{NULL}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{estimate.dates}}
}
\examples{
dates <- as.Date(.leap.seconds)
tr <- rtree(length(dates))
plotTreeTime(tr, dates)

## handling NA's:
dates[11:26] <- NA
plotTreeTime(tr, dates)

## dates can be on an arbitrary scale, e.g., [-1, 1]:
plotTreeTime(tr, runif(Ntip(tr), -1, 1))
}
\keyword{hplot}
