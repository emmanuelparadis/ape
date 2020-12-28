\name{add.scale.bar}
\alias{add.scale.bar}
\title{Add a Scale Bar to a Phylogeny Plot}
\usage{
add.scale.bar(x, y, length = NULL, ask = FALSE,
              lwd = 1, lcol = "black", ...)
}
\arguments{
  \item{x}{x location of the bar (can be left missing).}
  \item{y}{y location of the bar (can be left missing).}
  \item{length}{a numeric value giving the length of the scale bar. If
    none is supplied, a value is calculated from the data.}
  \item{ask}{a logical; if \code{TRUE} the user is asked to click where
    to draw the bar. The default is \code{FALSE}.}
  \item{lwd}{the width of the bar.}
  \item{lcol}{the colour of the bar (use \code{col} for the colour of
    the text).}
  \item{\dots}{further arguments to be passed to \code{text}.}
}
\description{
  This function adds a horizontal bar giving the scale of the branch
  lengths to a plot of a phylogenetic tree on the current graphical
  device.
}
\details{
  By default, the bar is placed in a corner of the graph depending on
  the direction of the tree. Otherwise both \code{x} and \code{y} must
  be specified (if only one is given it is ignored).

  The further arguments (\code{\dots}) are used to format the text. They
  may be \code{font}, \code{cex}, \code{col}, and so on (see examples
  below, and the help page on \code{\link[graphics]{text}}).

  The function \code{\link[graphics]{locator}}  may be used to
  determine the \code{x} and \code{y} arguments.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{axisPhylo}},
  \code{\link[graphics]{locator}}
}
\examples{
tr <- rtree(10)
layout(matrix(1:2, 2, 1))
plot(tr)
add.scale.bar()
plot(tr)
add.scale.bar(cex = 0.7, font = 2, col = "red")
layout(1)
}
\keyword{aplot}
