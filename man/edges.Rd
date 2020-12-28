\name{edges}
\alias{edges}
\alias{fancyarrows}
\title{Draw Additional Edges on a Plotted Tree}
\description{
  \code{edges} draws edges on a plotted tree. \code{fancyarrows}
  enhances \code{\link[graphics]{arrows}} with triangle and harpoon
  heads; it can be called from \code{edges}.
}
\usage{
edges(nodes0, nodes1, arrows = 0, type = "classical", ...)
fancyarrows(x0, y0, x1, y1, length = 0.25, angle = 30, code = 2,
            col = par("fg"), lty = par("lty"), lwd = par("lwd"),
            type = "triangle", ...)
}
\arguments{
  \item{nodes0, nodes1}{vectors of integers giving the tip and/or node
    numbers where to start and to end the edges (eventually recycled).}
  \item{arrows}{an integer between 0 and 3; 0: lines (the default); 1:
    an arrow head is drawn at \code{nodes0}; 2: at \code{nodes1}; 3:
    both.}
  \item{type}{if the previous argument is not 0, the type of arrow head:
    \code{"classical"} (just lines, the default), \code{"triangle"},
    \code{"harpoon"}, or any unambiguous abbreviations of these. For
    \code{fancyarrows} only the last two are available.}
  \item{x0, y0, x1, y1}{the coordinates of the start and end points for
    \code{fancyarrows} (these are not recycled and so should be vectors
    of the same length).}
  \item{length, angle, code, col, lty, lwd}{default options similar to
    those of \code{\link[graphics]{arrows}}.}
  \item{\dots}{further arguments passed to \code{\link[graphics]{segments}}.}
}
\details{
  The first function is helpful when drawing reticulations on a phylogeny,
  especially if computed from the edge matrix.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{nodelabels}}
}
\examples{
set.seed(2)
tr <- rcoal(6)
plot(tr, "c")
edges(10, 9, col = "red", lty = 2)
edges(10:11, 8, col = c("blue", "green")) # recycling of 'nodes1'
edges(1, 2, lwd = 2, type = "h", arrows = 3, col = "green")
nodelabels()
}
\keyword{aplot}
