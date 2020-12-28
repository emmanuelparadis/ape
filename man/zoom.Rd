\name{zoom}
\alias{zoom}
\title{Zoom on a Portion of a Phylogeny}
\description{
  This function plots simultaneously a whole phylogenetic tree
  (supposedly large) and a portion of it.
}
\usage{
zoom(phy, focus, subtree = FALSE, col = rainbow, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{focus}{a vector, either numeric or character, or a list of
    vectors specifying the tips to be focused on.}
  \item{subtree}{a logical indicating whether to show the context of the
    extracted subtrees.}
  \item{col}{a vector of colours used to show where the subtrees are in
    the main tree, or a function .}
  \item{\dots}{further arguments passed to \code{plot.phylo}.}
}
\details{
  This function aims at exploring very large trees. The main argument is
  a phylogenetic tree, and the second one is a vector or a list of
  vectors specifying the tips to be focused on. The vector(s) can be
  either numeric and thus taken as the indices of the tip labels, or
  character in which case it is taken as the corresponding tip labels.

  The whole tree is plotted on the left-hand side in a narrower
  sub-window (about a quarter of the device) without tip labels. The
  subtrees consisting of the tips in `focus' are extracted and plotted
  on the right-hand side starting from the top left corner and
  successively column-wise.

  If the argument `col' is a vector of colours, as many colours as the
  number of subtrees must be given. The alternative is to give a
  function that will create colours or grey levels from the number of
  subtrees: see \code{\link[grDevices]{rainbow}} for some possibilities
  with colours.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{drop.tip}},
  \code{\link[graphics]{layout}}, \code{\link[grDevices]{rainbow}},
  \code{\link[grDevices]{grey}}
}
\examples{
\dontrun{
data(chiroptera)
zoom(chiroptera, 1:20, subtree = TRUE)
zoom(chiroptera, grep("Plecotus", chiroptera$tip.label))
zoom(chiroptera, list(grep("Plecotus", chiroptera$tip.label),
                      grep("Pteropus", chiroptera$tip.label)))
}
}
\keyword{hplot}
