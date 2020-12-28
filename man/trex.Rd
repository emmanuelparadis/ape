\name{trex}
\alias{trex}
\title{Tree Explorer With Multiple Devices}
\description{
  This function requires a plotted tree: the user is invited to click
  close to a node and the corresponding subtree (or clade) is plotted on
  a new window.
}
\usage{
trex(phy, title = TRUE, subbg = "lightyellow3",
     return.tree = FALSE, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{title}{a logical or a character string (see details).}
  \item{subbg}{a character string giving the background colour for the
    subtree.}
  \item{return.tree}{a logical: if \code{TRUE}, the subtree is returned
    after being plotted and the operation is stopped.}
  \item{\dots}{further arguments to pass to \code{plot.phylo}.}
}
\details{
  This function works with a tree (freshly) plotted on an interactive
  graphical device (i.e., not a file). After calling \code{trex}, the
  user clicks close to a node of the tree, then the clade from this node
  is plotted on a \emph{new} window. The user can click as many times on
  the main tree: the clades are plotted successively on the \emph{same}
  new window. The process is stopped by a right-click. If the user clicks
  too close to the tips, a message ``Try again!'' is printed.

  Each time \code{trex} is called, the subtree is plotted on a new
  window without closing or deleting those possibly already
  plotted. They may be distinguished with the options \code{title}
  and/or \code{subbg}.

  In all cases, the device where \code{phy} is plotted is the active
  window after the operation. It should \emph{not} be closed during the
  whole process.

  If \code{title = TRUE}, a default title is printed on the new window
  using the node label, or the node number if there are no node labels
  in the tree. If \code{title = FALSE}, no title is printed. If
  \code{title} is a character string, it is used for the title.
}
\value{
  an object of class \code{"phylo"} if \code{return.tree = TRUE}
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{identify.phylo}}
}
\examples{
\dontrun{
tr <- rcoal(1000)
plot(tr, show.tip.label = FALSE)
trex(tr) # left-click as many times as you want, then right-click
tr <- makeNodeLabel(tr)
trex(tr, subbg = "lightgreen") # id.

## generate a random colour with control on the darkness:
rRGB <- function(a, b)
    rgb(runif(1, a, b), runif(1, a, b), runif(1, a, b))

### with a random pale background:
trex(tr, subbg = rRGB(0.8, 1))
## the above can be called many times...
graphics.off() # close all graphical devices
}}
\keyword{hplot}
