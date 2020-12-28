\name{evonet}
\alias{evonet}
\alias{as.evonet}
\alias{as.evonet.phylo}
\alias{read.evonet}
\alias{write.evonet}
\alias{print.evonet}
\alias{plot.evonet}
\alias{reorder.evonet}
\alias{as.phylo.evonet}
\alias{as.networx.evonet}
\alias{as.network.evonet}
\alias{as.igraph.evonet}
\alias{Nedge.evonet}
\title{Evolutionary Networks}
\description{
  \code{evonet} builds a network from a tree of class
  \code{"phylo"}. There are \code{print}, \code{plot}, and
  \code{reorder} methods as well as a few conversion functions.
}
\usage{
evonet(phy, from, to = NULL)
\method{print}{evonet}(x, ...)
\method{plot}{evonet}(x, col = "blue", lty = 1, lwd = 1, alpha = 0.5,
              arrows = 0, arrow.type = "classical", ...)
\method{Nedge}{evonet}(phy)
\method{reorder}{evonet}(x, order = "cladewise", index.only = FALSE, ...)
\method{as.phylo}{evonet}(x, ...)
\method{as.networx}{evonet}(x, weight = NA, ...)
\method{as.network}{evonet}(x, directed = TRUE, ...)
\method{as.igraph}{evonet}(x, directed = TRUE, use.labels = TRUE, ...)

as.evonet(x, ...)
\method{as.evonet}{phylo}(x, ...)

read.evonet(file = "", text = NULL, comment.char = "", ...)
write.evonet(x, file = "", ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{x}{an object of class \code{"evonet"}.}
  \item{from}{a vector (or a matrix if \code{to = NULL}) giving the node
    or tip numbers involved in the reticulations.}
  \item{to}{a vector of the same length than \code{from}.}
  \item{col, lty, lwd}{colors, line type and width of the reticulations
    (recycled if necessary).}
  \item{alpha}{a value between 0 and 1 specifying the transparency of
    the reticulations.}
  \item{arrows, arrow.type}{see \code{\link{fancyarrows}}.}
  \item{order, index.only}{see \code{\link{reorder.phylo}}.}
  \item{weight}{a numeric vector giving the weights for the
    reticulations when converting to the class \code{"networx"}
    (recycled or shortened if needed).}
  \item{directed}{a logical: should the network be considered as
    directed? \code{TRUE} by default.}
  \item{use.labels}{a logical specifying whether to use the tip and node
    labels when building the network of class \code{"igraph"}.}
  \item{file, text, comment.char}{see \code{\link{read.tree}}.}
  \item{\dots}{arguments passed to other methods.}
}
\details{
  \code{evonet} is a constructor function that checks the arguments.

  The classes \code{"networx"}, \code{"network"}, and \code{"igraph"}
  are defined in the packages \pkg{phangorn}, \pkg{network}, and
  \pkg{igraph}, respectively.

  \code{read.evonet} reads networks from files in extended newick format
  (Cardona et al. 2008).
}
\value{
  an object of class \code{c("evonet", "phylo")} which is made of an
  object of class \code{"\link{phylo}"} plus an element
  \code{reticulation} coding additional edges among nodes and uses the
  same coding rules than the \code{edge} matrix.

  The conversion functions return an object of the appropriate class.
}
\author{Emmanuel Paradis, Klaus Schliep}
\seealso{
  \code{\link[phangorn]{as.networx}} in package \pkg{phangorn}
}
\references{
  Cardona, G., Rossell, F., and Valiente, G. (2008) Extended Newick: it
  is time for a standard representation of phylogenetic
  networks. \emph{BMC Bioinformatics}, \bold{9}, 532.
}
\examples{
tr <- rcoal(5)
(x <- evonet(tr, 6:7, 8:9))
plot(x)
## simple example of extended Newick format:
(enet <- read.evonet(text = "((a:2,(b:1)#H1:1):1,(#H1,c:1):2);"))
plot(enet, arrows=1)
## from Fig. 2 in Cardona et al. 2008:
z <- read.evonet(text =
"((1,((2,(3,(4)Y#H1)g)e,(((Y#H1, 5)h,6)f)X#H2)c)a,((X#H2,7)d,8)b)r;")
z
plot(z)
\dontrun{
if (require(igraph)) {
    plot(as.igraph(z))
}}}
\keyword{manip}
\keyword{hplot}
