\name{cophyloplot}
\alias{cophyloplot}
\title{Plots two phylogenetic trees face to face with links between the tips.}
\description{
  This function plots two trees face to face with the links if specified. It is possible to rotate the branches of each tree around the nodes by clicking.
}
\usage{
cophyloplot(x, y, assoc = NULL, use.edge.length = FALSE, space = 0,
       length.line = 1, gap = 2, type = "phylogram", rotate = FALSE,
       col = par("fg"), lwd = par("lwd"), lty = par("lty"),
       show.tip.label = TRUE, font = 3, \dots)
}

\arguments{
  \item{x, y}{two objects of class \code{"phylo"}.}
  \item{assoc}{a matrix with 2 columns specifying the associations
    between the tips. If NULL, no links will be drawn.}
  \item{use.edge.length}{a logical indicating whether the branch lengths
    should be used to plot the trees; default is FALSE.}
  \item{space}{a positive value that specifies the distance between the
    two trees.}
  \item{length.line}{a positive value that specifies the length of the
    horizontal line associated to each taxa. Default is 1.}
  \item{gap}{a value specifying the distance between the tips of the
    phylogeny and the lines.}
  \item{type}{a character string specifying the type of phylogeny to be
    drawn; it must be one of "phylogram" (the default) or "cladogram".}
  \item{rotate}{a logical indicating whether the nodes of the phylogeny
    can be rotated by clicking. Default is FALSE.}
  \item{col}{a character vector indicating the color to be used for the
    links; recycled as necessary.}
  \item{lwd}{id. for the width.}
  \item{lty}{id. for the line type.}
  \item{show.tip.label}{a logical indicating whether to show the tip
    labels on the phylogeny (defaults to 'TRUE', i.e. the labels are
    shown).}
  \item{font}{an integer specifying the type of font for the
    labels: 1 (plain text), 2 (bold), 3 (italic, the default), or 4
    (bold italic).}
  \item{\dots}{(unused)}
}
\details{
The aim of this function is to plot simultaneously two phylogenetic trees with associated taxa. The two trees do not necessarily have the same number of tips and more than one tip in one phylogeny can be associated with a tip in the other.

The association matrix used to draw the links has to be a matrix with two columns containing the names of the tips. One line in the matrix represents one link on the plot. The first column of the matrix has to contain tip labels of the first tree (\code{phy1}) and the second column of the matrix, tip labels of the second tree (\code{phy2}). There is no limit (low or high) for the number of lines in the matrix. A matrix with two colums and one line will give a plot with one link.

Arguments \code{gap}, \code{length.line} and \code{space} have to be changed to get a nice plot of the two phylogenies. Note that the function takes into account the length of the character strings corresponding to the names at the tips, so that the lines do not overwrite those names.

The \code{rotate} argument can be used to transform both phylogenies in order to get the more readable plot (typically by decreasing the number of crossing lines). This can be done by clicking on the nodes. The escape button or right click take back to the console.
}
\author{Damien de Vienne \email{damien.de-vienne@u-psud.fr}}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{rotate}}, \code{\link{rotateConstr}}
}
\examples{
#two random trees
tree1 <- rtree(40)
tree2 <- rtree(20)

#creation of the association matrix:
association <- cbind(tree2$tip.label, tree2$tip.label)

cophyloplot(tree1, tree2, assoc = association,
            length.line = 4, space = 28, gap = 3)

#plot with rotations
\dontrun{
cophyloplot(tree1, tree2, assoc=association, length.line=4, space=28, gap=3, rotate=TRUE)
}
}
\keyword{hplot}
