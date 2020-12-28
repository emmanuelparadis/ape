\name{rotate}
\alias{rotate}
\alias{rotateConstr}
\title{Swapping Sister Clades}
\description{
  For a given node, \code{rotate} exchanges the position of two clades
  descending from this node. It can handle dichotomies as well as
  polytomies. In the latter case, two clades from the polytomy are
  selected for swapping.

  \code{rotateConstr} rotates internal branches giving a constraint on
  the order of the tips.
}
\usage{
rotate(phy, node, polytom = c(1, 2))
rotateConstr(phy, constraint)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{node}{a vector of mode numeric or character specifying the
    number of the node.}
  \item{polytom}{a vector of mode numeric and length two specifying the
    two clades that should be exchanged in a polytomy.}
  \item{constraint}{a vector of mode character specifying the order of
    the tips as they should appear when plotting the tree (from bottom
    to top).}
}
\details{
  \code{phy} can be either rooted or unrooted, contain polytomies and lack
  branch lengths. In the presence of very short branch lengths it is
  convenient to plot the phylogenetic tree without branch lengths in order
  to identify the number of the node in question.

  \code{node} can be any of the interior nodes of a phylogenetic tree
  including the root node. Number of the nodes can be identified by the
  nodelabels function. Alternatively, you can specify a vector of length
  two that contains either the number or the names of two tips that
  coalesce in the node of interest.

  If the node subtends a polytomy, any two clades of the the polytomy
  can be chosen by polytom. On a plotted phylogeny, the clades are
  numbered from bottom to top and polytom is used to index the two
  clades one likes to swop.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Christoph Heibl \email{heibl@lmu.de}, Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{nodelabels}},
  \code{\link{root}}, \code{\link{drop.tip}}}
\examples{
# create a random tree:
tre <- rtree(25)

# visualize labels of internal nodes:
plot(tre, use.edge.length=FALSE)
nodelabels()

# rotate clades around node 30:
tre.new <- rotate(tre, 30)

# compare the results:
par(mfrow=c(1,2)) # split graphical device
plot(tre) # plot old tre
plot(tre.new) # plot new tree

# visualize labels of terminal nodes:
plot(tre)
tiplabels()

# rotate clades containing nodes 12 and 20:
tre.new <- rotate(tre, c(12, 21))

# compare the results:
par(mfrow=c(1,2)) # split graphical device
plot(tre) # plot old tre
plot(tre.new) # plot new tree

# or you migth just specify tiplabel names:
tre.new <- rotate(tre, c("t3", "t14"))

# compare the results:
par(mfrow=c(1,2)) # devide graphical device
plot(tre) # plot old tre
plot(tre.new) # plot new tree

# a simple example for rotateConstr:
A <- read.tree(text = "((A,B),(C,D));")
B <- read.tree(text = "(((D,C),B),A);")
B <- rotateConstr(B, A$tip.label)
plot(A); plot(B, d = "l")

# something more interesting (from ?cophyloplot):
tr1 <- rtree(40)
## drop 20 randomly chosen tips:
tr2 <- drop.tip(tr1, sample(tr1$tip.label, size = 20))
## rotate the root and reorder the whole:
tr2 <- rotate(tr2, 21)
tr2 <- read.tree(text = write.tree(tr2))
X <- cbind(tr2$tip.label, tr2$tip.label) # association matrix
cophyloplot(tr1, tr2, assoc = X, space = 28)
## before reordering tr2 we have to find the constraint:
co <- tr2$tip.label[order(match(tr2$tip.label, tr1$tip.label))]
newtr2 <- rotateConstr(tr2, co)
cophyloplot(tr1, newtr2, assoc = X, space = 28)
}
\keyword{manip}
