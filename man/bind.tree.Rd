\name{bind.tree}
\alias{bind.tree}
\alias{+.phylo}
\title{Binds Trees}
\usage{
bind.tree(x, y, where = "root", position = 0, interactive = FALSE)
\special{x + y}
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{y}{an object of class \code{"phylo"}.}
  \item{where}{an integer giving the number of the node or tip of the
    tree \code{x} where the tree \code{y} is binded (\code{"root"} is a
    short-cut for the root).}
  \item{position}{a numeric value giving the position from the tip or
    node given by \code{node} where the tree \code{y} is binded;
    negative values are ignored.}
  \item{interactive}{if \code{TRUE} the user is asked to choose the tip
    or node of \code{x} by clicking on the tree which must be plotted.}
}
\description{
  This function binds together two phylogenetic trees to give a single
  object of class \code{"phylo"}.
}
\details{
  The argument \code{x} can be seen as the receptor tree, whereas
  \code{y} is the donor tree. The root of \code{y} is then grafted on a
  location of \code{x} specified by \code{where} and, possibly,
  \code{position}. If \code{y} has a root edge, this is added as in
  internal branch in the resulting tree.

  \code{x + y} is a shortcut for:

  \preformatted{
    bind.tree(x, y, position = if (is.null(x$root.edge)) 0 else
    x$root.edge)
  }

  If only one of the trees has no branch length, the branch lengths of
  the other one are ignored with a warning.

  If one (or both) of the trees has no branch length, it is possible to
  specify a value of 'position' to graft 'y' below the node of 'x'
  specified by 'where'. In this case, the exact value of 'position' is
  not important as long as it is greater than zero. The new node will be
  multichotomous if 'y' has no root edge. This can be solved by giving
  an arbitrary root edge to 'y' beforehand (e.g., \code{y$root.edge <-
  1}): it will be deleted during the binding operation.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{drop.tip}}, \code{\link{root}}
}
\examples{
### binds the two clades of bird orders
cat("((Struthioniformes:21.8,Tinamiformes:21.8):4.1,",
    "((Craciformes:21.6,Galliformes:21.6):1.3,Anseriformes:22.9):3.0):2.1;",
    file = "ex1.tre", sep = "\n")
cat("(Turniciformes:27.0,(Piciformes:26.3,((Galbuliformes:24.4,",
    "((Bucerotiformes:20.8,Upupiformes:20.8):2.6,",
    "(Trogoniformes:22.1,Coraciiformes:22.1):1.3):1.0):0.6,",
    "(Coliiformes:24.5,(Cuculiformes:23.7,(Psittaciformes:23.1,",
    "(((Apodiformes:21.3,Trochiliformes:21.3):0.6,",
    "(Musophagiformes:20.4,Strigiformes:20.4):1.5):0.6,",
    "((Columbiformes:20.8,(Gruiformes:20.1,Ciconiiformes:20.1):0.7):0.8,",
    "Passeriformes:21.6):0.9):0.6):0.6):0.8):0.5):1.3):0.7):1.0;",
    file = "ex2.tre", sep = "\n")
tree.bird1 <- read.tree("ex1.tre")
tree.bird2 <- read.tree("ex2.tre")
unlink(c("ex1.tre", "ex2.tre")) # clean-up
(birds <- tree.bird1 + tree.bird2)
layout(matrix(c(1, 2, 3, 3), 2, 2))
plot(tree.bird1)
plot(tree.bird2)
plot(birds)

### examples with random trees
x <- rtree(4, tip.label = LETTERS[1:4])
y <- rtree(4, tip.label = LETTERS[5:8])
x <- makeNodeLabel(x, prefix = "x_")
y <- makeNodeLabel(y, prefix = "y_")
x$root.edge <- y$root.edge <- .2

z <- bind.tree(x, y, po=.2)
plot(y, show.node.label = TRUE, font = 1, root.edge = TRUE)
title("y")
plot(x, show.node.label = TRUE, font = 1, root.edge = TRUE)
title("x")
plot(z, show.node.label = TRUE, font = 1, root.edge = TRUE)
title("z <- bind.tree(x, y, po=.2)")

## make sure the terminal branch length is long enough:
x$edge.length[x$edge[, 2] == 2] <- 0.2

z <- bind.tree(x, y, 2, .1)
plot(y, show.node.label = TRUE, font = 1, root.edge = TRUE)
title("y")
plot(x, show.node.label = TRUE, font = 1, root.edge = TRUE)
title("x")
plot(z, show.node.label = TRUE, font = 1, root.edge = TRUE)
title("z <- bind.tree(x, y, 2, .1)")

x <- rtree(50)
y <- rtree(50)
x$root.edge <- y$root.edge <- .2
z <- x + y
plot(y, show.tip.label = FALSE, root.edge = TRUE); axisPhylo()
title("y")
plot(x, show.tip.label = FALSE, root.edge = TRUE); axisPhylo()
title("x")
plot(z, show.tip.label = FALSE, root.edge = TRUE); axisPhylo()
title("z <- x + y")
layout(1)
}
\keyword{manip}
