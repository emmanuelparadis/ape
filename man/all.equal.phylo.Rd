\encoding{utf8}
\name{all.equal.phylo}
\alias{all.equal.phylo}
\title{Global Comparison of two Phylogenies}
\usage{
\method{all.equal}{phylo}(target, current, use.edge.length = TRUE,
                   use.tip.label = TRUE, index.return = FALSE,
                   tolerance = .Machine$double.eps ^ 0.5,
                   scale = NULL, \dots)
}
\arguments{
  \item{target}{an object of class \code{"phylo"}.}
  \item{current}{an object of class \code{"phylo"}.}
  \item{use.edge.length}{if \code{FALSE} only the topologies are
    compared; the default is \code{TRUE}.}
  \item{use.tip.label}{if \code{FALSE} the unlabelled trees are
    compared; the default is \code{TRUE}.}
  \item{index.return}{if \code{TRUE} the function returns a two-column
    matrix giving the correspondence between the nodes of both trees.}
  \item{tolerance}{the numeric tolerance used to compare the branch
    lengths.}
  \item{scale}{a positive number, comparison of branch lengths is made
    after scaling (i.e., dividing) them by this number.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  This function makes a global comparison of two phylogenetic trees.
}
\details{
  This function is meant to be an adaptation of the generic function
  \code{all.equal} for the comparison of phylogenetic trees.

  A single phylogenetic tree may have several representations in the Newick
  format and in the \code{"phylo"} class of objects used in `ape'. One
  aim of the present function is to be able to identify whether two
  objects of class \code{"phylo"} represent the same phylogeny.
}
\note{
  The algorithm used here does not work correctly for the comparison of
  topologies (i.e., ignoring tip labels) of unrooted trees. This also
  affects \code{\link{unique.multiPhylo}} which calls the present function. See:

  \url{https://www.mail-archive.com/r-sig-phylo@r-project.org/msg01445.html}.
}
\value{
  A logical value, or a two-column matrix.
}
\author{\enc{Benoît}{Benoit} Durand \email{b.durand@alfort.AFSSA.FR}}
\seealso{
  \code{\link[base]{all.equal}} for the generic \R function, \code{\link{comparePhylo}}
}
\examples{
### maybe the simplest example of two representations
### for the same rooted tree...:
t1 <- read.tree(text = "(a:1,b:1);")
t2 <- read.tree(text = "(b:1,a:1);")
all.equal(t1, t2)
### ... compare with this:
identical(t1, t2)
### one just slightly more complicated...:
t3 <- read.tree(text = "((a:1,b:1):1,c:2);")
t4 <- read.tree(text = "(c:2,(a:1,b:1):1);")
all.equal(t3, t4) # == all.equal.phylo(t3, t4)
### ... here we force the comparison as lists:
all.equal.list(t3, t4)
}
\keyword{manip}
