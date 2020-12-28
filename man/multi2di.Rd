\name{multi2di}
\alias{multi2di}
\alias{multi2di.phylo}
\alias{multi2di.multiPhylo}
\alias{di2multi}
\alias{di2multi.phylo}
\alias{di2multi.multiPhylo}
\title{Collapse and Resolve Multichotomies}
\description{
  These two functions collapse or resolve multichotomies in phylogenetic
  trees.
}
\usage{
multi2di(phy, ...)
\method{multi2di}{phylo}(phy, random = TRUE, equiprob = TRUE, ...)
\method{multi2di}{multiPhylo}(phy, random = TRUE, equiprob = TRUE, ...)
di2multi(phy, ...)
\method{di2multi}{phylo}(phy, tol = 1e-08, ...)
\method{di2multi}{multiPhylo}(phy, tol = 1e-08, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{random}{a logical value specifying whether to resolve the
    multichotomies randomly (the default) or in the order they appear in
    the tree (if \code{random = FALSE}).}
  \item{equiprob}{a logical value: should topologies generated in equal
    probabilities; see details in \code{\link{rtree}} (ignored if
    \code{random = FALSE}).}
  \item{tol}{a numeric value giving the tolerance to consider a branch
    length significantly greater than zero.}
  \item{\dots}{arguments passed among methods.}
}
\details{
  \code{multi2di} transforms all multichotomies into a series of
  dichotomies with one (or several) branch(es) of length zero.

  \code{di2multi} deletes all branches smaller than \code{tol} and
  collapses the corresponding dichotomies into a multichotomy.
}
\seealso{\code{\link{is.binary}}}
\author{Emmanuel Paradis}
\value{
  an object of the same class than the input.
}
\examples{
data(bird.families)
is.binary(bird.families)
is.binary(multi2di(bird.families))
all.equal(di2multi(multi2di(bird.families)), bird.families)
### To see the results of randomly resolving a trichotomy:
tr <- read.tree(text = "(a:1,b:1,c:1);")
layout(matrix(1:4, 2, 2))
for (i in 1:4)
  plot(multi2di(tr), use.edge.length = FALSE, cex = 1.5)
layout(1)
}
\keyword{manip}
