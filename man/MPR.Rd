\name{MPR}
\alias{MPR}
\title{Most Parsimonious Reconstruction}
\description{
  This function does ancestral character reconstruction by parsimony as
  described in Hanazawa et al. (1995) and modified by Narushima and
  Hanazawa (1997).
}
\usage{
MPR(x, phy, outgroup)
}
\arguments{
  \item{x}{a vector of integers.}
  \item{phy}{an object of class \code{"phylo"}; the tree must be
    unrooted and fully dichotomous.}
  \item{outgroup}{an integer or a character string giving the tip of
    \code{phy} used as outgroup.}
}
\details{
  Hanazawa et al. (1995) and Narushima and Hanazawa (1997) used Farris's
  (1970) and Swofford and Maddison's (1987) framework to reconstruct
  ancestral states using parsimony. The character is assumed to take
  integer values. The algorithm finds the sets of values for each node
  as intervals with lower and upper values.

  It is recommended to root the tree with the outgroup before the
  analysis, so plotting the values with \code{\link{nodelabels}} is
  simple.
}
\value{
  a matrix of integers with two columns named ``lower'' and ``upper''
  giving the lower and upper values of the reconstructed sets for each
  node.
}
\references{
  Farris, J. M. (1970) Methods for computing Wagner trees.
  \emph{Systematic Zoology}, \bold{19}, 83--92.

  Hanazawa, M., Narushima, H. and Minaka, N. (1995) Generating most
  parsimonious reconstructions on a tree: a generalization of the
  Farris--Swofford--Maddison method. \emph{Discrete Applied
    Mathematics}, \bold{56}, 245--265.

  Narushima, H. and Hanazawa, M. (1997) A more efficient algorithm for
  MPR problems in phylogeny. \emph{Discrete Applied Mathematics},
  \bold{80}, 231--238.

  Swofford, D. L. and Maddison, W. P. (1987) Reconstructing ancestral
  character states under Wagner parsimony. \emph{Mathematical
    Biosciences}, \bold{87}, 199--229.
}\author{Emmanuel Paradis}
\seealso{
  \code{\link{ace}}, \code{\link{root}}, \code{\link{nodelabels}}
}
\examples{
## the example in Narushima and Hanazawa (1997):
tr <- read.tree(text = "(((i,j)c,(k,l)b)a,(h,g)e,f)d;")
x <- c(1, 3, 0, 6, 5, 2, 4)
names(x) <- letters[6:12]
(o <- MPR(x, tr, "f"))
plot(tr)
nodelabels(paste0("[", o[, 1], ",", o[, 2], "]"))
tiplabels(x[tr$tip.label], adj = -2)
## some random data:
x <- rpois(30, 1)
tr <- rtree(30, rooted = FALSE)
MPR(x, tr, "t1")
}
\keyword{models}
