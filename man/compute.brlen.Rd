\name{compute.brlen}
\alias{compute.brlen}
\title{Branch Lengths Computation}
\usage{
compute.brlen(phy, method = "Grafen", power = 1, ...)
}
\arguments{
  \item{phy}{an object of class \code{phylo} representing the tree.}
  \item{method}{the method to be used to compute the branch lengths;
    this must be one of the followings: (i) \code{"Grafen"} (the
    default), (ii) a numeric vector, or (iii) a function.}
  \item{power}{The power at which heights must be raised (see below).}
  \item{\dots}{further argument(s) to be passed to \code{method} if it is
    a function.}
}
\description{
  This function computes branch lengths of a tree using different
  methods.
}
\details{
  Grafen's (1989) computation of branch lengths: each node is given a
  `height', namely the number of leaves of the subtree minus one, 0 for
  leaves. Each height is scaled so that root height is 1, and then
  raised at power 'rho' (> 0). Branch lengths are then computed as the
  difference between height of lower node and height of upper node.

  If one or several numeric values are provided as \code{method}, they
  are recycled if necessary. If a function is given instead, further
  arguments are given in place of \code{...} (they must be named, see
  examples).

  Zero-length branches are not treated as multichotomies, and thus may
  need to be collapsed (see \code{\link{di2multi}}).
}
\value{
  An object of class \code{phylo} with branch lengths.
}
\author{Julien Dutheil \email{dutheil@evolbio.mpg.de} and
  Emmanuel Paradis}
\references{
  Grafen, A. (1989) The phylogenetic regression. \emph{Philosophical
    Transactions of the Royal society of London. Series B. Biological
    Sciences}, \bold{326}, 119--157.
}
\seealso{
  \code{\link{read.tree}} for a description of \code{phylo} objects,
  \code{\link{di2multi}}, \code{\link{multi2di}}
}
\examples{
data(bird.orders)
plot(compute.brlen(bird.orders, 1))
plot(compute.brlen(bird.orders, runif, min = 0, max = 5))
layout(matrix(1:4, 2, 2))
plot(compute.brlen(bird.orders, power=1), main=expression(rho==1))
plot(compute.brlen(bird.orders, power=3), main=expression(rho==3))
plot(compute.brlen(bird.orders, power=0.5), main=expression(rho==0.5))
plot(compute.brlen(bird.orders, power=0.1), main=expression(rho==0.1))
layout(1)
}
\keyword{manip}
