\name{rTraitCont}
\alias{rTraitCont}
\title{Continuous Character Simulation}
\usage{
rTraitCont(phy, model = "BM", sigma = 0.1, alpha = 1, theta = 0,
           ancestor = FALSE, root.value = 0, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{model}{a character (either \code{"BM"} or \code{"OU"}) or a
    function specifying the model (see details).}
  \item{sigma}{a numeric vector giving the standard-deviation of the
    random component for each branch (can be a single value).}
  \item{alpha}{if \code{model = "OU"}, a numeric vector giving the
    strength of the selective constraint for each branch (can be a
    single value).}
  \item{theta}{if \code{model = "OU"}, a numeric vector giving the
    optimum for each branch (can be a single value).}
  \item{ancestor}{a logical value specifying whether to return the
    values at the nodes as well (by default, only the values at the tips
    are returned).}
  \item{root.value}{a numeric giving the value at the root.}
  \item{\dots}{further arguments passed to \code{model} if it is a
    function.}
}
\description{
  This function simulates the evolution of a continuous character along a
  phylogeny. The calculation is done recursively from the root. See
  Paradis (2012, pp. 232 and 324) for an introduction.
}
\details{
  There are three possibilities to specify \code{model}:

\itemize{
  \item{\code{"BM"}:}{a Browian motion model is used. If the arguments
  \code{sigma} has more than one value, its length must be equal to the
  the branches of the tree. This allows to specify a model with variable
  rates of evolution. You must be careful that branch numbering is done
  with the tree in ``postorder'' order: to see the order of the branches
  you can use: \code{tr <- reorder(tr, "po"); plor(tr); edgelabels()}.
  The arguments \code{alpha} and \code{theta} are ignored.}

  \item{\code{"OU"}:}{an Ornstein-Uhlenbeck model is used. The above
  indexing rule is used for the three parameters \code{sigma},
  \code{alpha}, and \code{theta}. This may be interesting for the last
  one to model varying phenotypic optima. The exact updating formula
  from Gillespie (1996) are used which are reduced to BM formula if
  \code{alpha = 0}.}

  \item{A function:}{it must be of the form \code{foo(x, l)} where
  \code{x} is the trait of the ancestor and \code{l} is the branch
  length. It must return the value of the descendant. The arguments
  \code{sigma}, \code{alpha}, and \code{theta} are ignored.}
}}
\value{
  A numeric vector with names taken from the tip labels of
  \code{phy}. If \code{ancestor = TRUE}, the node labels are used if
  present, otherwise, ``Node1'', ``Node2'', etc.
}
\references{
  Gillespie, D. T. (1996) Exact numerical simulation of the
  Ornstein-Uhlenbeck process and its integral. \emph{Physical Review E},
  \bold{54}, 2084--2091.

  Paradis, E. (2012) \emph{Analysis of Phylogenetics and Evolution with
    R (Second Edition).} New York: Springer.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{rTraitDisc}}, \code{\link{rTraitMult}}, \code{\link{ace}}
}
\examples{
data(bird.orders)
rTraitCont(bird.orders) # BM with sigma = 0.1
### OU model with two optima:
tr <- reorder(bird.orders, "postorder")
plot(tr)
edgelabels()
theta <- rep(0, Nedge(tr))
theta[c(1:4, 15:16, 23:24)] <- 2
## sensitive to 'alpha' and 'sigma':
rTraitCont(tr, "OU", theta = theta, alpha=.1, sigma=.01)
### an imaginary model with stasis 0.5 time unit after a node, then
### BM evolution with sigma = 0.1:
foo <- function(x, l) {
    if (l <= 0.5) return(x)
    x + (l - 0.5)*rnorm(1, 0, 0.1)
}
tr <- rcoal(20, br = runif)
rTraitCont(tr, foo, ancestor = TRUE)
### a cumulative Poisson process:
bar <- function(x, l) x + rpois(1, l)
(x <- rTraitCont(tr, bar, ancestor = TRUE))
plot(tr, show.tip.label = FALSE)
Y <- x[1:20]
A <- x[-(1:20)]
nodelabels(A)
tiplabels(Y)
}
\keyword{datagen}
