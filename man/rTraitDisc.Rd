\name{rTraitDisc}
\alias{rTraitDisc}
\title{Discrete Character Simulation}
\usage{
rTraitDisc(phy, model = "ER", k = if (is.matrix(model)) ncol(model) else 2,
           rate = 0.1, states = LETTERS[1:k], freq = rep(1/k, k),
           ancestor = FALSE, root.value = 1, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{model}{a character, a square numeric matrix, or a function
    specifying the model (see details).}
  \item{k}{the number of states of the character.}
  \item{rate}{the rate of change used if \code{model} is a character; it
    is \emph{not} recycled if \code{model = "ARD"} of \code{model =
      "SYM"}.}
  \item{states}{the labels used for the states; by default ``A'', ``B'',
    \dots}
  \item{freq}{a numeric vector giving the equilibrium relative
    frequencies of each state; by default the frequencies are equal.}
  \item{ancestor}{a logical value specifying whether to return the
    values at the nodes as well (by default, only the values at the tips
    are returned).}
  \item{root.value}{an integer giving the value at the root (by default,
    it's the first state). To have a random value, use \code{root.value
      = sample(k)}.}
  \item{\dots}{further arguments passed to \code{model} if it is a
    function.}
}
\description{
  This function simulates the evolution of a discrete character along a
  phylogeny. If \code{model} is a character or a matrix, evolution is
  simulated with a Markovian model; the transition probabilities are
  calculated for each branch with \eqn{P = e^{Qt}} where \eqn{Q} is the
  rate matrix given by \code{model} and \eqn{t} is the branch length.
  The calculation is done recursively from the root. See Paradis (2006,
  p. 101) for a general introduction applied to evolution.
}
\details{
  There are three possibilities to specify \code{model}:

\itemize{
  \item{A matrix:}{it must be a numeric square matrix; the diagonal is
    always ignored. The arguments \code{k} and \code{rate} are ignored.}

  \item{A character:}{these are the same short-cuts than in the function
  \code{\link{ace}}: \code{"ER"} is an equal-rates model, \code{"ARD"}
  is an all-rates-different model, and \code{"SYM"} is a symmetrical
  model. Note that the argument \code{rate} must be of the appropriate
  length, i.e., 1, \eqn{k(k - 1)}, or \eqn{k(k - 1)/2} for the three models,
  respectively. The rate matrix \eqn{Q} is then filled column-wise.}

  \item{A function:}{it must be of the form \code{foo(x, l)} where
    \code{x} is the trait of the ancestor and \code{l} is the branch
    length. It must return the value of the descendant as an integer.}
}}
\value{
  A factor with names taken from the tip labels of \code{phy}. If
  \code{ancestor = TRUE}, the node labels are used if present,
  otherwise, ``Node1'', ``Node2'', etc.
}
\references{
  Paradis, E. (2006) \emph{Analyses of Phylogenetics and Evolution with
    R.} New York: Springer.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{rTraitCont}}, \code{\link{rTraitMult}}, \code{\link{ace}}
}
\examples{
data(bird.orders)
### the two followings are the same:
rTraitDisc(bird.orders)
rTraitDisc(bird.orders, model = matrix(c(0, 0.1, 0.1, 0), 2))

### two-state model with irreversibility:
rTraitDisc(bird.orders, model = matrix(c(0, 0, 0.1, 0), 2))

### simple two-state model:
tr <- rcoal(n <- 40, br = runif)
x <- rTraitDisc(tr, ancestor = TRUE)
plot(tr, show.tip.label = FALSE)
nodelabels(pch = 19, col = x[-(1:n)])
tiplabels(pch = 19, col = x[1:n])

### an imaginary model with stasis 0.5 time unit after a node, then
### random evolution:
foo <- function(x, l) {
    if (l < 0.5) return(x)
    sample(2, size = 1)
}
tr <- rcoal(20, br = runif)
x <- rTraitDisc(tr, foo, ancestor = TRUE)
plot(tr, show.tip.label = FALSE)
co <- c("blue", "yellow")
cot <- c("white", "black")
Y <- x[1:20]
A <- x[-(1:20)]
nodelabels(A, bg = co[A], col = cot[A])
tiplabels(Y, bg = co[Y], col = cot[Y])
}
\keyword{datagen}
