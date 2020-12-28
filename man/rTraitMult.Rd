\name{rTraitMult}
\alias{rTraitMult}
\title{Multivariate Character Simulation}
\description{
  This function simulates the evolution of a multivariate set of traits
  along a phylogeny. The calculation is done recursively from the
  root.
}
\usage{
rTraitMult(phy, model, p = 1, root.value = rep(0, p), ancestor = FALSE,
           asFactor = NULL, trait.labels = paste("x", 1:p, sep = ""), ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{model}{a function specifying the model (see details).}
  \item{p}{an integer giving the number of traits.}
  \item{root.value}{a numeric vector giving the values at the root.}
  \item{ancestor}{a logical value specifying whether to return the
    values at the nodes as well (by default, only the values at the tips
    are returned).}
  \item{asFactor}{the indices of the traits that are returned as factors
    (discrete traits).}
  \item{trait.labels}{a vector of mode character giving the names of the
    traits.}
  \item{\dots}{further arguments passed to \code{model} if it is a
    function.}
}
\details{
  The model is specified with an \R function of the form \code{foo(x,
    l)} where \code{x} is a vector of the traits of the ancestor and
  \code{l} is the branch length. Other arguments may be added. The
  function must return a vector of length \code{p}.
}
\value{
  A data frame with \code{p} columns whose names are given by
  \code{trait.labels} and row names taken from the labels of the tree.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{rTraitCont}}, \code{\link{rTraitDisc}}, \code{\link{ace}}
}
\examples{
## correlated evolution of 2 continuous traits:
mod <- function(x, l) {
    y1 <- rnorm(1, x[1] + 0.5*x[2], 0.1)
    y2 <- rnorm(1, 0.5*x[1] + x[2], 0.1)
    c(y1, y2)
}
set.seed(11)
tr <- makeNodeLabel(rcoal(20))
x <- rTraitMult(tr, mod, 2, ancestor = TRUE)
op <- par(mfcol = c(2, 1))
plot(x, type = "n")
text(x, labels = rownames(x), cex = 0.7)
oq <- par(mar = c(0, 1, 0, 1), xpd = TRUE)
plot(tr, font = 1, cex = 0.7)
nodelabels(tr$node.label, cex = 0.7, adj = 1)
par(c(op, oq))
}
\keyword{datagen}
