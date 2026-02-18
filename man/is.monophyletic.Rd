\name{is.monophyletic}
\alias{is.monophyletic}
\title{Test If Group Monophyletic}
\description{
    This function tests whether a list of tip labels is monophyletic on
    a given tree.
}
\usage{
is.monophyletic(phy, tips, reroot = !is.rooted(phy), plot = FALSE, ...)
}
\arguments{
  \item{phy}{a phylogenetic tree description of class \code{"phylo"}.}
  \item{tips}{a vector of mode numeric or character specifying the tips
    to be tested.}
  \item{reroot}{a logical value. If \code{FALSE}, then the input tree is
    not unrooted before the test.}
  \item{plot}{a logical value. If \code{TRUE}, then the tree is plotted
    with the specified group \code{tips} highlighted.}
  \item{\dots}{further arguments passed to \code{plot}.}
}
\details{
    If \code{phy} is rooted, the test is done on the rooted tree, otherwise
    the tree is first unrooted, then arbitrarily rerooted, in order to be
    independent on the current position of the root. That is, the test
    asks if \code{tips} could be monophyletic given any favourably rooting
    of \code{phy}.

    If \code{phy} is unrooted the test is done on an unrooted tree, unless
    \code{reroot = FALSE} is specified.

    If labels in the list \code{tips} are given as characters, they need
    to be spelled as in the object \code{phy}.
}
\value{a logical value.}
\author{
    Johan Nylander \email{jnylander@users.sourceforge.net}
}
\seealso{
    \code{\link{which.edge}}, \code{\link{drop.tip}}, \code{\link{mrca}}.
}
\examples{
data(bird.orders)
## Test one monophyletic and one paraphyletic group:
is.monophyletic(phy = bird.orders, tips = c("Ciconiiformes", "Gruiformes"))
is.monophyletic(bird.orders, c("Passeriformes", "Ciconiiformes", "Gruiformes"))
}
\keyword{utilities}

