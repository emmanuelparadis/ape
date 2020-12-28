\name{is.monophyletic}
\alias{is.monophyletic}
\title{
  Is Group Monophyletic
}
\usage{
is.monophyletic(phy, tips, reroot = !is.rooted(phy), plot = FALSE, ...)
}
\description{
    This function tests whether a list of tip labels is monophyletic on a given tree.
}
\arguments{
    \item{phy}{
        a phylogenetic tree description of class \code{"phylo"}.
    }
    \item{tips}{
       a vector of mode numeric or character specifying the tips to be tested.
    }
    \item{reroot}{
       a logical. If \code{FALSE}, then the input tree is not unrooted before the test.
    }
    \item{plot}{
        a logical. If \code{TRUE}, then the tree is plotted with the specified group \code{tips} highlighted.
    }
    \item{\dots}{
       further arguments passed to \code{plot}.
    }
}
\details{
    If \code{phy} is rooted, the test is done on the rooted tree, otherwise
    the tree is first unrooted, then arbitrarily rerooted, in order to be
    independent on the current position of the root. That is, the test
    asks if \code{tips} could be monophyletic given any favourably rooting
    of \code{phy}.

    If \code{phy} is unrooted the test is done on an unrooted tree, unless
    \code{reroot = FALSE} is specified.

    If tip labels in the list \code{tips} are given as characters, they need
    to be spelled as in the object \code{phy}.
}
\value{
    \code{TRUE} or \code{FALSE}.
}
\author{
    Johan Nylander \email{jnylander@users.sourceforge.net}
}
\seealso{
    \code{\link{which.edge}}, \code{\link{drop.tip}}, \code{\link{mrca}}.
}
\examples{
    ## Test one monophyletic and one paraphyletic group on the bird.orders tree
    \dontrun{data("bird.orders")}
    \dontrun{is.monophyletic(phy = bird.orders, tips = c("Ciconiiformes", "Gruiformes"))}
    \dontrun{is.monophyletic(bird.orders, c("Passeriformes", "Ciconiiformes", "Gruiformes"))}
    \dontshow{\dontrun{rm(bird.orders)}}
}
\keyword{utilities}

