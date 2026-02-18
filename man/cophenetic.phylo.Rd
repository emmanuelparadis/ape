\name{cophenetic.phylo}
\alias{cophenetic.phylo}
\alias{dist.nodes}
\title{Pairwise Distances from a Phylogenetic Tree}
\usage{
\method{cophenetic}{phylo}(x)
dist.nodes(x, fail.if.no.length = FALSE)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
  \item{fail.if.no.length}{a logical values. If the tree has no branch
    lengths, these are all set equal to one (with a warning). If you
    prefer to catch the case of no branch lengths with an error, set
    this option to \code{TRUE}.}
}
\description{
  \code{cophenetic.phylo} computes the pairwise distances between the
  pairs of tips from a phylogenetic tree using its branch lengths.

  \code{dist.nodes} does the same but between all nodes, internal and
  terminal, of the tree.
}
\value{
  a numeric matrix with colnames and rownames set to the names of the
  tips (as given by the element \code{tip.label} of the argument
  \code{phy}), or, in the case of \code{dist.nodes}, the numbers of the
  tips and the nodes (as given by the element \code{edge}).
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.tree}} to read tree files in Newick format,
  \code{\link[stats]{cophenetic}} for the generic function
}
\keyword{manip}
