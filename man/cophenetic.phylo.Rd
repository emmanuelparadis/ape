\name{cophenetic.phylo}
\alias{cophenetic.phylo}
\alias{dist.nodes}
\title{Pairwise Distances from a Phylogenetic Tree}
\usage{
\method{cophenetic}{phylo}(x)
dist.nodes(x)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}.}
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
