\name{balance}
\alias{balance}
\title{Balance of a Dichotomous Phylogenetic Tree}
\usage{
balance(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
}
\description{
  This function computes the balance of a phylogenetic tree, that is for
  each node of the tree the numbers of descendants (i.e. tips) on each
  of its daughter-branch. The tree must be fully dichotomous.
}
\value{
  a numeric matrix with two columns and one row for each node of the
  tree. The columns give the numbers of descendants on each
  daughter-branches (the order of both columns being arbitrary). If the
  phylogeny \code{phy} has an element \code{node.label}, this is used as
  rownames for the returned matrix; otherwise the numbers (of mode
  character) of the matrix \code{edge} of \code{phy} are used as rownames.
}
\references{
  Aldous, D. J. (2001) Stochastic models and descriptive statistics for
  phylogenetic trees, from Yule to today. \emph{Statistical Science},
  \bold{16}, 23--34.
}
\author{Emmanuel Paradis}
\keyword{manip}
