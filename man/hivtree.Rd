\name{hivtree}
\alias{hivtree}
\alias{hivtree.newick}
\alias{hivtree.table}
\title{Phylogenetic Tree of 193 HIV-1 Sequences}
\description{
  This data set describes an estimated clock-like phylogeny of 193 HIV-1
  group M sequences sampled in the Democratic Republic of Congo.
}
\usage{
data(hivtree.newick)
data(hivtree.table)
}
\format{
  \code{hivtree.newick} is a string with the tree in Newick format.
  The data frame \code{hivtree.table} contains the corresponding internode
  distances.
}
\source{This is a data example from Strimmer and Pybus (2001).}
\references{
  Strimmer, K. and Pybus, O. G. (2001) Exploring the demographic history
  of DNA sequences using the generalized skyline plot. \emph{Molecular
    Biology and Evolution}, \bold{18}, 2298--2305.
}
\seealso{
  \code{\link{coalescent.intervals}}, \code{\link{collapsed.intervals}}
}
\keyword{datasets}
