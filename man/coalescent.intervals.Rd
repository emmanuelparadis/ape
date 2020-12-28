\name{coalescent.intervals}
\alias{coalescent.intervals}
\alias{coalescent.intervals.phylo}
\alias{coalescent.intervals.default}
\title{Coalescent Intervals}
\usage{
coalescent.intervals(x)
}
\arguments{
  \item{x}{either an ultra-metric phylogenetic tree (i.e. an object of
    class \code{"phylo"}) or, alternatively, a vector of interval lengths.}
}
\description{
 This function extracts or generates information about coalescent intervals
 (number of lineages, interval lengths, interval count, total depth) from
 a phylogenetic tree or a list of internode distances. The input tree
 needs to be ultra-metric (i.e. clock-like).
}
\value{
An object of class \code{"coalescentIntervals"} with the following entries:

  \item{lineages}{ A vector with the number of lineages at the start of each coalescent
    interval.}
  \item{interval.length}{ A vector with the length of each coalescent
    interval.}
  \item{interval.count}{ The total number of coalescent
    intervals.}
  \item{total.depth}{ The sum of the lengths of all coalescent
    intervals.}
}
\seealso{
\code{\link{branching.times}}, \code{\link{collapsed.intervals}},
\code{\link{read.tree}}.
}
\author{Korbinian Strimmer}
\examples{
data("hivtree.newick") # example tree in NH format
tree.hiv <- read.tree(text = hivtree.newick) # load tree
ci <- coalescent.intervals(tree.hiv) # from tree
ci
data("hivtree.table") # same tree, but in table format
ci <- coalescent.intervals(hivtree.table$size) # from vector of interval lengths
ci
}
\keyword{manip}
