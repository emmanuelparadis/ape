\name{summary.phylo}
\alias{summary.phylo}
\alias{Ntip}
\alias{Ntip.phylo}
\alias{Ntip.multiPhylo}
\alias{Nnode}
\alias{Nnode.phylo}
\alias{Nnode.multiPhylo}
\alias{Nedge}
\alias{Nedge.phylo}
\alias{Nedge.multiPhylo}
\title{Print Summary of a Phylogeny}
\usage{
\method{summary}{phylo}(object, \dots)

Ntip(phy)
\method{Ntip}{phylo}(phy)
\method{Ntip}{multiPhylo}(phy)

Nnode(phy, ...)
\method{Nnode}{phylo}(phy, internal.only = TRUE, ...)
\method{Nnode}{multiPhylo}(phy, internal.only = TRUE, ...)

Nedge(phy)
\method{Nedge}{phylo}(phy)
\method{Nedge}{multiPhylo}(phy)
}
\arguments{
  \item{object, phy}{an object of class \code{"phylo"} or
    \code{"multiPhylo"}.}
  \item{\dots}{further arguments passed to or from other methods.}
  \item{internal.only}{a logical indicating whether to return the number
    of internal nodes only (the default), or of internal and terminal
    (tips) nodes (if \code{FALSE}).}
}
\description{
  The first function prints a compact summary of a phylogenetic tree (an
  object of class \code{"phylo"}). The three other functions return the
  number of tips, nodes, or edges, respectively.
}
\details{
  The summary includes the numbers of tips and of nodes, summary
  statistics of the branch lengths (if they are available) with mean,
  variance, minimum, first quartile, median, third quartile, and
  maximum, listing of the first ten tip labels, and (if available) of
  the first ten node labels. It is also printed whether some of these
  optional elements (branch lengths, node labels, and root edge) are not
  found in the tree.

  \code{summary} simply prints its results on the standard output and is
  not meant for programming.
}
\value{
  A NULL value in the case of \code{summary}, a single numeric value for
  the three other functions.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.tree}}, \code{\link[base]{summary}} for the generic R
  function, \code{\link{multiphylo}}, \code{\link{c.phylo}}
}
\examples{
data(bird.families)
summary(bird.families)
Ntip(bird.families)
Nnode(bird.families)
Nedge(bird.families)
}
\keyword{manip}
