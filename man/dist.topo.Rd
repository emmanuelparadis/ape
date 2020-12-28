\name{dist.topo}
\alias{dist.topo}
\title{Topological Distances Between Two Trees}
\description{
  This function computes the topological distance between two
  phylogenetic trees or among trees in a list (if \code{y = NULL} using
  different methods.
}
\usage{
dist.topo(x, y = NULL, method = "PH85")
}
\arguments{
  \item{x}{an object of class \code{"phylo"} or of class
    \code{"multiPhylo"}.}
  \item{y}{an (optional) object of class \code{"phylo"}.}
  \item{method}{a character string giving the method to be used: either
    \code{"PH85"}, or \code{"score"}.}
}
\value{
  a single numeric value if both \code{x} and \code{y} are used, an
  object of class \code{"dist"} otherwise.
}
\details{
  Two methods are available: the one by Penny and Hendy (1985,
  originally from Robinson and Foulds 1981), and the branch length score
  by Kuhner and Felsenstein (1994). The trees are always considered as
  unrooted.

  The topological distance is defined as twice the number of internal
  branches defining different bipartitions of the tips (Robinson and
  Foulds 1981; Penny and Hendy 1985). Rzhetsky and Nei (1992) proposed a
  modification of the original formula to take multifurcations into
  account.

  The branch length score may be seen as similar to the previous
  distance but taking branch lengths into account. Kuhner and
  Felsenstein (1994) proposed to calculate the square root of the sum of
  the squared differences of the (internal) branch lengths defining
  similar bipartitions (or splits) in both trees.
}
\note{
  The geodesic distance of Billera et al. (2001) has been disabled: see
  the package \pkg{distory} on CRAN.
}
\references{
  Billera, L. J., Holmes, S. P. and Vogtmann, K. (2001) Geometry of the
  space of phylogenetic trees. \emph{Advances in Applied Mathematics},
  \bold{27}, 733--767.

  Kuhner, M. K. and Felsenstein, J. (1994) Simulation comparison of
  phylogeny algorithms under equal and unequal evolutionary rates.
  \emph{Molecular Biology and Evolution}, \bold{11}, 459--468.

  Nei, M. and Kumar, S. (2000) \emph{Molecular Evolution and
  Phylogenetics}. Oxford: Oxford University Press.

  Penny, D. and Hendy, M. D. (1985) The use of tree comparison
  metrics. \emph{Systemetic Zoology}, \bold{34}, 75--82.

  Robinson, D. F. and Foulds, L. R. (1981) Comparison of phylogenetic
  trees. \emph{Mathematical Biosciences}, \bold{53}, 131--147.

  Rzhetsky, A. and Nei, M. (1992) A simple method for estimating and
  testing minimum-evolution trees. \emph{Molecular Biology and
    Evolution}, \bold{9}, 945--967.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{cophenetic.phylo}}, \code{\link{prop.part}}
}
\examples{
ta <- rtree(30)
tb <- rtree(30)
dist.topo(ta, ta) # 0
dist.topo(ta, tb) # unlikely to be 0
}
\keyword{manip}
