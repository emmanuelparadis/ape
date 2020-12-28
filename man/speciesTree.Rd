\name{speciesTree}
\alias{speciesTree}
\title{Species Tree Estimation}
\description{
  This function calculates the species tree from a set of gene trees.
}
\usage{
speciesTree(x, FUN = min)
}
\arguments{
  \item{x}{a list of trees, e.g., an object of class
    \code{"multiPhylo"}.}
  \item{FUN}{a function used to compute the divergence times of each
    pair of tips.}
}
\details{
  For all trees in \code{x}, the divergence time of each pair of tips is
  calculated: these are then `summarized' with \code{FUN} to build a new
  distance matrix used to calculate the species tree with a
  single-linkage hierarchical clustering. The default for \code{FUN}
  computes the maximum tree (maxtree) of Liu et al. (2010). Using
  \code{FUN = mean} gives the shallowest divergence tree of Maddison and
  Knowles (2006).
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  Liu, L., Yu, L. and Pearl, D. K. (2010) Maximum tree: a consistent
  estimator of the species tree. \emph{Journal of Mathematical Biology},
  \bold{60}, 95--106.

  Maddison, W. P. and Knowles, L. L. (2006) Inferring phylogeny despite
  incomplete lineage sorting. \emph{Systematic Biology}, \bold{55}, 21--30.
}
\author{Emmanuel Paradis}
\examples{
### example in Liu et al. (2010):
tr1 <- read.tree(text = "(((B:0.05,C:0.05):0.01,D:0.06):0.04,A:0.1);")
tr2 <- read.tree(text = "(((A:0.07,C:0.07):0.02,D:0.09):0.03,B:0.12);")
TR <- c(tr1, tr2)
TSmax <- speciesTree(TR) # MAXTREE
TSsha <- speciesTree(TR, mean) # shallowest divergence

kronoviz(c(tr1, tr2, TSmax, TSsha), horiz = FALSE,
         type = "c", cex = 1.5, font = 1)
mtext(c("Gene tree 1", "Gene tree 2", "Species tree - MAXTREE"),
      at = -c(7.5, 4, 1))
mtext("Species tree - Shallowest Divergence")
layout(1)
}
\keyword{models}
