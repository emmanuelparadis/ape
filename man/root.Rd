\name{root}
\alias{root}
\alias{root.phylo}
\alias{root.multiPhylo}
\alias{unroot}
\alias{unroot.phylo}
\alias{unroot.multiPhylo}
\alias{is.rooted}
\alias{is.rooted.phylo}
\alias{is.rooted.multiPhylo}
\title{Roots Phylogenetic Trees}
\description{
  \code{root} reroots a phylogenetic tree with respect to the specified
  outgroup or at the node specified in \code{node}.

  \code{unroot} unroots a phylogenetic tree, or returns it unchanged if
  it is already unrooted.

  \code{is.rooted} tests whether a tree is rooted.
}
\usage{
root(phy, ...)
\method{root}{phylo}(phy, outgroup, node = NULL, resolve.root = FALSE,
     interactive = FALSE, edgelabel = FALSE, ...)
\method{root}{multiPhylo}(phy, outgroup, ...)

unroot(phy)
\method{unroot}{phylo}(phy)
\method{unroot}{multiPhylo}(phy)

is.rooted(phy)
\method{is.rooted}{phylo}(phy)
\method{is.rooted}{multiPhylo}(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{outgroup}{a vector of mode numeric or character specifying the
    new outgroup.}
  \item{node}{alternatively, a node number where to root the tree.}
  \item{resolve.root}{a logical specifying whether to resolve the new
    root as a bifurcating node.}
  \item{interactive}{if \code{TRUE} the user is asked to select the node
    by clicking on the tree which must be plotted.}
  \item{edgelabel}{a logical value specifying whether to treat node
    labels as edge labels and thus eventually switching them so that
    they are associated with the correct edges when using
    \code{\link{drawSupportOnEdges}} (see Czech et al. 2016).}
  \item{\dots}{arguments passed among methods (e.g., when rooting lists
    of trees).}
}
\details{
  The argument \code{outgroup} can be either character or numeric. In
  the first case, it gives the labels of the tips of the new outgroup;
  in the second case the numbers of these labels in the vector
  \code{phy$tip.label} are given.

  If \code{outgroup} is of length one (i.e., a single value), then the
  tree is rerooted using the node below this tip as the new root.

  If \code{outgroup} is of length two or more, the most recent common
  ancestor (MRCA) \emph{of the ingroup} is used as the new root. Note
  that the tree is unrooted before being rerooted, so that if
  \code{outgroup} is already the outgroup, then the returned tree is not
  the same than the original one (see examples). If \code{outgroup} is
  not monophyletic, the operation fails and an error message is issued.

  If \code{resolve.root = TRUE}, \code{root} adds a zero-length branch
  below the MRCA of the ingroup.

  A tree is considered rooted if either only two branches connect to the
  root, or if there is a \code{root.edge} element. In all other cases,
  \code{is.rooted} returns \code{FALSE}.
}
\note{
  The use of \code{resolve.root = TRUE} together with \code{node = }
  gives an error if the specified node is the current root of the
  tree. This is because there is an ambiguity when resolving a node in
  an unrooted tree with no explicit outgroup. If the node is not the
  current root, the ambiguity is solved arbitrarily by considering the
  clade on the right of \code{node} (when the tree is plotted by
  default) as the ingroup. See a detailed explanation there:

  \url{https://www.mail-archive.com/r-sig-phylo@r-project.org/msg03805.html}.
}
\value{
  an object of class \code{"phylo"} or \code{"multiPhylo"} for
  \code{root} and \code{unroot}; a logical vector for \code{is.rooted}.
}
\references{
  Czech, L., Huerta-Cepas, J. and Stamatakis, A. (2017) A critical
  review on the use of support values in tree viewers and bioinformatics
  toolkits. \emph{Molecular Biology and Evolution}, \bold{34},
  1535--1542. \url{https://doi.org/10.1093/molbev/msx055}
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{bind.tree}}, \code{\link{drop.tip}},
  \code{\link{nodelabels}}, \code{\link{identify.phylo}}
}
\examples{
data(bird.orders)
plot(root(bird.orders, 1))
plot(root(bird.orders, 1:5))

tr <- root(bird.orders, 1)
is.rooted(bird.orders) # yes
is.rooted(tr)          # no
### This is because the tree has been unrooted first before rerooting.
### You can delete the outgroup...
is.rooted(drop.tip(tr, "Struthioniformes"))
### ... or resolve the basal trichotomy in two ways:
is.rooted(multi2di(tr))
is.rooted(root(bird.orders, 1, r = TRUE))
### To keep the basal trichotomy but forcing the tree as rooted:
tr$root.edge <- 0
is.rooted(tr)

x <- setNames(rmtree(10, 10), LETTERS[1:10])
is.rooted(x)
}
\keyword{manip}
