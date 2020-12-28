\name{boot.phylo}
\alias{boot.phylo}
\alias{prop.part}
\alias{prop.clades}
\alias{print.prop.part}
\alias{summary.prop.part}
\alias{plot.prop.part}
\title{Tree Bipartition and Bootstrapping Phylogenies}
\description{
  These functions analyse bipartitions found in a series of trees.

  \code{prop.part} counts the number of bipartitions found in a series
  of trees given as \code{\dots}. If a single tree is passed, the
  returned object is a list of vectors with the tips descending from
  each node (i.e., clade compositions indexed by node number).

  \code{prop.clades} counts the number of times the bipartitions present
  in \code{phy} are present in a series of trees given as \code{\dots} or
  in the list previously computed and given with \code{part}.

  \code{boot.phylo} performs a bootstrap analysis.
}
\usage{
boot.phylo(phy, x, FUN, B = 100, block = 1,
           trees = FALSE, quiet = FALSE,
           rooted = is.rooted(phy), jumble = TRUE,
            mc.cores = 1)
prop.part(..., check.labels = TRUE)
prop.clades(phy, ..., part = NULL, rooted = FALSE)
\method{print}{prop.part}(x, ...)
\method{summary}{prop.part}(object, ...)
\method{plot}{prop.part}(x, barcol = "blue", leftmar = 4, col = "red", ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{x}{in the case of \code{boot.phylo}: a taxa (rows) by characters
    (columns) matrix; in the case of \code{print} and \code{plot}: an
    object of class \code{"prop.part"}.}
  \item{FUN}{the function used to estimate \code{phy} (see details).}
  \item{B}{the number of bootstrap replicates.}
  \item{block}{the number of columns in \code{x} that will be resampled
    together (see details).}
  \item{trees}{a logical specifying whether to return the bootstraped
    trees (\code{FALSE} by default).}
  \item{quiet}{a logical: a progress bar is displayed by default.}
  \item{rooted}{a logical specifying whether the trees should be treated
    as rooted or not.}
  \item{jumble}{a logical value. By default, the rows of \code{x} are
    randomized to avoid artificially too large bootstrap values
    associated with very short branches.}
  \item{mc.cores}{the number of cores (CPUs) to be used (passed to
    \pkg{parallel}).}
  \item{\dots}{either (i) a single object of class \code{"phylo"}, (ii) a
    series of such objects separated by commas, or (iii) a list
    containing such objects. In the case of \code{plot} further
    arguments for the plot (see details).}
  \item{check.labels}{a logical specifying whether to check the labels
    of each tree. If \code{FALSE}, it is assumed that all trees have the
    same tip labels, and that they are in the same order (see details).}
  \item{part}{a list of partitions as returned by \code{prop.part}; if
    this is used then \code{\dots} is ignored.}
  \item{object}{an object of class \code{"prop.part"}.}
  \item{barcol}{the colour used for the bars displaying the number of
    partitions in the upper panel.}
  \item{leftmar}{the size of the margin on the left to display the tip
    labels.}
  \item{col}{the colour used to visualise the bipartitions.}
}
\details{
  The argument \code{FUN} in \code{boot.phylo} must be the function used
  to estimate the tree from the original data matrix. Thus, if the tree
  was estimated with neighbor-joining (see \code{nj}), one maybe wants
  something like \code{FUN = function(xx) nj(dist.dna(xx))}.

  \code{block} in \code{boot.phylo} specifies the number of columns to
  be resampled altogether. For instance, if one wants to resample at the
  codon-level, then \code{block = 3} must be used.

  Using \code{check.labels = FALSE} in \code{prop.part} decreases
  computing times. This requires that (i) all trees have the same tip
  labels, \emph{and} (ii) these labels are ordered similarly in all
  trees (in other words, the element \code{tip.label} are identical in
  all trees).

  The plot function represents a contingency table of the different
  partitions (on the \emph{x}-axis) in the lower panel, and their observed
  numbers in the upper panel. Any further arguments (\dots) are used to
  change the aspects of the points in the lower panel: these may be
  \code{pch}, \code{col}, \code{bg}, \code{cex}, etc. This function
  works only if there is an attribute \code{labels} in the object.

  The print method displays the partitions and their numbers. The
  summary method extracts the numbers only.
}
\note{
  \code{prop.clades} calls internally \code{prop.part} with the option
  \code{check.labels = TRUE}, which may be very slow. If the trees
  passed as \code{\dots} fulfills conditions (i) and (ii) above, then it
  might be faster to first call, e.g., \code{pp <- prop.part(...)}, then
  use the option \code{part}: \code{prop.clades(phy, part = pp)}.

  Since \pkg{ape} 3.5, \code{prop.clades} should return sensible results
  for all values of \code{rooted}: if \code{FALSE}, the numbers of
  bipartitions (or splits); if \code{TRUE}, the number of clades (of
  hopefully rooted trees).
}
\value{
  \code{prop.part} returns an object of class \code{"prop.part"} which
  is a list with an attribute \code{"number"}. The elements of this list
  are the observed clades, and the attribute their respective
  numbers. If the default \code{check.labels = FALSE} is used, an
  attribute \code{"labels"} is added, and the vectors of the returned
  object contains the indices of these labels instead of the labels
  themselves.

  \code{prop.clades} and \code{boot.phylo} return a numeric vector
  which \emph{i}th element is the number associated to the \emph{i}th
  node of \code{phy}. If \code{trees = TRUE}, \code{boot.phylo} returns
  a list whose first element (named \code{"BP"}) is like before, and the
  second element (\code{"trees"}) is a list with the bootstraped
  trees.

  \code{summary} returns a numeric vector.
}
\references{
  Efron, B., Halloran, E. and Holmes, S. (1996) Bootstrap confidence
  levels for phylogenetic trees. \emph{Proceedings of the National
    Academy of Sciences USA}, \bold{93}, 13429--13434.

  Felsenstein, J. (1985) Confidence limits on phylogenies: an approach
  using the bootstrap. \emph{Evolution}, \bold{39}, 783--791.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{as.bitsplits}}, \code{\link{dist.topo}},
  \code{\link{consensus}}, \code{\link{nodelabels}}
}
\examples{
data(woodmouse)
f <- function(x) nj(dist.dna(x))
tr <- f(woodmouse)
### Are bootstrap values stable?
for (i in 1:5)
  print(boot.phylo(tr, woodmouse, f, quiet = TRUE))
### How many partitions in 100 random trees of 10 labels?...
TR <- rmtree(100, 10)
pp10 <- prop.part(TR)
length(pp10)
### ... and in 100 random trees of 20 labels?
TR <- rmtree(100, 20)
pp20 <- prop.part(TR)
length(pp20)
plot(pp10, pch = "x", col = 2)
plot(pp20, pch = "x", col = 2)

set.seed(2)
tr <- rtree(10) # rooted
## the following used to return a wrong result with ape <= 3.4:
prop.clades(tr, tr)
prop.clades(tr, tr, rooted = TRUE)
tr <- rtree(10, rooted = FALSE)
prop.clades(tr, tr) # correct

### an illustration of the use of prop.clades with bootstrap trees:

fun <- function(x) as.phylo(hclust(dist.dna(x), "average")) # upgma() in phangorn
tree <- fun(woodmouse)
## get 100 bootstrap trees:
bstrees <- boot.phylo(tree, woodmouse, fun, trees = TRUE)$trees
## get proportions of each clade:
clad <- prop.clades(tree, bstrees, rooted = TRUE)
## get proportions of each bipartition:
boot <- prop.clades(tree, bstrees)
layout(1)
par(mar = rep(2, 4))
plot(tree, main = "Bipartition vs. Clade Support Values")
drawSupportOnEdges(boot)
nodelabels(clad)
legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
       pt.bg = c("green", "lightblue"), pt.cex = 2.5)

\dontrun{
## an example of double bootstrap:
nrep1 <- 100
nrep2 <- 100
p <- ncol(woodmouse)
DB <- 0

for (b in 1:nrep1) {
    X <- woodmouse[, sample(p, p, TRUE)]
    DB <- DB + boot.phylo(tr, X, f, nrep2, quiet = TRUE)
}
DB
## to compare with:
boot.phylo(tr, woodmouse, f, 1e4)
}
}
\keyword{manip}
\keyword{htest}
