\name{pic}
\alias{pic}
\title{Phylogenetically Independent Contrasts}
\description{
  Compute the phylogenetically independent contrasts using the method
  described by Felsenstein (1985).
}
\usage{
pic(x, phy, scaled = TRUE, var.contrasts = FALSE,
    rescaled.tree = FALSE)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{scaled}{logical, indicates whether the contrasts should be
    scaled with their expected variances (default to \code{TRUE}).}
  \item{var.contrasts}{logical, indicates whether the expected
    variances of the contrasts should be returned (default to
    \code{FALSE}).}
  \item{rescaled.tree}{logical, if \code{TRUE} the rescaled tree is
    returned together with the main results.}
}
\details{
  If \code{x} has names, its values are matched to the tip labels of
  \code{phy}, otherwise its values are taken to be in the same order
  than the tip labels of \code{phy}.

  The user must be careful here since the function requires that both
  series of names perfectly match. If both series of names do not match,
  the values in the \code{x} are taken to be in the same order than the
  tip labels of \code{phy}, and a warning message is issued.
}
\value{
  either a vector of phylogenetically independent contrasts (if
  \code{var.contrasts = FALSE}), or a two-column matrix with the
  phylogenetically independent contrasts in the first column and their
  expected variance in the second column (if \code{var.contrasts =
  TRUE}). If the tree has node labels, these are used as labels of the
  returned object.

  If \code{rescaled.tree = TRUE}, a list is returned with two elements
  named ``contr'' with the above results and ``rescaled.tree'' with the
  tree and its rescaled branch lengths (see Felsenstein 1985).
}
\references{
  Felsenstein, J. (1985) Phylogenies and the comparative method.
  \emph{American Naturalist}, \bold{125}, 1--15.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.tree}}, \code{\link{compar.gee}},
  \code{\link{compar.lynch}}, \code{\link{pic.ortho}},
  \code{\link{varCompPhylip}}
}
\examples{
### The example in Phylip 3.5c (originally from Lynch 1991)
cat("((((Homo:0.21,Pongo:0.21):0.28,",
   "Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);",
   file = "ex.tre", sep = "\n")
tree.primates <- read.tree("ex.tre")
X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
names(X) <- names(Y) <- c("Homo", "Pongo", "Macaca", "Ateles", "Galago")
pic.X <- pic(X, tree.primates)
pic.Y <- pic(Y, tree.primates)
cor.test(pic.X, pic.Y)
lm(pic.Y ~ pic.X - 1) # both regressions
lm(pic.X ~ pic.Y - 1) # through the origin
unlink("ex.tre") # delete the file "ex.tre"
}
\keyword{regression}
