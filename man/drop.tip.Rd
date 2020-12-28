\name{drop.tip}
\alias{drop.tip}
\alias{keep.tip}
\alias{extract.clade}
\title{Remove Tips in a Phylogenetic Tree}
\description{
  \code{drop.tip} removes the terminal branches of a phylogenetic tree,
  possibly removing the corresponding internal branches. \code{keep.tip}
  does the opposite operation (i.e., returns the induced tree).

  \code{extract.clade} does the inverse operation: it keeps all the tips
  from a given node, and deletes all the other tips.
}
\usage{
drop.tip(phy, tip, trim.internal = TRUE, subtree = FALSE,
         root.edge = 0, rooted = is.rooted(phy), collapse.singles = TRUE,
         interactive = FALSE)
keep.tip(phy, tip)
extract.clade(phy, node, root.edge = 0, collapse.singles = TRUE,
              interactive = FALSE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{tip}{a vector of mode numeric or character specifying the tips
    to delete.}
  \item{trim.internal}{a logical specifying whether to delete the
    corresponding internal branches.}
  \item{subtree}{a logical specifying whether to output in the tree how
    many tips have been deleted and where.}
  \item{root.edge}{an integer giving the number of internal branches to
    be used to build the new root edge. This has no effect if
    \code{trim.internal = FALSE}.}
  \item{rooted}{a logical indicating whether the tree must be treated as
    rooted or not. This allows to force the tree to be considered as
    unrooted (see examples).}
  \item{collapse.singles}{a logical specifying whether to delete the
    internal nodes of degree 2.}
  \item{node}{a node number or label.}
  \item{interactive}{if \code{TRUE} the user is asked to select the tips
    or the node by clicking on the tree which must be plotted.}
}
\details{
  The argument \code{tip} can be either character or numeric. In the
  first case, it gives the labels of the tips to be deleted; in the
  second case the numbers of these labels in the vector
  \code{phy$tip.label} are given.

  This also applies to \code{node}, but if this argument is character
  and the tree has no node label, this results in an error. If more than
  one value is given with \code{node} (i.e., a vector of length two or
  more), only the first one is used with a warning.

  If \code{trim.internal = FALSE}, the new tips are given \code{"NA"} as
  labels, unless there are node labels in the tree in which case they
  are used.

  If \code{subtree = TRUE}, the returned tree has one or several
  terminal branches indicating how many tips have been removed (with a
  label \code{"[x_tips]"}). This is done for as many monophyletic groups
  that have been deleted.

  Note that \code{subtree = TRUE} implies \code{trim.internal = TRUE}.

  To undestand how the option \code{root.edge} works, see the examples
  below.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis, Klaus Schliep, Joseph Brown}
\seealso{
  \code{\link{bind.tree}}, \code{\link{root}}
}
\examples{
data(bird.families)
tip <- c(
"Eopsaltriidae", "Acanthisittidae", "Pittidae", "Eurylaimidae",
"Philepittidae", "Tyrannidae", "Thamnophilidae", "Furnariidae",
"Formicariidae", "Conopophagidae", "Rhinocryptidae", "Climacteridae",
"Menuridae", "Ptilonorhynchidae", "Maluridae", "Meliphagidae",
"Pardalotidae", "Petroicidae", "Irenidae", "Orthonychidae",
"Pomatostomidae", "Laniidae", "Vireonidae", "Corvidae",
"Callaeatidae", "Picathartidae", "Bombycillidae", "Cinclidae",
"Muscicapidae", "Sturnidae", "Sittidae", "Certhiidae",
"Paridae", "Aegithalidae", "Hirundinidae", "Regulidae",
"Pycnonotidae", "Hypocoliidae", "Cisticolidae", "Zosteropidae",
"Sylviidae", "Alaudidae", "Nectariniidae", "Melanocharitidae",
"Paramythiidae","Passeridae", "Fringillidae")
plot(drop.tip(bird.families, tip))
plot(drop.tip(bird.families, tip, trim.internal = FALSE))
data(bird.orders)
plot(drop.tip(bird.orders, 6:23, subtree = TRUE))
plot(drop.tip(bird.orders, c(1:5, 20:23), subtree = TRUE))
plot(drop.tip(bird.orders, c(1:20, 23), subtree = TRUE))
plot(drop.tip(bird.orders, c(1:20, 23), subtree = TRUE, rooted = FALSE))
### Examples of the use of `root.edge'
tr <- read.tree(text = "(A:1,(B:1,(C:1,(D:1,E:1):1):1):1):1;")
drop.tip(tr, c("A", "B"), root.edge = 0) # = (C:1,(D:1,E:1):1);
drop.tip(tr, c("A", "B"), root.edge = 1) # = (C:1,(D:1,E:1):1):1;
drop.tip(tr, c("A", "B"), root.edge = 2) # = (C:1,(D:1,E:1):1):2;
drop.tip(tr, c("A", "B"), root.edge = 3) # = (C:1,(D:1,E:1):1):3;
}
\keyword{manip}
