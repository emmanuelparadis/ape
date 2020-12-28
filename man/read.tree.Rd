\name{read.tree}
\alias{read.tree}
\alias{phylo}
\title{Read Tree File in Parenthetic Format}
\usage{
read.tree(file = "", text = NULL, tree.names = NULL, skip = 0,
    comment.char = "", keep.multi = FALSE, ...)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{file = ""} (the default) then the
    tree is input on the keyboard, the entry being terminated with a
    blank line.}
  \item{text}{alternatively, the name of a variable of mode character
    which contains the tree(s) in parenthetic format. By default, this
    is ignored (set to \code{NULL}, meaning that the tree is read in a
    file); if \code{text} is not \code{NULL}, then the argument
    \code{file} is ignored.}
  \item{tree.names}{if there are several trees to be read, a vector of
    mode character that gives names to the individual trees; if
    \code{NULL} (the default), the trees are named \code{"tree1"},
    \code{"tree2"}, ...}
  \item{skip}{the number of lines of the input file to skip before
    beginning to read data (this is passed directly to\code{ scan()}).}
  \item{comment.char}{a single character, the remaining of the line
    after this character is ignored (this is passed directly to
    \code{scan()}).}
  \item{keep.multi}{if \code{TRUE} and \code{tree.names = NULL} then
    single trees are returned in \code{"multiPhylo"} format, with any
    name that is present (see details). Default is \code{FALSE}.}
  \item{\dots}{further arguments to be passed to \code{scan()}.}
}
\description{
  This function reads a file which contains one or several trees in
  parenthetic format known as the Newick or New Hampshire format.
}
\details{
  The default option for \code{file} allows to type directly the tree on
  the keyboard (or possibly to copy from an editor and paste in R's
  console) with, e.g., \code{mytree <- read.tree()}.

  `read.tree' tries to represent correctly trees with a badly
  represented root edge (i.e. with an extra pair of parentheses). For
  instance, the tree "((A:1,B:1):10);" will be read like "(A:1,B:1):10;"
  but a warning message will be issued in the former case as this is
  apparently not a valid Newick format. If there are two root edges
  (e.g., "(((A:1,B:1):10):10);"), then the tree is not read and an error
  message is issued.

  If there are any characters preceding the first "(" in a line then
  this is assigned to the name. This is returned when a "multiPhylo"
  object is returned and \code{tree.names = NULL}.

  Until \pkg{ape} 4.1, the default of \code{comment.char} was \code{"#"}
  (as in \code{scan}). This has been changed so that extended Newick
  files can be read.
}
\value{
  an object of class \code{"phylo"} with the following components:
  \item{edge}{a two-column matrix of mode numeric where each row
    represents an edge of the tree; the nodes and the tips are
    symbolized with numbers; the tips are numbered 1, 2, \dots, and the
    nodes are numbered after the tips. For each row, the first column
    gives the ancestor.}
  \item{edge.length}{(optional) a numeric vector giving the lengths of the
    branches given by \code{edge}.}
  \item{tip.label}{a vector of mode character giving the names of the
    tips; the order of the names in this vector corresponds to the
    (positive) number in \code{edge}.}
  \item{Nnode}{the number of (internal) nodes.}
  \item{node.label}{(optional) a vector of mode character giving the
    names of the nodes.}
  \item{root.edge}{(optional) a numeric value giving the length of the
    branch at the root if it exists.}

  If several trees are read in the file, the returned object is of class
  \code{"multiPhylo"}, and is a list of objects of class \code{"phylo"}.
  The name of each tree can be specified by \code{tree.names}, or can be
  read from the file (see details).
}
\references{
  Felsenstein, J. The Newick tree format.
  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}

  Olsen, G. Interpretation of the "Newick's 8:45" tree format standard.
  \url{http://evolution.genetics.washington.edu/phylip/newick_doc.html}

  Paradis, E. (2008) Definition of Formats for Coding Phylogenetic Trees
  in R. \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}

  Paradis, E. (2012) \emph{Analysis of Phylogenetics and Evolution with
    R (Second Edition).} New York: Springer.
}

\author{Emmanuel Paradis and Daniel Lawson \email{dan.lawson@bristol.ac.uk}}
\seealso{
  \code{\link{write.tree}}, \code{\link{read.nexus}},
  \code{\link{write.nexus}}, \code{\link[base]{scan}} for the basic R
  function to read data in a file
}
\examples{
### An extract from Sibley and Ahlquist (1990)
s <- "owls(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
cat(s, file = "ex.tre", sep = "\n")
tree.owls <- read.tree("ex.tre")
str(tree.owls)
tree.owls
tree.owls <- read.tree("ex.tre", keep.multi = TRUE)
tree.owls
names(tree.owls)
unlink("ex.tre") # delete the file "ex.tre"
### Only the first three species using the option `text'
TREE <- "((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3);"
TREE
tree.owls.bis <- read.tree(text = TREE)
str(tree.owls.bis)
tree.owls.bis

## tree with singleton nodes:
ts <- read.tree(text="((((a))),d);")
plot(ts, node.depth = 2) # the default will overlap the singleton node with the tip
nodelabels()

## skeleton tree with a singleton node:
tx <- read.tree(text = "(((,)),);")
plot(tx, node.depth = 2)
nodelabels()

## a tree with single quoted labels (the 2nd label is not quoted):
z <- "(('a: France, Spain (Europe)',b),'c: Australia [Outgroup]');"
tz <- read.tree(text = z)
plot(tz, font = 1)
}
\keyword{manip}
\keyword{IO}
