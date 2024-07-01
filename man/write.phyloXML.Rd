\name{write.phyloXML}
\alias{write.phyloXML}
\title{Write Tree File in phyloXML Format}
\usage{
write.phyloXML(phy, file = "", tree.names = FALSE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{file = ""} (the default) then the
    tree is written on the standard output connection (i.e. the console).}
  \item{tree.names}{either a logical or a vector of mode character
    specifying whether or which tree names should be written to the file.}
}
\description{
  This function writes trees to a file of phyloXML format.
}
\details{
  If several trees are given, they will be represented as multiple
  <phylogeny> elements. Contrary to \code{\link{write.nexus}}, the trees
  need not have the same tip labels.

  When \code{tree.names} is \code{TRUE}, the tree names will be always
  added as <name> tags to each phylogeny element. If the \code{phy}
  object is unnamed, then the names will be automatically generated
  from the tree indices as "tree<index>" (e.g. tree1, tree2, ...). If
  \code{tree.names} is a character vector, the specified names will be
  used instead.

  Branch lengths, labels, and rootedness are preserved in the phyloXML
  file.
}
\value{
  None (invisible \code{NULL}).
}
\references{
  Han, M.V., Zmasek, C.M. (2009) phyloXML: XML for evolutionary biology and comparative genomics. \emph{BMC Bioinformatics}, \bold{10}, 356.
}

\author{Emmanuel Paradis and Federico Marotta}
\seealso{
  \code{\link{read.tree}}, \code{\link{write.tree}},
  \code{\link{read.nexus}}, \code{\link{write.nexus}},
  \code{\link{read.nexus.data}}, \code{\link{write.nexus.data}}
}
\keyword{manip}
\keyword{IO}
