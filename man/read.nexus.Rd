\name{read.nexus}
\alias{read.nexus}
\title{Read Tree File in Nexus Format}
\usage{
read.nexus(file, tree.names = NULL, force.multi = FALSE)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string.}
  \item{tree.names}{if there are several trees to be read, a vector of
    mode character giving names to the individual trees (by default,
    this uses the labels in the NEXUS file if these are present).}
  \item{force.multi}{a logical value; if \code{TRUE}, an object of class
    \code{"multiPhylo"} is always returned even if the file contains a
    single tree (see details).}
}
\description{
  This function reads one or several trees in a NEXUS file.
}
\details{
  The present implementation tries to follow as much as possible the
  NEXUS standard. Only the block ``TREES'' is read; the other data can
  be read with other functions (e.g., \code{\link{read.dna}},
  \code{\link[utils]{read.table}}, \dots).

  Until \pkg{ape} version 5.8-1, if a TRANSLATION table was present it
  was assumed that only the tip labels are translated with integers
  without gap and nodes labels were not looked for in the translation
  table. Recent versions of \pkg{ape} now conforms more closely to the
  NEXUS standard as described in Maddison \emph{et al} (1997). Note that
  \code{\link{write.nexus}} translates only the tip labels.

  Using \code{force.multi = TRUE} when the file contains a single tree
  makes possible to keep the tree name (as names of the list).

  `read.nexus' tries to represent correctly trees with a badly
  represented root edge (i.e. with an extra pair of parentheses). For
  instance, the tree "((A:1,B:1):10);" will be read like "(A:1,B:1):10;"
  but a warning message will be issued in the former case as this is
  apparently not a valid Newick format. If there are two root edges
  (e.g., "(((A:1,B:1):10):10);"), then the tree is not read and an error
  message is issued.
}
\value{
  an object of class \code{"phylo"} or \code{"multiPhylo"}.
}
\references{
  Maddison, D. R., Swofford, D. L. and Maddison, W. P. (1997) NEXUS: an
  extensible file format for systematic information. \emph{Systematic
    Biology}, \bold{46}, 590--621.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.tree}}, \code{\link{write.nexus}},
  \code{\link{write.tree}}, \code{\link{read.nexus.data}},
  \code{\link{write.nexus.data}}
}
\keyword{manip}
\keyword{IO}
