\name{as.phylo.formula}
\alias{as.phylo.formula}
\title{Conversion from Taxonomy Variables to Phylogenetic Trees}
\description{
  The function \code{as.phylo.formula} (short form \code{as.phylo})
  builds a phylogenetic tree (an object of class \code{phylo}) from
  a set of nested taxonomic variables.
}
\usage{
\method{as.phylo}{formula}(x, data = parent.frame(), collapse = TRUE, ...)
}
\arguments{
  \item{x}{a right-side formula describing the taxonomic relationship:
    \code{~C1/C2/.../Cn}.}
  \item{data}{the data.frame where to look for the variables (default
    to user's workspace).}
  \item{collapse}{a logical value specifying whether to collapse single
    nodes in the returned tree (see details).}
  \item{\dots}{further arguments to be passed from other methods.}
}
\details{
  Taxonomic variables must be nested and passed in the correct order:
  the higher clade must be on the left of the formula, for instance
  \code{~Order/Family/Genus/Species}. In most cases, the resulting tree
  will be unresolved and will contain polytomies.

  The option \code{collapse = FALSE} has for effect to add single nodes
  in the tree when a given higher level has only one element in the
  level below (e.g., a monospecific genus); see the example below.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Julien Dutheil \email{dutheil@evolbio.mpg.de} and Eric
  Marcon}
\seealso{
  \code{\link{as.phylo}}, \code{\link{read.tree}} for a
  description of \code{"phylo"} objects, \code{\link{multi2di}}
}
\examples{
data(carnivora)
frm <- ~SuperFamily/Family/Genus/Species
tr <- as.phylo(frm, data = carnivora)
plot(tr)
Nnode(tr)
## compare with:
Nnode(as.phylo(frm, data = carnivora, collapse = FALSE))
}
\keyword{manip}
