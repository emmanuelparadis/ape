\name{chiroptera}
\alias{chiroptera}
\title{Bat Phylogeny}
\description{
  This phylogeny of bats (Mammalia: Chiroptera) is a supertree (i.e. a
  composite phylogeny constructed from several sources; see source for
  details).
}
\usage{
data(chiroptera)
}
\format{
  The data are stored in RData (binary) format.
}
\source{
  Jones, K. E., Purvis, A., MacLarnon, A., Bininda-Emonds, O. R. P. and
  Simmons, N. B. (2002) A phylogenetic supertree of the bats (Mammalia:
  Chiroptera). \emph{Biological Reviews of the Cambridge Philosophical
    Society}, \bold{77}, 223--259.
}
\seealso{
  \code{\link{read.nexus}}, \code{\link{zoom}}
}
\examples{
data(chiroptera)
str(chiroptera)
op <- par()
par(cex = 0.3)
plot(chiroptera, type = "c")
par(cex = op$cex)
}
\keyword{datasets}
