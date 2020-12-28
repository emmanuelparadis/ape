\name{woodmouse}
\alias{woodmouse}
\title{Cytochrome b Gene Sequences of Woodmice}
\description{
  This is a set of 15 sequences of the mitochondrial gene cytochrome
  \emph{b} of the woodmouse (\emph{Apodemus sylvaticus}) which is a
  subset of the data analysed by Michaux et al. (2003). The full data
  set is available through GenBank (accession numbers AJ511877 to
  AJ511987).
}
\usage{
data(woodmouse)
}
\format{
  An object of class \code{"DNAbin"}.
}
\source{
  Michaux, J. R., Magnanou, E., Paradis, E., Nieberding, C. and Libois,
  R. (2003) Mitochondrial phylogeography of the Woodmouse
  (\emph{Apodemus sylvaticus}) in the Western Palearctic region.
  \emph{Molecular Ecology}, \bold{12}, 685--697.
}
\seealso{
  \code{\link{read.dna}}, \code{\link{DNAbin}}, \code{\link{dist.dna}}
}
\examples{
data(woodmouse)
str(woodmouse)
}
\keyword{datasets}
