\name{njs}
\alias{njs}
\alias{bionjs}
\title{Tree Reconstruction from Incomplete Distances With NJ* or bio-NJ*}
\description{
  Reconstructs a phylogenetic tree from a distance matrix with possibly
  missing values.
}
\usage{
njs(X, fs = 15)
bionjs(X, fs = 15)
}
\arguments{
  \item{X}{a distance matrix.}
  \item{fs}{argument \emph{s} of the agglomerative criterion: it is
    coerced as an integer and must at least equal to one.}
}
\details{
  Missing values represented by either \code{NA} or any negative number.

  Basically, the Q* criterion is applied to all the pairs of leaves, and
  the \emph{s} highest scoring ones are chosen for further analysis by
  the agglomeration criteria that better handle missing distances (see
  references for details).
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  Criscuolo, A., Gascuel, O. (2008) Fast NJ-like algorithms to deal with
  incomplete distance matrices. \emph{BMC Bioinformatics}, \bold{9},
  166.
}
\author{Andrei Popescu}
\seealso{
  \code{\link{nj}}, \code{\link{bionj}}, \code{\link{triangMtds}}
}
\examples{
data(woodmouse)
d <- dist.dna(woodmouse)
dm <- d
dm[sample(length(dm), size = 3)] <- NA
dist.topo(njs(dm), nj(d)) # often 0
dm[sample(length(dm), size = 10)] <- NA
dist.topo(njs(dm), nj(d)) # sometimes 0
}
\keyword{models}
