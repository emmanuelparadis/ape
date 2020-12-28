\name{mvr}
\alias{mvr}
\alias{mvrs}
\title{Minimum Variance Reduction}
\description{
  Phylogenetic tree construction based on the minimum variance reduction.
}
\usage{
mvr(X, V)
mvrs(X, V, fs = 15)
}
\arguments{
  \item{X}{a distance matrix.}
  \item{V}{a variance matrix.}
  \item{fs}{agglomeration criterion parameter: it is coerced as an
    integer and must at least equal to one.}
}
\details{
  The MVR method can be seen as a version of BIONJ which is not
  restricted to the Poison model of variance (Gascuel 2000).
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  Criscuolo, A. and Gascuel, O. (2008). Fast NJ-like algorithms to deal
  with incomplete distance matrices. \emph{BMC Bioinformatics}, 9.

  Gascuel, O. (2000). Data model and classification by trees: the
  minimum variance reduction (MVR) method. \emph{Journal of
    Classification}, \bold{17}, 67--99.
}
\author{Andrei Popescu}
\seealso{
  \code{\link{bionj}}, \code{\link{fastme}}, \code{\link{njs}},
  \code{\link{SDM}}
}
\examples{
data(woodmouse)
rt <- dist.dna(woodmouse, variance = TRUE)
v <- attr(rt, "variance")
tr <- mvr(rt, v)
plot(tr, "u")
}
\keyword{models}
