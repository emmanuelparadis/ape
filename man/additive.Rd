\name{additive}
\alias{additive}
\alias{ultrametric}
\title{Incomplete Distance Matrix Filling}
\description{
  Fills missing entries from incomplete distance matrix using the
  additive or the ultrametric procedure (see reference for details).
}
\usage{
additive(X)
ultrametric(X)
}
\arguments{
  \item{X}{a distance matrix or an object of class \code{"dist"}.}
}
\value{
  a distance matrix.
}
\references{
  Makarenkov, V. and Lapointe, F.-J. (2004) A weighted least-squares
  approach for inferring phylogenies from incomplete distance
  matrices. \emph{Bioinformatics}, \bold{20}, 2113--2121.
}
\author{Andrei Popescu}
\keyword{manip}
