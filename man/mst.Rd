\name{mst}
\alias{mst}
\alias{plot.mst}
\title{Minimum Spanning Tree}
\usage{
mst(X)
\method{plot}{mst}(x, graph = "circle", x1 = NULL, x2 = NULL, \dots)
}
\arguments{
  \item{X}{either a matrix that can be interpreted as a distance matrix,
    or an object of class \code{"dist"}.}
  \item{x}{an object of class \code{"mst"} (e.g. returned by \code{mst()}).}
  \item{graph}{a character string indicating the type of graph to plot
    the minimum spanning tree; two choices are possible: \code{"circle"} where
    the observations are plotted regularly spaced on a circle, and
    \code{"nsca"} where the two first axes of a non-symmetric correspondence
    analysis are used to plot the observations (see Details below). If
    both arguments \code{x1} and \code{x2} are given, the argument
    \code{graph} is ignored.}
  \item{x1}{a numeric vector giving the coordinates of the observations
    on the \emph{x}-axis. Both \code{x1} and \code{x2} must be specified
    to be used.}
  \item{x2}{a numeric vector giving the coordinates of the observations
    on the \emph{y}-axis. Both \code{x1} and \code{x2} must be specified
    to be used.}
  \item{\dots}{further arguments to be passed to \code{plot()}.}
}
\description{
  The function \code{mst} finds the minimum spanning tree between a set of
  observations using a matrix of pairwise distances.

  The \code{plot} method plots the minimum spanning tree showing the
  links where the observations are identified by their numbers.
}
\details{
  These functions provide two ways to plot the minimum spanning tree which
  try to space as much as possible the observations in order to show as
  clearly as possible the links. The option \code{graph = "circle"}
  simply plots regularly the observations on a circle, whereas
  \code{graph = "nsca"} uses a non-symmetric correspondence analysis
  where each observation is represented at the centroid of its neighbours.

  Alternatively, the user may use any system of coordinates for the
  obsevations, for instance a principal components analysis (PCA) if the
  distances were computed from an original matrix of continous variables.
}
\value{
  an object of class \code{"mst"} which is a square numeric matrix of size
  equal to the number of observations with either \code{1} if a link
  between the corresponding observations was found, or \code{0}
  otherwise. The names of the rows  and columns of the distance matrix,
  if available, are given as rownames and colnames to the returned object.
}
\author{
  Yvonnick Noel \email{noel@univ-lille3.fr},
  Julien Claude \email{julien.claude@umontpellier.fr} and
  Emmanuel Paradis
}
\seealso{
  \code{\link{dist.dna}}, \code{\link{dist.gene}},
  \code{\link[stats]{dist}}, \code{\link[graphics]{plot}}
}
\examples{
require(stats)
X <- matrix(runif(200), 20, 10)
d <- dist(X)
PC <- prcomp(X)
M <- mst(d)
opar <- par()
par(mfcol = c(2, 2))
plot(M)
plot(M, graph = "nsca")
plot(M, x1 = PC$x[, 1], x2 = PC$x[, 2])
par(opar)
}
\keyword{multivariate}
