\name{plot.varcomp}
\alias{plot.varcomp}
\title{Plot Variance Components}
\description{
  Plot previously estimated variance components.
}
\usage{
\method{plot}{varcomp}(x, xlab = "Levels", ylab = "Variance", type = "b", ...)
}
\arguments{
  \item{x}{ A \var{varcomp} object}
  \item{xlab}{ x axis label}
  \item{ylab}{ y axis label }
  \item{type}{ plot type ("l", "p" or "b", see \code{\link{plot}})}
  \item{\dots}{Further argument sent to the \code{\link[lattice]{xyplot}} function.}
}
\value{
  The same as \code{\link[lattice]{xyplot}}.
}
\author{Julien Dutheil \email{dutheil@evolbio.mpg.de}}
\seealso{\code{\link{varcomp}}}
\keyword{hplot}
