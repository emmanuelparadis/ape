\name{delta.plot}
\alias{delta.plot}
\title{Delta Plots}
\usage{
delta.plot(X, k = 20, plot = TRUE, which = 1:2)
}
\arguments{
  \item{X}{a distance matrix, may be an object of class ``dist''.}
  \item{k}{an integer giving the number of intervals in the plot.}
  \item{plot}{a logical specifying whether to draw the
    \eqn{\delta}{delta} plot (the default).}
  \item{which}{a numeric vector indicating which plots are done; 1: the
    histogram of the \eqn{\delta_q}{delta_q} values, 2: the plot of the
    individual \eqn{\bar{\delta}}{delta.bar} values. By default, both
    plots are done.}
}
\description{
  This function makes a \eqn{\delta}{delta} plot following Holland et
  al. (2002).
}
\details{
  See Holland et al. (2002) for details and interpretation.

  The computing time of this function is proportional to the fourth
  power of the number of observations (\eqn{O(n^4)}), so calculations
  may be very long with only a slight increase in sample size.
}
\value{
  This function returns invisibly a named list with two components:

  \itemize{
    \item{counts}{the counts for the histogram of
    \eqn{\delta_q}{delta_q} values}
    \item{delta.bar}{the mean \eqn{\delta}{delta} value for each
    observation}
  }
}
\references{
  Holland, B. R., Huber, K. T., Dress, A. and Moulton, V. (2002) Delta
  plots: a tool for analyzing phylogenetic distance data.
  \emph{Molecular Biology and Evolution}, \bold{12}, 2051--2059.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{dist.dna}}
}
\examples{
data(woodmouse)
d <- dist.dna(woodmouse)
delta.plot(d)
layout(1)
delta.plot(d, 40, which = 1)
}
\keyword{hplot}
