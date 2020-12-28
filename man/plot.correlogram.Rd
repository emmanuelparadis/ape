\name{plot.correlogram}
\alias{plot.correlogram}
\alias{plot.correlogramList}
\title{Plot a Correlogram}
\usage{
  \method{plot}{correlogram}(x, legend = TRUE, test.level = 0.05,
                col = c("grey", "red"), type = "b", xlab = "",
                ylab = "Moran's I", pch = 21, cex = 2, ...)
  \method{plot}{correlogramList}(x, lattice = TRUE, legend = TRUE,
                test.level = 0.05, col = c("grey", "red"),
                xlab = "", ylab = "Moran's I",
                type = "b", pch = 21, cex = 2, ...)
}
\arguments{
  \item{x}{an object of class \code{"correlogram"} or of class
    \code{"correlogramList"} (both produced by
    \code{\link{correlogram.formula}}).}
  \item{legend}{should a legend be added on the plot?}
  \item{test.level}{the level used to discriminate the plotting symbols
    with colours considering the P-values.}
  \item{col}{two colours for the plotting symbols: the first one is used
    if the P-value is greater than or equal to \code{test.level}, the
    second one otherwise.}
  \item{type}{the type of plot to produce (see
    \code{\link[graphics]{plot}} for possible choices).}
  \item{xlab}{an optional character string for the label on the x-axis
    (none by default).}
  \item{ylab}{the default label on the y-axis.}
  \item{pch}{the type of plotting symbol.}
  \item{cex}{the default size for the plotting symbols.}
  \item{lattice}{when plotting several correlograms, should they be
    plotted in trellis-style with lattice (the default), or together on
    the same plot?}
  \item{\dots}{other parameters passed to the \code{plot} or \code{lines}
    function.}
}
\description{
  These functions plot correlagrams previously computed with
  \code{\link{correlogram.formula}}.
}
\details{
  When plotting several correlograms with lattice, some options have no
  effect: \code{legend}, \code{type}, and \code{pch} (\code{pch=19} is
  always used in this situation).

  When using \code{pch} between 1 and 20 (i.e., non-filled symbols, the
  colours specified in \code{col} are also used for the lines joining
  the points. To keep black lines, it is better to leave \code{pch}
  between 21 and 25.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{correlogram.formula}}, \code{\link{Moran.I}}
}
\keyword{hplot}
