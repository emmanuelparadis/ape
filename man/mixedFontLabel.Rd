\name{mixedFontLabel}
\alias{mixedFontLabel}
\title{Mixed Font Labels for Plotting}
\description{
  This function helps to format labels with bits of text in different
  font shapes (italics, bold, or bolditalics) and different
  separators. The output is intended to be used for plotting.
}
\usage{
mixedFontLabel(..., sep = " ", italic = NULL, bold = NULL,
               parenthesis = NULL,
               always.upright = c("sp.", "spp.", "ssp."))
}
\arguments{
  \item{\dots}{vectors of mode character to be formatted. They may be
    of different lengths in which case the shortest ones are
    recycled.}
  \item{sep}{a vector of mode character giving the separators to be
    printed between the elements in \code{\dots}.}
  \item{italic}{a vector of integers specifying the elements in
    \code{\dots} to be printed in italics.}
  \item{bold}{id. in boldface.}
  \item{parenthesis}{id. within parentheses.}
  \item{always.upright}{of vector of mode character giving the strings
    to not print in italics. Use \code{always.upright = ""} to cancel
    this option.}
}
\details{
  The idea is to have different bits of text in different vectors that
  are put together to make a vector of \R expressions. This vector is
  interpreted by graphical functions to format the text. A simple use
  may be \code{mixedFontLabel(genus, species, italic = 1:2)}, but it is
  more interesting when mixing fonts (see examples).

  To have an element in bolditalics, its number must given in both
  \code{italic} and \code{bold}.

  The vector returned by this function may be assigned as the
  \code{tip.label} element of a tree of class \code{"phylo"}, or even as
  its \code{node.label} element.
}
\value{
  A vector of mode expression.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{makeLabel}}, \code{\link{makeNodeLabel}},
  \code{\link{label2table}}, \code{\link{updateLabel}},
  \code{\link{checkLabel}}
}
\examples{
tr <- read.tree(text = "((a,(b,c)),d);")
genus <- c("Gorilla", "Pan", "Homo", "Pongo")
species <- c("gorilla", "spp.", "sapiens", "pygmaeus")
geo <- c("Africa", "Africa", "World", "Asia")
tr$tip.label <- mixedFontLabel(genus, species, geo, italic = 1:2,
  parenthesis = 3)
layout(matrix(c(1, 2), 2))
plot(tr)
tr$tip.label <- mixedFontLabel(genus, species, geo, sep = c(" ", " | "),
  italic = 1:2, bold = 3)
plot(tr)
layout(1)
}
\keyword{manip}
