\name{makeLabel}
\alias{makeLabel}
\alias{makeLabel.character}
\alias{makeLabel.phylo}
\alias{makeLabel.multiPhylo}
\alias{makeLabel.DNAbin}
\title{Label Management}
\description{
  This is a generic function with methods for character vectors, trees
  of class \code{"phylo"}, lists of trees of class \code{"multiPhylo"},
  and DNA sequences of class \code{"DNAbin"}. All options for the class
  character may be used in the other methods.
}
\usage{
makeLabel(x, ...)
\method{makeLabel}{character}(x, len = 99, space = "_", make.unique = TRUE,
          illegal = "():;,[]", quote = FALSE, ...)
\method{makeLabel}{phylo}(x, tips = TRUE, nodes = TRUE, ...)
\method{makeLabel}{multiPhylo}(x, tips = TRUE, nodes = TRUE, ...)
\method{makeLabel}{DNAbin}(x, ...)
}
\arguments{
  \item{x}{a vector of mode character or an object for which labels are
    to be changed.}
  \item{len}{the maximum length of the labels: those longer than `len'
    will be truncated.}
  \item{space}{the character to replace spaces, tabulations, and
    linebreaks.}
  \item{make.unique}{a logical specifying whether duplicate labels
    should be made unique by appending numerals; \code{TRUE} by
    default.}
  \item{illegal}{a string specifying the characters to be deleted.}
  \item{quote}{a logical specifying whether to quote the labels;
    \code{FALSE} by default.}
  \item{tips}{a logical specifying whether tip labels are to be
    modified; \code{TRUE} by default.}
  \item{nodes}{a logical specifying whether node labels are to be
    modified; \code{TRUE} by default.}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
  The option \code{make.unique} does not work exactly in the same way
  then the function of the same name: numbers are suffixed to all labels
  that are identical (without separator). See the examples.

  If there are 10--99 identical labels, the labels returned are "xxx01",
  "xxx02", etc, or "xxx001", "xxx002", etc, if they are 100--999, and so
  on. The number of digits added preserves the option `len'.

  The default for `len' makes labels short enough to be read by
  PhyML. Clustal accepts labels up to 30 character long.
}
\note{
  The current version does not perform well when trying to make very
  short unique labels (e.g., less than 5 character long).
}
\value{
  An object of the appropriate class.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{makeNodeLabel}}, \code{\link[base]{make.unique}},
  \code{\link[base]{make.names}}, \code{\link[base]{abbreviate}},
  \code{\link{mixedFontLabel}}, \code{\link{label2table}},
  \code{\link{updateLabel}}, \code{\link{checkLabel}}
}
\examples{
x <- rep("a", 3)
makeLabel(x)
make.unique(x) # <- from R's base
x <- rep("aaaaa", 2)
makeLabel(x, len = 3) # made unique and of length 3
makeLabel(x, len = 3, make.unique = FALSE)
}
\keyword{manip}
