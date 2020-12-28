\name{DNAbin}
\alias{DNAbin}
\alias{print.DNAbin}
\alias{[.DNAbin}
\alias{rbind.DNAbin}
\alias{cbind.DNAbin}
\alias{as.matrix.DNAbin}
\alias{c.DNAbin}
\alias{as.list.DNAbin}
\alias{labels.DNAbin}
\title{Manipulate DNA Sequences in Bit-Level Format}
\description{
  These functions help to manipulate DNA sequences coded in the
  bit-level coding scheme.
}
\usage{
\method{print}{DNAbin}(x, printlen = 6, digits = 3, \dots)
\method{rbind}{DNAbin}(\dots)
\method{cbind}{DNAbin}(\dots, check.names = TRUE, fill.with.gaps = FALSE,
             quiet = FALSE)
\method{[}{DNAbin}(x, i, j, drop = FALSE)
\method{as.matrix}{DNAbin}(x, \dots)
\method{c}{DNAbin}(\dots, recursive = FALSE)
\method{as.list}{DNAbin}(x, \dots)
\method{labels}{DNAbin}(object, \dots)
}
\arguments{
  \item{x, object}{an object of class \code{"DNAbin"}.}
  \item{\dots}{either further arguments to be passed to or from other
    methods in the case of \code{print}, \code{as.matrix}, and
    \code{labels}, or a series of objects of class \code{"DNAbin"} in the
    case of \code{rbind}, \code{cbind}, and \code{c}.}
  \item{printlen}{the number of labels to print (6 by default).}
  \item{digits}{the number of digits to print (3 by default).}
  \item{check.names}{a logical specifying whether to check the rownames
    before binding the columns (see details).}
  \item{fill.with.gaps}{a logical indicating whether to keep all
    possible individuals as indicating by the rownames, and eventually
    filling the missing data with insertion gaps (ignored if
    \code{check.names = FALSE}).}
  \item{quiet}{a logical to switch off warning messages when some rows
    are dropped.}
  \item{i, j}{indices of the rows and/or columns to select or to drop.
    They may be numeric, logical, or character (in the same way than for
    standard \R objects).}
  \item{drop}{logical; if \code{TRUE}, the returned object is of the
    lowest possible dimension.}
  \item{recursive}{for compatibility with the generic (unused).}
}
\details{
  These are all `methods' of generic functions which are here applied to
  DNA sequences stored as objects of class \code{"DNAbin"}. They are
  used in the same way than the standard \R functions to manipulate
  vectors, matrices, and lists. Additionally, the operators \code{[[}
  and \code{$} may be used to extract a vector from a list. Note that
  the default of \code{drop} is not the same than the generic operator:
  this is to avoid dropping rownames when selecting a single sequence.

  These functions are provided to manipulate easily DNA sequences coded
  with the bit-level coding scheme. The latter allows much faster
  comparisons of sequences, as well as storing them in less memory
  compared to the format used before \pkg{ape} 1.10.

  For \code{cbind}, the default behaviour is to keep only individuals
  (as indicated by the rownames) for which there are no missing data. If
  \code{fill.with.gaps = TRUE}, a `complete' matrix is returned,
  enventually with insertion gaps as missing data. If \code{check.names
  = TRUE} (the default), the rownames of each matrix are checked, and
  the rows are reordered if necessary (if some rownames are duplicated,
  an error is returned). If \code{check.names = FALSE}, the matrices
  must all have the same number of rows, and are simply binded; the
  rownames of the first matrix are used. See the examples.

  \code{as.matrix} may be used to convert DNA sequences (of the same
  length) stored in a list into a matrix while keeping the names and the
  class. \code{as.list} does the reverse operation.
}
\value{
  an object of class \code{"DNAbin"} in the case of \code{rbind},
  \code{cbind}, and \code{[}.
}
\references{
  Paradis, E. (2007) A Bit-Level Coding Scheme for Nucleotides.
  \url{http://ape-package.ird.fr/misc/BitLevelCodingScheme_20April2007.pdf}

  Paradis, E. (2012) \emph{Analysis of Phylogenetics and Evolution with
  R (Second Edition).} New York: Springer.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{as.DNAbin}}, \code{\link{read.dna}},
  \code{\link{read.GenBank}}, \code{\link{write.dna}},
  \code{\link{image.DNAbin}},\code{\link{AAbin}}

  The corresponding generic functions are documented in the package
  \pkg{base}.
}
\examples{
data(woodmouse)
woodmouse
print(woodmouse, 15, 6)
print(woodmouse[1:5, 1:300], 15, 6)
### Just to show how distances could be influenced by sampling:
dist.dna(woodmouse[1:2, ])
dist.dna(woodmouse[1:3, ])
### cbind and its options:
x <- woodmouse[1:2, 1:5]
y <- woodmouse[2:4, 6:10]
as.character(cbind(x, y)) # gives warning
as.character(cbind(x, y, fill.with.gaps = TRUE))
\dontrun{
as.character(cbind(x, y, check.names = FALSE)) # gives an error
}
}
\keyword{manip}
