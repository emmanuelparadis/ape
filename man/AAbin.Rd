\name{AAbin}
\alias{AAbin}
\alias{print.AAbin}
\alias{[.AAbin}
\alias{as.character.AAbin}
\alias{labels.AAbin}
\alias{image.AAbin}
\alias{as.AAbin}
\alias{as.AAbin.AAString}
\alias{as.AAbin.AAStringSet}
\alias{as.AAbin.AAMultipleAlignment}
\alias{as.AAbin.character}
\alias{as.phyDat.AAbin}
\alias{dist.aa}
\alias{AAsubst}
\alias{c.AAbin}
\alias{cbind.AAbin}
\alias{rbind.AAbin}
\alias{as.AAbin.list}
\alias{as.list.AAbin}
\alias{as.matrix.AAbin}
\title{Amino Acid Sequences}
\description{
  These functions help to create and manipulate AA sequences.
}
\usage{
\method{print}{AAbin}(x, \dots)

\method{[}{AAbin}(x, i, j, drop = FALSE)

\method{c}{AAbin}(..., recursive = FALSE)

\method{rbind}{AAbin}(\dots)
\method{cbind}{AAbin}(\dots, check.names = TRUE, fill.with.Xs = FALSE,
             quiet = FALSE)

\method{as.character}{AAbin}(x, \dots)

\method{labels}{AAbin}(object, \dots)

\method{image}{AAbin}(x, what, col, bg = "white", xlab = "", ylab = "",
      show.labels = TRUE, cex.lab = 1, legend = TRUE, grid = FALSE,
      show.aa = FALSE, aa.cex = 1, aa.font = 1, aa.col = "black",\dots)

as.AAbin(x, \dots)
\method{as.AAbin}{character}(x, \dots)

\method{as.AAbin}{list}(x, ...)

\method{as.AAbin}{AAString}(x, ...)

\method{as.AAbin}{AAStringSet}(x, ...)

\method{as.AAbin}{AAMultipleAlignment}(x, ...)

\method{as.list}{AAbin}(x, ...)

\method{as.matrix}{AAbin}(x, ...)

\method{as.phyDat}{AAbin}(x, \dots)

dist.aa(x, pairwise.deletion = FALSE, scaled = FALSE)
AAsubst(x)
}
\arguments{
\item{x, object}{an object of class \code{"AAbin"} (or else depending on
the function).}
\item{i, j}{indices of the rows and/or columns to select or to
drop. They may be numeric, logical, or character (in the same way than
for standard \R objects).}
\item{drop}{logical; if \code{TRUE}, the returned object is of the
lowest possible dimension.}
\item{recursive}{logical; whether to go down lists and concatenate its
elements.}
  \item{check.names}{a logical specifying whether to check the rownames
    before binding the columns (see details).}
  \item{fill.with.Xs}{a logical indicating whether to keep all
    possible individuals as indicating by the rownames, and eventually
    filling the missing data with insertion gaps (ignored if
    \code{check.names = FALSE}).}
  \item{quiet}{a logical to switch off warning messages when some rows
    are dropped.}
\item{what}{a vector of characters specifying the amino acids to
visualize. Currently, the only possible choice is to show the three
categories hydrophobic, small, and hydrophilic.}
\item{col}{a vector of colours. If missing, this is set to ``red'',
``yellow'' and ``blue''.}
\item{bg}{the colour used for AA codes not among \code{what} (typically
X and *).}
\item{xlab}{the label for the \emph{x}-axis; none by default.}
\item{ylab}{Idem for the \emph{y}-axis. Note that by default, the labels
of the sequences are printed on the \emph{y}-axis (see next option).}
\item{show.labels}{a logical controlling whether the sequence labels are
printed (\code{TRUE} by default).}
\item{cex.lab}{a single numeric controlling the size of the sequence
labels. Use \code{cex.axis} to control the size of the annotations on
the \emph{x}-axis.}
\item{legend}{a logical controlling whether the legend is plotted
(\code{TRUE} by default).}
\item{grid}{a logical controlling whether to draw a grid (\code{FALSE} by
default).}
\item{show.aa}{a logical controlling whether to show the AA symbols (\code{FALSE} by
default).}
\item{aa.cex, aa.font, aa.col}{control the aspect of the AA symbols
(ignored if the previous is \code{FALSE}).}
\item{pairwise.deletion}{a logical indicating whether to delete the
sites with missing data in a pairwise way. The default is to delete the
sites with at least one missing data for all sequences.}
\item{scaled}{a logical value specifying whether to scale the number of
AA differences by the sequence length.}
\item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
These functions help to manipulate amino acid sequences of class
\code{"AAbin"}. These objects are stored in vectors, matrices, or lists
which can be manipulated with the usual \code{[} operator.

There is a conversion function to and from characters.

The function \code{dist.aa} computes the number of AA differences
between each pair of sequences in a matrix; this can be scaled by the
sequence length. See the function \code{\link[phangorn]{dist.ml}} in
\pkg{phangorn} for evolutionary distances with AA sequences.

The function \code{AAsubst} returns the indices of the polymorphic sites
(similar to \code{\link{seg.sites}} for DNA sequences; see examples below).

The two functions \code{cbind.AAbin} and \code{rbind.AAbin} work in the
same way than the similar methods for the class \code{"DNAbin"}: see
\code{\link{cbind.DNAbin}} for more explanations about their respective
behaviours.
}
\value{
an object of class \code{"AAbin"}, \code{"character"},
\code{"dist"}, or \code{"numeric"}, depending on the function.
}
\author{Emmanuel Paradis, Franz Krah}
\seealso{
\code{\link{read.FASTA}}, \code{\link{trans}}, \code{\link{alview}}
}
\examples{
data(woodmouse)
AA <- trans(woodmouse, 2)
seg.sites(woodmouse)
AAsubst(AA)
}
\keyword{manip}
