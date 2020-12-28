\name{write.dna}
\alias{write.dna}
\alias{write.FASTA}
\title{Write DNA Sequences in a File}
\usage{
write.dna(x, file, format = "interleaved", append = FALSE,
          nbcol = 6, colsep = " ", colw = 10, indent = NULL,
          blocksep = 1)
write.FASTA(x, file, header = NULL, append = FALSE)
}
\arguments{
  \item{x}{a list or a matrix of DNA sequences, or of AA sequences for
    \code{write.FASTA}.}
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string.}
  \item{format}{a character string specifying the format of the DNA
    sequences. Three choices are possible: \code{"interleaved"},
    \code{"sequential"}, or \code{"fasta"}, or any unambiguous
    abbreviation of these.}
  \item{append}{a logical, if \code{TRUE} the data are appended to the
    file without erasing the data possibly existing in the file,
    otherwise the file (if it exists) is overwritten (\code{FALSE} the
    default).}
  \item{nbcol}{a numeric specifying the number of columns per row (6 by
    default); may be negative implying that the nucleotides are printed
    on a single line.}
  \item{colsep}{a character used to separate the columns (a single
    space by default).}
  \item{colw}{a numeric specifying the number of nucleotides per column
    (10 by default).}
  \item{indent}{a numeric or a character specifying how the blocks of
    nucleotides are indented (see details).}
  \item{blocksep}{a numeric specifying the number of lines between the
    blocks of nucleotides (this has an effect only if `format =
    "interleaved"').}
  \item{header}{a vector of mode character giving the header to be
    written in the FASTA file before the sequences. By default, there is
    no header.}
}
\description{
  These functions write in a file a list of DNA sequences in sequential,
  interleaved, or FASTA format. \code{write.FASTA} can write either DNA
  or AA sequences.
}
\details{
  Three formats are supported in the present function: see the help page
  of \code{\link{read.dna}} and the references below for a description.

  If the sequences have no names, then they are given "1", "2", ... as
  labels in the file.

  With the interleaved and sequential formats, the sequences must be all
  of the same length. The names of the sequences are not truncated.

  The argument \code{indent} specifies how the rows of nucleotides are
  indented. In the interleaved and sequential formats, the rows with
  the taxon names are never indented; the subsequent rows are indented
  with 10 spaces by default (i.e., if \code{indent = NULL}). In the FASTA
  format, the rows are not indented by default. This default behaviour
  can be modified by specifying a value to \code{indent}: the rows are then
  indented with ``indent'' (if it is a character) or `indent' spaces (if
  it is a numeric). For example, specifying \code{indent = "   "} or
  \code{indent = 3} will have the same effect (use \code{indent = "\\t"}
  for a tabulation).

  The different options are intended to give flexibility in formatting
  the sequences. For instance, if the sequences are very long it may be
  judicious to remove all the spaces beween columns (colsep = ""), in
  the margins (indent = 0), and between the blocks (blocksep = 0) to
  produce a smaller file.

  \code{write.dna(, format = "fasta")} can be very slow if the sequences
  are long (> 10 kb). \code{write.FASTA} is much faster in this
  situation but the formatting is not flexible: each sequence is printed
  on a single line, which is OK for big files that are not intended to
  be open with a text editor.
}
\note{
  Specifying a negative value for `nbcol' (meaning that the nucleotides
  are printed on a single line) gives the same output for the
  interleaved and sequential formats.

  The names of the sequences can be truncated with the function
  \code{\link{makeLabel}}. In particular, Clustal is limited to 30
  characters, and PHYML seems limited to 99 characters.
}
\value{
  None (invisible `NULL').
}
\author{Emmanuel Paradis}
\references{
  Anonymous. FASTA format. \url{https://en.wikipedia.org/wiki/FASTA_format}

  Felsenstein, J. (1993) Phylip (Phylogeny Inference Package) version
  3.5c. Department of Genetics, University of Washington.
  \url{http://evolution.genetics.washington.edu/phylip/phylip.html}
}
\seealso{
  \code{\link{read.dna}}, \code{\link{read.GenBank}},
  \code{\link{makeLabel}}
}
\keyword{IO}
