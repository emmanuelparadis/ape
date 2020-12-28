\name{read.dna}
\alias{read.dna}
\alias{read.FASTA}
\alias{read.fastq}
\title{Read DNA Sequences in a File}
\description{
  These functions read DNA sequences in a file, and returns a matrix or a
  list of DNA sequences with the names of the taxa read in the file as
  rownames or names, respectively. By default, the sequences are stored
  in binary format, otherwise (if \code{as.character = TRUE}) in
  lowercase.
}
\usage{
read.dna(file, format = "interleaved", skip = 0,
         nlines = 0, comment.char = "#",
         as.character = FALSE, as.matrix = NULL)
read.FASTA(file, type = "DNA")
read.fastq(file, offset = -33)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string. Can also be a \link{connection} (which
    will be opened for reading if necessary, and if so
    \code{\link{close}}d (and hence destroyed) at the end of the
    function call).}
  \item{format}{a character string specifying the format of the DNA
    sequences. Four choices are possible: \code{"interleaved"},
    \code{"sequential"}, \code{"clustal"}, or \code{"fasta"}, or any
    unambiguous abbreviation of these.}
  \item{skip}{the number of lines of the input file to skip before
    beginning to read data (ignored for FASTA files; see below).}
  \item{nlines}{the number of lines to be read (by default the file is
    read untill its end; ignored for FASTA files)).}
  \item{comment.char}{a single character, the remaining of the line
    after this character is ignored (ignored for FASTA files).}
  \item{as.character}{a logical controlling whether to return the
    sequences as an object of class \code{"DNAbin"} (the default).}
  \item{as.matrix}{(used if \code{format = "fasta"}) one of the three
    followings: (i) \code{NULL}: returns the sequences in a matrix if
    they are of the same length, otherwise in a list; (ii) \code{TRUE}:
    returns the sequences in a matrix, or stops with an error if they
    are of different lengths; (iii) \code{FALSE}: always returns the
    sequences in a list.}
  \item{type}{a character string giving the type of the sequences: one of
    \code{"DNA"} or \code{"AA"} (case-independent, can be abbreviated).}
  \item{offset}{the value to be added to the quality scores (the default
    applies to the Sanger format and should work for most recent FASTQ
    files).}
}
\details{
  \code{read.dna} follows the interleaved and sequential formats defined
  in PHYLIP (Felsenstein, 1993) but with the original feature than there
  is no restriction on the lengths of the taxa names. For these two
  formats, the first line of the file must contain the dimensions of the
  data (the numbers of taxa and the numbers of nucleotides); the
  sequences are considered as aligned and thus must be of the same
  lengths for all taxa. For the FASTA and FASTQ formats, the conventions
  defined in the references are followed; the sequences are taken as
  non-aligned. For all formats, the nucleotides can be arranged in any
  way with blanks and line-breaks inside (with the restriction that the
  first ten nucleotides must be contiguous for the interleaved and
  sequential formats, see below). The names of the sequences are read in
  the file. Particularities for each format are detailed below.

\itemize{
  \item{Interleaved:}{the function starts to read the sequences after it
    finds one or more spaces (or tabulations). All characters before the
    sequences are taken as the taxa names after removing the leading and
    trailing spaces (so spaces in taxa names are not allowed). It is
    assumed that the taxa names are not repeated in the subsequent
    blocks of nucleotides.}

  \item{Sequential:}{the same criterion than for the interleaved format
    is used to start reading the sequences and the taxa names; the
    sequences are then read until the number of nucleotides specified in
    the first line of the file is reached. This is repeated for each taxa.}

  \item{Clustal:}{this is the format output by the Clustal programs
    (.aln). It is close to the interleaved format: the differences are
    that the dimensions of the data are not indicated in the file, and
    the names of the sequences are repeated in each block.}

  \item{FASTA:}{This looks like the sequential format but the taxa names
    (or a description of the sequence) are on separate lines beginning
    with a `greater than' character `>' (there may be leading spaces
    before this character). These lines are taken as taxa names after
    removing the `>' and the possible leading and trailing spaces. All
    the data in the file before the first sequence are ignored.}
}

The FASTQ format is explained in the references.

Compressed files must be read through connections (see examples).
\code{read.fastq} can read compressed files directly (see
examples).
}
\value{
  a matrix or a list (if \code{format = "fasta"}) of DNA sequences
  stored in binary format, or of mode character (if \code{as.character =
    "TRUE"}).

  \code{read.FASTA} always returns a list of class \code{"DNAbin"} or
  \code{"AAbin"}.

  \code{read.fastq} returns a list of class \code{"DNAbin"} with an
  atrribute \code{"QUAL"} (see examples).
}
\references{
  Anonymous. FASTA format. \url{https://en.wikipedia.org/wiki/FASTA_format}

  Anonymous. FASTQ format. \url{https://en.wikipedia.org/wiki/FASTQ_format}

  Felsenstein, J. (1993) Phylip (Phylogeny Inference Package) version
  3.5c. Department of Genetics, University of Washington.
  \url{http://evolution.genetics.washington.edu/phylip/phylip.html}
}
\seealso{
  \code{\link{read.GenBank}}, \code{\link{write.dna}},
  \code{\link{DNAbin}}, \code{\link{dist.dna}}, \code{\link{woodmouse}}
}
\author{Emmanuel Paradis and RJ Ewing}
\examples{
## a small extract from data(woddmouse) in sequential format:
cat("3 40",
"No305     NTTCGAAAAACACACCCACTACTAAAANTTATCAGTCACT",
"No304     ATTCGAAAAACACACCCACTACTAAAAATTATCAACCACT",
"No306     ATTCGAAAAACACACCCACTACTAAAAATTATCAATCACT",
file = "exdna.txt", sep = "\n")
ex.dna <- read.dna("exdna.txt", format = "sequential")
str(ex.dna)
ex.dna
## the same data in interleaved format...
cat("3 40",
"No305     NTTCGAAAAA CACACCCACT",
"No304     ATTCGAAAAA CACACCCACT",
"No306     ATTCGAAAAA CACACCCACT",
"          ACTAAAANTT ATCAGTCACT",
"          ACTAAAAATT ATCAACCACT",
"          ACTAAAAATT ATCAATCACT",
file = "exdna.txt", sep = "\n")
ex.dna2 <- read.dna("exdna.txt")
## ... in clustal format...
cat("CLUSTAL (ape) multiple sequence alignment", "",
"No305     NTTCGAAAAACACACCCACTACTAAAANTTATCAGTCACT",
"No304     ATTCGAAAAACACACCCACTACTAAAAATTATCAACCACT",
"No306     ATTCGAAAAACACACCCACTACTAAAAATTATCAATCACT",
"           ************************** ******  ****",
file = "exdna.txt", sep = "\n")
ex.dna3 <- read.dna("exdna.txt", format = "clustal")
## ... and in FASTA format
cat(">No305",
"NTTCGAAAAACACACCCACTACTAAAANTTATCAGTCACT",
">No304",
"ATTCGAAAAACACACCCACTACTAAAAATTATCAACCACT",
">No306",
"ATTCGAAAAACACACCCACTACTAAAAATTATCAATCACT",
file = "exdna.fas", sep = "\n")
ex.dna4 <- read.dna("exdna.fas", format = "fasta")
## They are the same:
identical(ex.dna, ex.dna2)
identical(ex.dna, ex.dna3)
identical(ex.dna, ex.dna4)

## How to read compressed files:

## create the ZIP file:
zip("exdna.fas.zip", "exdna.fas")
## create the GZ file with a connection:
con <- gzfile("exdna.fas.gz", "wt")
cat(">No305", "NTTCGAAAAACACACCCACTACTAAAANTTATCAGTCACT",
    ">No304", "ATTCGAAAAACACACCCACTACTAAAAATTATCAACCACT",
    ">No306", "ATTCGAAAAACACACCCACTACTAAAAATTATCAATCACT",
    file = con, sep = "\n")
close(con)

ex.dna5 <- read.dna(unz("exdna.fas.zip", "exdna.fas"), "fasta")
ex.dna6 <- read.dna(gzfile("exdna.fas.gz"), "fasta")
identical(ex.dna5, ex.dna4)
identical(ex.dna6, ex.dna4)

unlink("exdna.txt")
unlink("exdna.fas")
unlink("exdna.fas.zip")
unlink("exdna.fas.gz")

## read a FASTQ file from 1000 Genomes:
\dontrun{
a <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/"
b <- "SRR062641.filt.fastq.gz"
URL <- paste0(a, b)
download.file(URL, b)
X <- read.fastq(b)
X # 109,811 sequences
## get the qualities of the first sequence:
(qual1 <- attr(X, "QUAL")[[1]])
## the corresponding probabilities:
10^(-qual1/10)
## get the mean quality for each sequence:
mean.qual <- sapply(attr(X, "Q"), mean)
## can do the same for var, sd, ...
}}
\keyword{IO}
