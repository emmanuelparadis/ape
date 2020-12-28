\name{read.GenBank}
\alias{read.GenBank}
\title{Read DNA Sequences from GenBank via Internet}
\usage{
read.GenBank(access.nb, seq.names = access.nb, species.names = TRUE,
             as.character = FALSE, chunk.size = 400, quiet = TRUE)
}
\description{
  This function connects to the GenBank database, and reads nucleotide
  sequences using accession numbers given as arguments.
}
\arguments{
  \item{access.nb}{a vector of mode character giving the accession numbers.}
  \item{seq.names}{the names to give to each sequence; by default the
    accession numbers are used.}
  \item{species.names}{a logical indicating whether to attribute the
    species names to the returned object.}
  \item{as.character}{a logical controlling whether to return the
    sequences as an object of class \code{"DNAbin"} (the default).}
  \item{chunk.size}{the number of sequences downloaded together (see
    details).}
  \item{quiet}{a logical value indicating whether to show the progress
    of the downloads. If \code{TRUE}, will also print the (full) name of
    the FASTA file containing the downloaded sequences.}
}
\details{
  The function uses the site \url{https://www.ncbi.nlm.nih.gov/} from
  where the sequences are retrieved.

  If \code{species.names = TRUE}, the returned list has an attribute
  \code{"species"} containing the names of the species taken from the
  field ``ORGANISM'' in GenBank.

  Since \pkg{ape} 3.6, this function retrieves the sequences in FASTA
  format: this is more efficient and more flexible (scaffolds and
  contigs can be read) than what was done in previous versions. The
  option \code{gene.names} has been removed in \pkg{ape} 5.4; this
  information is also present in the description.

  Setting \code{species.names = FALSE} is much faster (could be useful
  if you read a series of scaffolds or contigs, or if you already have
  the species names).

  The argument \code{chunk.size} is set by default to 400 which is
  likely to work in many cases. If an error occurs such as ``Cannot open
  file \dots'' showing the list of the accession numbers, then you may
  try decreasing \code{chunk.size} to 200 or 300.

  If \code{quiet = FALSE}, the display is done chunk by chunk, so the
  message ``Downloading sequences: 400 / 400 ...'' means that the
  download from sequence 1 to sequence 400 is under progress (it is not
  possible to display a more accurate message because the download
  method depends on the platform).
}
\value{
  A list of DNA sequences made of vectors of class \code{"DNAbin"}, or
  of single characters (if \code{as.character = TRUE}) with two
  attributes (species and description).
}
\seealso{
  \code{\link{read.dna}}, \code{\link{write.dna}},
  \code{\link{dist.dna}}, \code{\link{DNAbin}}
}
\author{Emmanuel Paradis}
\examples{
## This won't work if your computer is not connected
## to the Internet

## Get the 8 sequences of tanagers (Ramphocelus)
## as used in Paradis (1997)
ref <- c("U15717", "U15718", "U15719", "U15720",
         "U15721", "U15722", "U15723", "U15724")
## Copy/paste or type the following commands if you
## want to try them.
\dontrun{
Rampho <- read.GenBank(ref)
## get the species names:
attr(Rampho, "species")
## build a matrix with the species names and the accession numbers:
cbind(attr(Rampho, "species"), names(Rampho))
## print the first sequence
## (can be done with `Rampho$U15717' as well)
Rampho[[1]]
## the description from each FASTA sequence:
attr(Rampho, "description")
}
}
\keyword{IO}
