\name{read.gff}
\alias{read.gff}
\title{Read GFF Files}
\description{
  This function reads a file in general feature format version 3 (GFF3)
  and returns a data frame.
}
\usage{
read.gff(file, na.strings = c(".", "?"), GFF3 = TRUE)
}
\arguments{
  \item{file}{a file name specified by a character string.}
  \item{na.strings}{the strings in the GFF file that will be converted
    as NA's (missing values).}
  \item{GFF3}{a logical value specifying whether if the file is
    formatted according to version 3 of GFF.}
}
\details{
  The returned data frame has its (column) names correctly set (see
  References) and the categorical variables (seqid, source, type,
  strand, and phase) set as factors.

  This function should be more efficient than using \code{read.delim}.

  GFF2 (aka GTF) files can also be read: use \code{GFF3 = FALSE} to have
  the correct field names. Note that GFF2 files and GFF3 files have the
  same structure, although some fields are slightly different (see
  reference).

  The file can be gz-compressed (see examples), but not zipped.
}
\value{NULL}
\author{Emmanuel Paradis}
\references{
  \url{https://en.wikipedia.org/wiki/General_feature_format}
}
\examples{
\dontrun{
## requires to be connected on Internet
d <- "ftp://ftp.ensembl.org/pub/release-86/gff3/homo_sapiens/"
f <- "Homo_sapiens.GRCh38.86.chromosome.MT.gff3.gz"
download.file(paste0(d, f), "mt_gff3.gz")
gff.mito <- read.gff("mt_gff3.gz")
## the lengths of the sequence features:
gff.mito$end - (gff.mito$start - 1)
table(gff.mito$type)
## where the exons start:
gff.mito$start[gff.mito$type == "exon"]
}
}
\keyword{IO}
