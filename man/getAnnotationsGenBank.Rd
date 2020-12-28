\name{getAnnotationsGenBank}
\alias{getAnnotationsGenBank}
\title{Read Annotations from GenBank}
\description{
  This function connects to the GenBank database and reads sequence annotations
  using accession number(s) given as argument.
}
\usage{
getAnnotationsGenBank(access.nb, quiet = TRUE)
}
\arguments{
  \item{access.nb}{a vector of mode character giving the accession
    numbers.}
  \item{quiet}{a logical value indicating whether to show the progress
    of the downloads.}
}
\details{
  The sequence annotations (a.k.a. feature list) are returned in a data
  frame with five or six columns: start, end, type, product, others, and
  gene (the last being optional). This is the same information that can
  be downloaded from NCBI's Web interface by clicking on `Send to:',
  `File', and then selecting `Feature Table' under `Format'.

  A warning is given if some features are incomplete (this information
  is then dropped from the returned object).

  A warning is given if some accession numbers are not found on GenBank.
}
\value{
  On of the followings: (i) a data frame if \code{access.nb} contains a
  single accession number; (ii) a list of data frames if
  \code{access.nb} contains several accession numbers, the names are set
  with \code{access.nb} (if some accession numbers are not found on
  GenBank, the corresponding entries are set to \code{NULL}); (iii)
  \code{NULL} if all accession numbers are not found on GenBank.
}
\seealso{\code{\link{read.GenBank}}, \code{\link{read.gff}},
  \code{\link{DNAbin}}
}
\author{Emmanuel Paradis}
\references{\url{https://www.ncbi.nlm.nih.gov/Sequin/table.html}}
\examples{
## The 8 sequences of tanagers (Ramphocelus):
ref <- c("U15717", "U15718", "U15719", "U15720",
         "U15721", "U15722", "U15723", "U15724")
## Copy/paste or type the following commands if you
## want to try them.
\dontrun{
annot.rampho <- getAnnotationsGenBank(ref)
annot.rampho
## check all annotations are the same:
unique(do.call(rbind, annot.rampho)[, -5])
}
}
\keyword{IO}
