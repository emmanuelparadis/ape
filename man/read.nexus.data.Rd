\name{read.nexus.data}
\alias{read.nexus.data}
\alias{nexus2DNAbin}
\title{
  Read Character Data In NEXUS Format
}
\description{
  \code{read.nexus.data} reads a file with sequences in the NEXUS
  format. \code{nexus2DNAbin} is a helper function to convert the output
  from the previous function into the class \code{"DNAbin"}.
}
\usage{
read.nexus.data(file)
nexus2DNAbin(x)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode
    character, or a double-quoted string.}
  \item{x}{an object output by \code{read.nexus.data}.}
}
\details{
  This parser tries to read data from a file written in a
  \emph{restricted} NEXUS format (see examples below).

  Please see files \file{data.nex} and \file{taxacharacters.nex} for
  examples of formats that will work.

  Some noticeable exceptions from the NEXUS standard (non-exhaustive
  list):

  \itemize{
    \item{\bold{I}}{Comments must be either on separate lines or at the
      end of lines. Examples:\cr
      \code{[Comment]} \bold{--- OK}\cr
      \code{Taxon ACGTACG [Comment]} \bold{--- OK}\cr
      \code{[Comment line 1}

      \code{Comment line 2]} \bold{--- NOT OK!}\cr
      \code{Tax[Comment]on ACG[Comment]T} \bold{--- NOT OK!}}
    \item{\bold{II}}{No spaces (or comments) are allowed in the
      sequences. Examples:\cr
      \code{name ACGT} \bold{--- OK}\cr
      \code{name AC GT} \bold{--- NOT OK!}}
    \item{\bold{III}}{No spaces are allowed in taxon names, not even if
      names are in single quotes. That is, single-quoted names are not
      treated as such by the parser. Examples:\cr
      \code{Genus_species} \bold{--- OK}\cr
      \code{'Genus_species'} \bold{--- OK}\cr
      \code{'Genus species'} \bold{--- NOT OK!}}
    \item{\bold{IV}}{The trailing \code{end} that closes the
      \code{matrix} must be on a separate line. Examples:\cr
      \code{taxon AACCGGT}

      \code{end;} \bold{--- OK}\cr
      \code{taxon AACCGGT;}

      \code{end;} \bold{--- OK}\cr
      \code{taxon AACCCGT; end;} \bold{--- NOT OK!}}
    \item{\bold{V}}{Multistate characters are not allowed. That is,
      NEXUS allows you to specify multiple character states at a
      character position either as an uncertainty, \code{(XY)}, or as an
      actual appearance of multiple states, \code{\{XY\}}. This is
      information is not handled by the parser. Examples:\cr
      \code{taxon 0011?110} \bold{--- OK}\cr
      \code{taxon 0011{01}110} \bold{--- NOT OK!}\cr
      \code{taxon 0011(01)110} \bold{--- NOT OK!}}
    \item{\bold{VI}}{The number of taxa must be on the same line as
      \code{ntax}. The same applies to \code{nchar}. Examples:\cr
      \code{ntax = 12} \bold{--- OK}\cr
      \code{ntax =}

      \code{12} \bold{--- NOT OK!}}
    \item{\bold{VII}}{The word \dQuote{matrix} can not occur anywhere in
      the file before the actual \code{matrix} command, unless it is in
      a comment. Examples:\cr
      \code{BEGIN CHARACTERS;}

      \code{TITLE 'Data in file "03a-cytochromeB.nex"';}

      \code{DIMENSIONS  NCHAR=382;}

      \code{FORMAT DATATYPE=Protein GAP=- MISSING=?;}

      \code{["This is The Matrix"]} \bold{--- OK}

      \code{MATRIX}\cr

      \code{BEGIN CHARACTERS;}

      \code{TITLE 'Matrix in file "03a-cytochromeB.nex"';} \bold{--- NOT OK!}

      \code{DIMENSIONS  NCHAR=382;}

      \code{FORMAT DATATYPE=Protein GAP=- MISSING=?;}

      \code{MATRIX}}
  }
}
\value{
  A list of sequences each made of a single vector of mode character
  where each element is a (phylogenetic) character state.
}
\references{
  Maddison, D. R., Swofford, D. L. and Maddison, W. P. (1997) NEXUS: an
  extensible file format for systematic information. \emph{Systematic
    Biology}, \bold{46}, 590--621.
}
\author{Johan Nylander, Thomas Guillerme, and Klaus Schliep}
\seealso{
  \code{\link{read.nexus}}, \code{\link{write.nexus}},
  \code{\link{write.nexus.data}}
}
\examples{
## Use read.nexus.data to read a file in NEXUS format into object x
\dontrun{x <- read.nexus.data("file.nex")}
}
\keyword{file}
