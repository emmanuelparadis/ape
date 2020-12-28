\name{write.nexus.data}
\alias{write.nexus.data}
\title{Write Character Data in NEXUS Format}
\description{
  This function writes in a file a list of data in the NEXUS format. The
  names of the vectors of the list are used as taxon names.
}
\usage{
write.nexus.data(x, file, format = "dna", datablock = TRUE,
                 interleaved = TRUE, charsperline = NULL,
                 gap = NULL, missing = NULL)
}
\arguments{
  \item{x}{a matrix or a list of data each made of a single vector
    of mode character where each element is a character state (e.g.,
    \dQuote{A}, \dQuote{C}, \dots) Objects of class of \dQuote{DNAbin}
    are accepted.}
  \item{file}{a file name specified by either a variable of mode
    character, or a double-quoted string.}
  \item{format}{a character string specifying the format of the
    sequences. Four choices are possible: \code{"dna"} (the default)
    \code{"protein"}, \code{"standard"} or \code{"continuous"} or any
    unambiguous abbreviation of these (case insensitive).}
  \item{datablock}{a logical, if \code{TRUE} the data are written in a
    single DATA block. If \code{FALSE}, the data are written in TAXA and
    CHARACTER blocks. Default is \code{TRUE}.}
  \item{interleaved}{a logical, if \code{TRUE} the data is written in
    interleaved format with number of characters per line as specified
    with \code{charsperline = numerical_value}. If \code{FALSE}, the
    data are written in sequential format. Default is \code{TRUE}.}
  \item{charsperline}{a numeric value specifying the number of
    characters per line when used with \code{interleaved =
      TRUE}. Default is 80.}
  \item{gap}{a character specifying the symbol for gap. Default is
    \dQuote{\code{-}}.}
  \item{missing}{a character specifying the symbol for missing
    data. Default is \dQuote{\code{?}}.}
}
\details{
  If the sequences have no names, then they are given \dQuote{1},
  \dQuote{2}, ..., as names in the file.

  Sequences must be all of the same length.
}
\value{
  None (invisible \sQuote{NULL}).
}
\references{
  Maddison, D. R., Swofford, D. L. and Maddison, W. P. (1997) NEXUS: an
  extensible file format for systematic information. \emph{Systematic
    Biology}, \bold{46}, 590--621.
}
\author{Johan Nylander \email{nylander@scs.fsu.edu} and Thomas Guillerme}
\seealso{
  \code{\link{read.nexus}},\code{\link{write.nexus}},
  \code{\link{read.nexus.data}}
}
\examples{
\dontrun{
## Write interleaved DNA data with 100 characters per line in a DATA block
data(woodmouse)
write.nexus.data(woodmouse, file= "wood.ex.nex", interleaved = TRUE, charsperline = 100)
## Write sequential DNA data in TAXA and CHARACTERS blocks
data(cynipids)
write.nexus.data(cynipids, file = "cyn.ex.nex", format = "protein",
                 datablock = FALSE, interleaved = FALSE)
unlink(c("wood.ex.nex", "cyn.ex.nex"))
}}
\keyword{file}
