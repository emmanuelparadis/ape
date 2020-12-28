\name{write.nexus}
\alias{write.nexus}
\title{Write Tree File in Nexus Format}
\usage{
write.nexus(..., file = "", translate = TRUE)
}
\arguments{
  \item{\dots}{either (i) a single object of class \code{"phylo"}, (ii) a
    series of such objects separated by commas, or (iii) a list
    containing such objects.}
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string; if \code{file = ""} (the default) then the
    tree is written on the standard output connection.}
  \item{translate}{a logical, if \code{TRUE} (the default) a translation
    of the tip labels is done which are replaced in the parenthetic
    representation with tokens.}
}
\description{
  This function writes trees in a file with the NEXUS format.
}
\details{
  If several trees are given, they must all have the same tip labels.

  If among the objects given some are not trees of class \code{"phylo"},
  they are simply skipped and not written in the file.

  See \code{\link{write.tree}} for details on how tip (and node) labels
  are checked before being printed.
}
\value{
  None (invisible `NULL').
}
\references{
  Maddison, D. R., Swofford, D. L. and Maddison, W. P. (1997) NEXUS: an
  extensible file format for systematic information. \emph{Systematic
    Biology}, \bold{46}, 590--621.
}

\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.nexus}}, \code{\link{read.tree}},
  \code{\link{write.tree}}, \code{\link{read.nexus.data}},
  \code{\link{write.nexus.data}}
}
\keyword{manip}
\keyword{IO}
