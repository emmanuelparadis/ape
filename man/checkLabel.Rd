\name{checkLabel}
\alias{checkLabel}
\title{Checking Labels}
\description{
  Checking and correcting character strings, particularly before writing
  a Newick tree.
}
\usage{checkLabel(x)}
\arguments{\item{x}{a vector of mode character.}}
\details{
  This function deletes the leading and trailing spaces (including
  tabulations, new lines, and left or right parentheses at the beginning
  or end of the strings), substitutes the spaces inside the strings by
  underscores, and substitutes commas, colons, semicolons, and
  parentheses inside the strings by dashes.
}
\value{a vector of mode character.}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{makeLabel}}, \code{\link{makeNodeLabel}},
  \code{\link{mixedFontLabel}}, \code{\link{stripLabel}},
  \code{\link{updateLabel}}
}
\examples{
checkLabel(" Homo sapiens\t(Primates; World)   ")
}
\keyword{manip}
