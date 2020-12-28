\name{is.compatible}
\alias{is.compatible}
\alias{is.compatible.bitsplits}
\alias{arecompatible}
\title{Check Compatibility of Splits}
\description{
  \code{is.compatible} is a generic function with a method for the class
  \code{"bitsplits"}. It checks whether a set of splits is compatible
  using the \code{arecompatible} function.
}
\usage{
is.compatible(obj)
\method{is.compatible}{bitsplits}(obj)
arecompatible(x, y, n)
}
\arguments{
  \item{obj}{an object of class \code{"bitsplits"}.}
  \item{x, y}{a vector of mode raw\code{}.}
  \item{n}{the number of taxa in the splits.}
}
\value{
  \code{TRUE} if the splits are compatible, \code{FALSE} otherwise.
}
\author{Andrei Popescu}
\seealso{\code{\link{as.bitsplits}}}
\keyword{manip}
