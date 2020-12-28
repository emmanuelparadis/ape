\name{Initialize.corPhyl}
\alias{Initialize.corPhyl}
\title{Initialize a `corPhyl' Structure Object}
\usage{
	\method{Initialize}{corPhyl}(object, data, ...)
}
\arguments{
	\item{object}{An object inheriting from class \code{corPhyl}.}
	\item{data}{The data to use. If it contains rownames, they are matched with the tree tip labels, otherwise data are supposed to be in the same order than tip labels and a warning is sent.}
	\item{\dots}{some methods for this generic require additional arguments. None are used in this method.}
}
\description{
	Initialize a \code{corPhyl} correlation structure object.
	Does the same as \code{Initialize.corStruct}, but also checks the row names of data and builds an index.
}
\value{
	An initialized object of same class as \code{object}.
}
\author{Julien Dutheil \email{dutheil@evolbio.mpg.de}}
\seealso{
  \code{\link{corClasses}}, \code{\link[nlme]{Initialize.corStruct}}.
}
\keyword{models}
\keyword{manip}
