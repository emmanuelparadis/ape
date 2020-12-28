\name{multiphylo}
\alias{multiphylo}
\alias{[.multiPhylo}
\alias{[[.multiPhylo}
\alias{$.multiPhylo}
\alias{[<-.multiPhylo}
\alias{[[<-.multiPhylo}
\alias{$<-.multiPhylo}
\title{Manipulating Lists of Trees}
\description{
  These are extraction and replacement operators for lists of trees
  stored in the class \code{"multiPhylo"}.
}
\usage{
\method{[}{multiPhylo}(x, i)
\method{[[}{multiPhylo}(x, i)
\method{$}{multiPhylo}(x, name)
\method{[}{multiPhylo}(x, ...) <- value
\method{[[}{multiPhylo}(x, ...) <- value
\method{$}{multiPhylo}(x, ...) <- value
}
\arguments{
  \item{x, value}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
  \item{i}{index(ices) of the tree(s) to select from a list; this may be a
    vector of integers, logicals, or names.}
  \item{name}{a character string specifying the tree to be extracted.}
  \item{\dots}{index(ices) of the tree(s) to replace; this may be a
    vector of integers, logicals, or names.}
}
\details{
  The subsetting operator \code{[} keeps the class correctly
  (\code{"multiPhylo"}).

The replacement operators check the labels of \code{value} if \code{x}
has a single vector of tip labels for all trees (see examples).
}
\value{
  An object of class \code{"phylo"} (\code{[[}, \code{$}) or of class
  \code{"multiPhylo"} (\code{[} and the replacement operators).
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{summary.phylo}}, \code{\link{c.phylo}}
}
\examples{
x <- rmtree(10, 20)
names(x) <- paste("tree", 1:10, sep = "")
x[1:5]
x[1] # subsetting
x[[1]] # extraction
x$tree1 # same than above
x[[1]] <- rtree(20)

y <- .compressTipLabel(x)
## up to here 'x' and 'y' have exactly the same information
## but 'y' has a unique vector of tip labels for all the trees
x[[1]] <- rtree(10) # no error
try(y[[1]] <- rtree(10)) # error

try(x[1] <- rtree(20)) # error
## use instead one of the two:
x[1] <- list(rtree(20))
x[1] <- c(rtree(20))

x[1:5] <- rmtree(5, 20) # replacement
x[11:20] <- rmtree(10, 20) # elongation
x # 20 trees
}
\keyword{manip}
