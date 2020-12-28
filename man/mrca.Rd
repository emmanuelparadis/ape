\name{mrca}
\alias{mrca}
\alias{getMRCA}
\title{Find Most Recent Common Ancestors Between Pairs}
\description{
  \code{mrca} returns for each pair of tips (and nodes) its most
  recent common ancestor (MRCA).

  \code{getMRCA} returns the MRCA of two or more tips.
}
\usage{
mrca(phy, full = FALSE)
getMRCA(phy, tip)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{full}{a logical indicating whether to return the MRCAs among
    all tips and nodes (if \code{TRUE}); the default is to return only
    the MRCAs among tips.}
  \item{tip}{a vector of mode numeric or character specifying the tips;
    can also be node numbers.}
}
\details{
  For \code{mrca}, the diagonal is set to the number of the tips (and
  nodes if \code{full = TRUE}). If \code{full = FALSE}, the colnames and
  rownames are set with the tip labels of the tree; otherwise the
  numbers are given as names.

  For \code{getMRCA}, if \code{tip} is of length one or zero then
  \code{NULL} is returned.
}
\value{
  a matrix of mode numeric (\code{mrca}) or a single numeric value
  (\code{getMRCA}).
}
\author{Emmanuel Paradis, Klaus Schliep, Joseph W. Brown}
\keyword{manip}
