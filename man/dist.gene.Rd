\name{dist.gene}
\alias{dist.gene}
\title{Pairwise Distances from Genetic Data}
\usage{
dist.gene(x, method = "pairwise", pairwise.deletion = FALSE,
          variance = FALSE)
}
\arguments{
  \item{x}{a matrix or a data frame (will be coerced as a matrix).}
  \item{method}{a character string specifying the method used to compute
    the distances; two choices are available: \code{"pairwise"} and
    \code{"percentage"}, or any unambiguous abbreviation of these.}
  \item{pairwise.deletion}{a logical indicating whether to delete the
    columns with missing data on a pairwise basis. The default is to
    delete the columns with at least one missing observation.}
  \item{variance}{a logical, indicates whether the variance of the
    distances should be returned (default to \code{FALSE}).}
}
\description{
  This function computes a matrix of distances between pairs of
  individuals from a matrix or a data frame of genetic data.
}
\details{
  This function is meant to be very general and accepts different kinds
  of data (alleles, haplotypes, SNP, DNA sequences, \dots). The rows of
  the data matrix represent the individuals, and the columns the loci.

  In the case of the pairwise method, the distance \eqn{d} between two
  individuals is the number of loci for which they differ, and the
  associated variance is \eqn{d(L - d)/L}, where \eqn{L} is the number
  of loci.

  In the case of the percentage method, this distance is divided by \eqn{L},
  and the associated variance is \eqn{d(1 - d)/L}.

  For more elaborate distances with DNA sequences, see the function
  \code{dist.dna}.
}
\note{
  Missing data (\code{NA}) are coded and treated in R's usual way.
}
\value{
  an object of class \code{dist}. If \code{variance = TRUE} an
  attribute called \code{"variance"} is given to the returned object.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{dist.dna}}, \code{\link{cophenetic.phylo}},
  \code{\link[stats]{dist}}
}
\keyword{manip}

