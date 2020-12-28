\name{gammaStat}
\alias{gammaStat}
\title{Gamma-Statistic of Pybus and Harvey}
\usage{
gammaStat(phy)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
}
\description{
  This function computes the gamma-statistic which summarizes the
  information contained in the inter-node intervals of a phylogeny. It
  is assumed that the tree is ultrametric. Note that the function does
  not check that the tree is effectively ultrametric, so if it is not,
  the returned result may not be meaningful.
}
\value{
  a numeric vector of length one.
}
\details{
  The gamma-statistic is a summary of the information contained in the
  inter-node intervals of a phylogeny; it follows, under the assumption
  that the clade diversified with constant rates, a normal distribution
  with mean zero and standard-deviation unity (Pybus and Harvey
  2000). Thus, the null hypothesis that the clade diversified with
  constant rates may be tested with \code{2*(1 -
    pnorm(abs(gammaStat(phy))))} for a two-tailed test, or \code{1 -
    pnorm(abs(gammaStat(phy)))} for a one-tailed test, both returning
  the corresponding P-value.
}
\references{
  Pybus, O. G. and Harvey, P. H. (2000) Testing macro-evolutionary
  models using incomplete molecular phylogenies. \emph{Proceedings of
    the Royal Society of London. Series B. Biological Sciences},
  \bold{267}, 2267--2272.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{branching.times}}, \code{\link{ltt.plot}}, \code{\link{skyline}}
}
\keyword{univar}
