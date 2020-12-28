\name{mcconwaysims.test}
\alias{mcconwaysims.test}
\title{McConway-Sims Test of Homogeneous Diversification}
\description{
  This function performs the McConway--Sims test that a trait or
  variable does not affect diversification rate.
}
\usage{
mcconwaysims.test(x)
}
\arguments{
  \item{x}{a matrix or a data frame with at least two columns: the first
    one gives the number of species in clades with a trait supposed to
    increase or decrease diversification rate, and the second one the number of
    species in the sister-clades without the trait. Each
    row represents a pair of sister-clades.}
}
\details{
  The McConway--Sims test compares a series of sister-clades where one
  of the two is characterized by a trait supposed to affect
  diversification rate. The null hypothesis is that the trait does not
  affect diversification. The alternative hypothesis is that
  diversification rate is increased or decreased by the trait (by
  contrast to the Slowinski--Guyer test). The test is a likelihood-ratio
  of a null Yule model and an alternative model with two parameters.
}
\value{
  a data frame with the \eqn{\chi^2}{chi2}, the number of degrees of
  freedom, and the \emph{P}-value.
}
\references{
  McConway, K. J. and Sims, H. J. (2004) A likelihood-based method for
  testing for nonstochastic variation of diversification rates in
  phylogenies. \emph{Evolution}, \bold{58}, 12--23.

  Paradis, E. (2012) Shift in diversification in sister-clade
  comparisons: a more powerful test. \emph{Evolution}, \bold{66},
  288--295.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{balance}}, \code{\link{slowinskiguyer.test}},
  \code{\link[geiger]{rc}} in \pkg{geiger},
  \code{\link[apTreeshape]{shift.test}} in \pkg{apTreeshape}
}
\examples{
### simulate 10 clades with lambda = 0.1 and mu = 0.09:
n0 <- replicate(10, balance(rbdtree(.1, .09, Tmax = 35))[1])
### simulate 10 clades with lambda = 0.15 and mu = 0.1:
n1 <- replicate(10, balance(rbdtree(.15, .1, Tmax = 35))[1])
x <- cbind(n1, n0)
mcconwaysims.test(x)
slowinskiguyer.test(x)
richness.yule.test(x, 35)
}
\keyword{htest}
