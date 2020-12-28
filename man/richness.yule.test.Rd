\name{richness.yule.test}
\alias{richness.yule.test}
\title{Test of Diversification-Shift With the Yule Process}
\description{
  This function performs a test of shift in diversification rate using
  probabilities from the Yule process.
}
\usage{
richness.yule.test(x, t)
}
\arguments{
  \item{x}{a matrix or a data frame with at least two columns: the first
    one gives the number of species in clades with a trait supposed to
    increase or decrease diversification rate, and the second one the number of
    species in the sister-clades without the trait. Each
    row represents a pair of sister-clades.}
  \item{t}{a numeric vector giving the divergence times of each pair of
    clades in \code{x}.}
}
\value{
  a data frame with the \eqn{\chi^2}{chi2}, the number of degrees of
  freedom (= 1), and the \emph{P}-value.
}
\references{
  Paradis, E. (2012) Shift in diversification in sister-clade
  comparisons: a more powerful test. \emph{Evolution}, \bold{66},
  288--295.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{slowinskiguyer.test}}, \code{\link{mcconwaysims.test}},
  \code{\link{diversity.contrast.test}}
}
\examples{
### see example(mcconwaysims.test)
}
\keyword{htest}
