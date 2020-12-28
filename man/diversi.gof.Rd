\encoding{utf8}
\name{diversi.gof}
\alias{diversi.gof}
\title{Tests of Constant Diversification Rates}
\usage{
diversi.gof(x, null = "exponential", z = NULL)
}
\arguments{
  \item{x}{a numeric vector with the branching times.}
  \item{null}{a character string specifying the null distribution for
    the branching times. Only two choices are possible: either
    \code{"exponential"}, or \code{"user"}.}
  \item{z}{used if \code{null = "user"}; gives the expected distribution
    under the model.}
}
\description{
  This function computes two tests of the distribution of branching
  times using the \enc{Cramér}{Cramer}--von Mises and Anderson--Darling
  goodness-of-fit tests. By default, it is assumed that the
  diversification rate is constant, and an exponential distribution is
  assumed for the branching times. In this case, the expected
  distribution under this model is computed with a rate estimated from
  the data. Alternatively, the user may specify an expected cumulative
  density function (\code{z}): in this case, \code{x} and \code{z} must
  be of the same length. See the examples for how to compute the latter
  from a sample of expected branching times.
}
\details{
  The \enc{Cramér}{Cramer}--von Mises and Anderson--Darling tests
  compare the empirical density function (EDF) of the observations to an
  expected cumulative density function. By contrast to the
  Kolmogorov--Smirnov test where the greatest difference between these
  two functions is used, in both tests all differences are taken into
  account.

  The distributions of both test statistics depend on the null
  hypothesis, and on whether or not some parameters were estimated from
  the data. However, these distributions are not known precisely and
  critical values were determined by Stephens (1974) using
  simulations. These critical values were used for the present function.
}
\value{
  A NULL value is returned, the results are simply printed.
}
\references{
  Paradis, E. (1998) Testing for constant diversification rates using
  molecular phylogenies: a general approach based on statistical tests
  for goodness of fit. \emph{Molecular Biology and Evolution},
  \bold{15}, 476--479.

  Stephens, M. A. (1974) EDF statistics for goodness of fit and some
  comparisons. \emph{Journal of the American Statistical Association},
  \bold{69}, 730--737.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{branching.times}}, \code{\link{diversi.time}}
  \code{\link{ltt.plot}}, \code{\link{birthdeath}}, \code{\link{yule}},
  \code{\link{yule.cov}}
}
\examples{
data(bird.families)
x <- branching.times(bird.families)
### suppose we have a sample of expected branching times `y';
### for simplicity, take them from a uniform distribution:
y <- runif(500, 0, max(x) + 1) # + 1 to avoid A2 = Inf
### now compute the expected cumulative distribution:
x <- sort(x)
N <- length(x)
ecdf <- numeric(N)
for (i in 1:N) ecdf[i] <- sum(y <= x[i])/500
### finally do the test:
diversi.gof(x, "user", z = ecdf)
}
\keyword{univar}
