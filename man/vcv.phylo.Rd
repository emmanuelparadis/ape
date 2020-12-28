\name{vcv}
\alias{vcv}
\alias{vcv.phylo}
\alias{vcv.corPhyl}
\title{Phylogenetic Variance-covariance or Correlation Matrix}
\usage{
vcv(phy, ...)
\method{vcv}{phylo}(phy, model = "Brownian", corr = FALSE, ...)
\method{vcv}{corPhyl}(phy, corr = FALSE, ...)
}
\arguments{
  \item{phy}{an object of the correct class (see above).}
  \item{model}{a character giving the model used to compute the
    variances and covariances; only \code{"Brownian"} is available (for
    other models, a correlation structure may be used).}
  \item{corr}{a logical indicating whether the correlation matrix should
    be returned (\code{TRUE}); by default the variance-covariance matrix
    is returned (\code{FALSE}).}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\description{
  This function computes the expected variances and covariances of a
  continuous trait assuming it evolves under a given model.

  This is a generic function with methods for objects of class
  \code{"phylo"} and \code{"corPhyl"}.
}
\value{
  a numeric matrix with the names of the tips as colnames and rownames.
}
\references{
  Garland, T. Jr. and Ives, A. R. (2000) Using the past to predict the
  present: confidence intervals for regression equations in phylogenetic
  comparative methods. \emph{American Naturalist}, \bold{155}, 346--364.
}

\author{Emmanuel Paradis}
\note{
  Do not confuse this function with \code{\link[stats]{vcov}} which
  computes the variance-covariance matrix among parameters of a fitted
  model object.
}
\seealso{
  \code{\link{corBrownian}}, \code{\link{corMartins}},
  \code{\link{corGrafen}}, \code{\link{corPagel}},
  \code{\link{corBlomberg}}, \code{\link{vcv2phylo}}
}
\examples{
tr <- rtree(5)
## all are the same:
vcv(tr)
vcv(corBrownian(1, tr))
vcv(corPagel(1, tr))
}
\keyword{manip}
\keyword{multivariate}
