\name{corPagel}
\alias{corPagel}
\alias{coef.corPagel}
\alias{corMatrix.corPagel}
\title{Pagel's ``lambda'' Correlation Structure}
\usage{
corPagel(value, phy, form = ~1, fixed = FALSE)
\method{corMatrix}{corPagel}(object, covariate = getCovariate(object),
                   corr = TRUE, ...)
\method{coef}{corPagel}(object, unconstrained = TRUE, \dots)
}
\arguments{
  \item{value}{the (initial) value of the parameter
    \eqn{\lambda}{lambda}.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{form}{a one sided formula of the form ~ t, or ~ t | g, specifying the taxa covariate t and, optionally, a grouping factor g. A covariate for this correlation structure must be character valued, with entries matching the tip labels in the phylogenetic tree. When a grouping factor is present in form, the correlation structure is assumed to apply only to observations within the same grouping level; observations with different grouping levels are assumed to be uncorrelated. Defaults to ~ 1, which corresponds to using the order of the observations in the data as a covariate, and no groups.}
  \item{fixed}{a logical specifying whether \code{gls} should
    estimate \eqn{\lambda}{lambda} (the default) or keep it fixed.}
  \item{object}{an (initialized) object of class \code{"corPagel"}.}
  \item{covariate}{an optional covariate vector (matrix), or list of covariate vectors (matrices), at which values the correlation matrix, or list of correlation matrices, are to be evaluated. Defaults to getCovariate(object).}
  \item{corr}{a logical value specifying whether to return the
    correlation matrix (the default) or the variance-covariance matrix.}
  \item{unconstrained}{a logical value. If \code{TRUE} (the default),
    the coefficients are returned in unconstrained form (the same used
    in the optimization algorithm). If \code{FALSE} the coefficients are
    returned in ``natural'', possibly constrained, form.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  The correlation structure from the present model is derived from the
  Brownian motion model by multiplying the off-diagonal elements (i.e.,
  the covariances) by \eqn{\lambda}{lambda}. The variances are thus the
  same than for a Brownian motion model.
}
\value{
  an object of class \code{"corPagel"}, the coefficients from an object
  of this class, or the correlation matrix of an initialized object of
  this class. In most situations, only \code{corPagel} will be called
  by the user.
}
\author{Emmanuel Paradis}
\references{
  Freckleton, R. P., Harvey, P. H. and M. Pagel, M. (2002) Phylogenetic
  analysis and comparative data: a test and review of evidence.
  \emph{American Naturalist}, \bold{160}, 712--726.

  Pagel, M. (1999) Inferring the historical patterns of biological
  evolution. \emph{Nature}, \bold{401},877--884.
}
\keyword{models}
