\name{corBrownian}
\alias{corBrownian}
\alias{coef.corBrownian}
\alias{corMatrix.corBrownian}
\title{Brownian Correlation Structure}
\usage{
corBrownian(value=1, phy, form=~1)
\method{coef}{corBrownian}(object, unconstrained = TRUE, ...)
\method{corMatrix}{corBrownian}(object, covariate = getCovariate(object), corr = TRUE, ...)
}
\arguments{
  \item{value}{The \eqn{\gamma}{gamma} parameter (default to 1)}
  \item{phy}{An object of class \code{phylo} representing the phylogeny
    (with branch lengths) to consider}
  \item{object}{An (initialized) object of class \code{corBrownian}}
  \item{corr}{a logical value. If 'TRUE' the function returns the
    correlation matrix, otherwise it returns the variance/covariance matrix.}
  \item{form}{a one sided formula of the form ~ t, or ~ t | g, specifying the taxa covariate t and, optionally, a grouping factor g. A covariate for this correlation structure must be character valued, with entries matching the tip labels in the phylogenetic tree. When a grouping factor is present in form, the correlation structure is assumed to apply only to observations within the same grouping level; observations with different grouping levels are assumed to be uncorrelated. Defaults to ~ 1, which corresponds to using the order of the observations in the data as a covariate, and no groups.}
  \item{covariate}{an optional covariate vector (matrix), or list of covariate vectors (matrices), at which values the correlation matrix, or list of correlation matrices, are to be evaluated. Defaults to getCovariate(object).}
  \item{unconstrained}{a logical value. If 'TRUE' the coefficients are returned
    in unconstrained form (the same used in the optimization
    algorithm). If 'FALSE' the coefficients are returned in
    "natural", possibly constrained, form. Defaults to 'TRUE'}
  \item{\dots}{some methods for these generics require additional arguments.
    None are used in these methods.}
}
\description{
  Expected covariance under a Brownian model (Felsenstein 1985,	Martins
  and Hansen 1997)

  \deqn{V_{ij} = \gamma \times t_a}{Vij = gamma . ta}

  where \eqn{t_a}{ta} is the distance on the phylogeny between the root
  and the most recent common ancestor of taxa \eqn{i}{i} and \eqn{j}{j}
  and \eqn{\gamma}{gamma} is a constant.
}
\value{
  An object of class \code{corBrownian}, or the coefficient from an
  object of this class (actually sends \code{numeric(0)}), or the
  correlation matrix of an initialized object of this class.
}
\author{Julien Dutheil \email{dutheil@evolbio.mpg.de}}
\seealso{
  \code{\link{corClasses}}
}
\references{
  Felsenstein, J. (1985) Phylogenies and the comparative method.
  \emph{American Naturalist}, \bold{125}, 1--15.

  Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the comparative
  method: a general approach to incorporating phylogenetic information
  into the analysis of interspecific data. \emph{American Naturalist},
  \bold{149}, 646--667.
}
\keyword{models}
