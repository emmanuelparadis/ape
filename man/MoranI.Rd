\name{Moran.I}
\alias{Moran.I}
\title{Moran's I Autocorrelation Index}
\usage{
  Moran.I(x, weight, scaled = FALSE, na.rm = FALSE,
          alternative = "two.sided")
}
\arguments{
  \item{x}{a numeric vector.}
  \item{weight}{a matrix of weights.}
  \item{scaled}{a logical indicating whether the coefficient should be
    scaled so that it varies between -1 and +1 (default to
    \code{FALSE}).}
  \item{na.rm}{a logical indicating whether missing values should be
    removed.}
  \item{alternative}{a character string specifying the alternative
    hypothesis that is tested against the null hypothesis of no
    phylogenetic correlation; must be of one "two.sided", "less", or
    "greater", or any unambiguous abbrevation of these.}
}
\description{
  This function computes Moran's I autocorrelation coefficient of
  \code{x} giving a matrix of weights using the method described by
  Gittleman and Kot (1990).
}
\details{
  The matrix \code{weight} is used as ``neighbourhood'' weights, and
  Moran's I coefficient is computed using the formula:
  \deqn{I = \frac{n}{S_0} \frac{\sum_{i=1}^n\sum_{j=1}^n w_{i,j}(y_i -
      \overline{y})(y_j - \overline{y})}{\sum_{i=1}^n {(y_i -
	\overline{y})}^2}}{\code{I = n/S0 * (sum\{i=1..n\} sum\{j=1..n\} wij(yi - ym))(yj - ym)
      / (sum\{i=1..n\} (yi - ym)^2)}}
  with
  \itemize{
    \item \eqn{y_i}{yi} = observations
    \item \eqn{w_{i,j}}{wij} = distance weight
    \item \eqn{n} = number of observations
    \item \eqn{S_0}{S0} = \eqn{\sum_{i=1}^n\sum_{j=1}^n wij}{\code{sum_{i=1..n} sum{j=1..n} wij}}
  }

  The null hypothesis of no phylogenetic correlation is tested assuming
  normality of I under this null hypothesis. If the observed value
  of I is significantly greater than the expected value, then the values
  of \code{x} are positively autocorrelated, whereas if Iobserved <
  Iexpected, this will indicate negative autocorrelation.
}
\value{
  A list containing the elements:

  \item{observed}{the computed Moran's I.}
  \item{expected}{the expected value of I under the null hypothesis.}
  \item{sd}{the standard deviation of I under the null hypothesis.}
  \item{p.value}{the P-value of the test of the null hypothesis against
    the alternative hypothesis specified in \code{alternative}.}
}
\references{
  Gittleman, J. L. and Kot, M. (1990) Adaptation: statistics and a null
  model for estimating phylogenetic effects. \emph{Systematic Zoology},
  \bold{39}, 227--241.
}
\author{Julien Dutheil \email{dutheil@evolbio.mpg.de} and
  Emmanuel Paradis}
\seealso{\code{\link{weight.taxo}}}
\examples{
tr <- rtree(30)
x <- rnorm(30)
## weights w[i,j] = 1/d[i,j]:
w <- 1/cophenetic(tr)
## set the diagonal w[i,i] = 0 (instead of Inf...):
diag(w) <- 0
Moran.I(x, w)
Moran.I(x, w, alt = "l")
Moran.I(x, w, alt = "g")
Moran.I(x, w, scaled = TRUE) # usualy the same
}
\keyword{models}
\keyword{regression}
