\name{compar.ou}
\alias{compar.ou}
\title{Ornstein--Uhlenbeck Model for Continuous Characters}
\usage{
compar.ou(x, phy, node = NULL, alpha = NULL)
}
\arguments{
  \item{x}{a numeric vector giving the values of a continuous
    character.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{node}{a vector giving the number(s) of the node(s) where the
    parameter `theta' (the trait optimum) is assumed to change. The
    node(s) can be specified with their labels if \code{phy} has node
    labels. By default there is no change (same optimum thoughout lineages).}
  \item{alpha}{the value of \eqn{\alpha}{alpha} to be used when fitting
    the model. By default, this parameter is estimated (see details).}
}
\description{
  This function fits an Ornstein--Uhlenbeck model giving a phylogenetic
  tree, and a continuous character. The user specifies the node(s) where
  the optimum changes. The parameters are estimated by maximum
  likelihood; their standard-errors are computed assuming normality of
  these estimates.
}
\details{
  The Ornstein--Uhlenbeck (OU) process can be seen as a generalization
  of the Brownian motion process. In the latter, characters are assumed
  to evolve randomly under a random walk, that is change is equally
  likely in any direction. In the OU model, change is more likely
  towards the direction of an optimum (denoted \eqn{\theta}{theta}) with
  a strength controlled by a parameter denoted \eqn{\alpha}{alpha}.

  The present function fits a model where the optimum parameter
  \eqn{\theta}{theta}, is allowed to vary throughout the tree. This is
  specified with the argument \code{node}: \eqn{\theta}{theta} changes
  after each node whose number is given there. Note that the optimum
  changes \emph{only} for the lineages which are descendants of this
  node.

  Hansen (1997) recommends to not estimate \eqn{\alpha}{alpha} together
  with the other parameters. The present function allows this by giving
  a numeric value to the argument \code{alpha}. By default, this
  parameter is estimated, but this seems to yield very large
  standard-errors, thus validating Hansen's recommendation. In practice,
  a ``poor man estimation'' of \eqn{\alpha}{alpha} can be done by
  repeating the function call with different values of \code{alpha}, and
  selecting the one that minimizes the deviance (see Hansen 1997 for an
  example).

  If \code{x} has names, its values are matched to the tip labels of
  \code{phy}, otherwise its values are taken to be in the same order
  than the tip labels of \code{phy}.

  The user must be careful here since the function requires that both
  series of names perfectly match, so this operation may fail if there
  is a typing or syntax error. If both series of names do not match, the
  values in the \code{x} are taken to be in the same order than the tip
  labels of \code{phy}, and a warning message is issued.
}
\note{
  The inversion of the variance-covariance matrix in the likelihood
  function appeared as somehow problematic. The present implementation
  uses a Cholevski decomposition with the function
  \code{\link[base]{chol2inv}} instead of the usual function
  \code{\link[base]{solve}}.
}
\value{
  an object of class \code{"compar.ou"} which is list with the following
  components:

  \item{deviance}{the deviance (= -2 * loglik).}
  \item{para}{a data frame with the maximum likelihood estimates and
    their standard-errors.}
  \item{call}{the function call.}
}
\references{
  Hansen, T. F. (1997) Stabilizing selection and the comparative
  analysis of adaptation. \emph{Evolution}, \bold{51}, 1341--1351.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{ace}}, \code{\link{compar.lynch}},
  \code{\link{corBrownian}}, \code{\link{corMartins}}, \code{\link{pic}}
}
\examples{
data(bird.orders)
### This is likely to give you estimates close to 0, 1, and 0
### for alpha, sigma^2, and theta, respectively:
compar.ou(x <- rnorm(23), bird.orders)
### Much better with a fixed alpha:
compar.ou(x, bird.orders, alpha = 0.1)
### Let us 'mimick' the effect of different optima
### for the two clades of birds...
x <- c(rnorm(5, 0), rnorm(18, 5))
### ... the model with two optima:
compar.ou(x, bird.orders, node = 25, alpha = .1)
### ... and the model with a single optimum:
compar.ou(x, bird.orders, node = NULL, alpha = .1)
### => Compare both models with the difference in deviances
##     which follows a chi^2 with df = 1.
}
\keyword{models}
