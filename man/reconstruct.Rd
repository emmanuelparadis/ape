\name{reconstruct}
\alias{reconstruct}
\title{Continuous Ancestral Character Estimation}
\description{
  This function estimates ancestral character states, and the associated
  uncertainty, for continuous characters. It mainly works as the ace
  function, from which it differs, first, in the fact that computations
  are not performed by numerical optimisation but through matrix
  calculus. Second, besides classical Brownian-based reconstruction
  methods, it reconstructs ancestral states under Arithmetic Brownian
  Motion (ABM, i.e. Brownian with linear trend) and Ornstein-Uhlenbeck
  process (OU, i.e. Brownian with an attractive optimum).
}
\usage{
reconstruct(x, phyInit, method = "ML", alpha = NULL, CI = TRUE)
}
\arguments{
  \item{x}{a numerical vector.}
  \item{phyInit}{an object of class \code{"phylo"}.}
  \item{method}{a character specifying the method used for
    estimation. Six choices are possible: \code{"ML"}, \code{"REML"},
    \code{"GLS"}, \code{"GLS_ABM"}, \code{"GLS_OU"} or
    \code{"GLS_OUS"}.}
  \item{alpha}{a numerical value which accounts for the attractive
    strength parameter of \code{"GLS_OU"} or \code{"GLS_OUS"} (used only
    in these cases). If alpha = NULL (the default), then it is estimated
    by maximum likelihood using \code{optim} which may lead to
    convergence issue.}
  \item{CI}{a logical specifying whether to return the 95\% confidence
    intervals of the ancestral state estimates.}
}
\details{
  For \code{"ML"}, \code{"REML"} and \code{"GLS"}, the default model is
  Brownian motion. This model can be fitted by maximum likelihood
  (\code{method = "ML"}, Felsenstein 1973, Schluter et al. 1997) - the
  default, residual maximum likelihood (\code{method = "REML"}), or
  generalized least squares (\code{method = "GLS"}, Martins and Hansen
  1997, Garland T and Ives AR 2000).

  \code{"GLS_ABM"} is based on Brownian motion with trend model. Both
  \code{"GLS_OU"} and \code{"GLS_OUS"} are based on Ornstein-Uhlenbeck
  model.

  \code{"GLS_OU"} and \code{"GLS_OUS"} differs in the fact that
  \code{"GLS_OUS"} assume that the process starts from the optimum,
  while the root state has to be estimated for \code{"GLS_OU"}, which
  may rise some issues (see Royer-Carenzi and Didier, 2016). Users may
  provide the attractive strength parameter \code{alpha}, for these two
  models.

  \code{"GLS_ABM"}, \code{"GLS_OU"} and \code{"GLS_OUS"} are all fitted
  by generalized least squares (Royer-Carenzi and Didier, 2016).
}
\value{
  an object of class \code{"ace"} with the following elements:

  \item{ace}{the estimates of the ancestral character values.}
  \item{CI95}{the estimated 95\% confidence intervals.}
  \item{sigma2}{if \code{method = "ML"}, the maximum likelihood estimate
    of the Brownian parameter.}
  \item{loglik}{if \code{method = "ML"}, the maximum log-likelihood.}
}
\references{
  Felsenstein, J. (1973) Maximum likelihood estimation of evolutionary
  trees from continuous characters. \emph{American Journal of Human
  Genetics}, \bold{25}, 471--492.

  Garland T. and Ives A.R. (2000) Using the past to predict the present:
  confidence intervals for regression equations in phylogenetic
  comparative methods. \emph{American Naturalist}, \bold{155},
  346--364.

  Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the
  comparative method: a general approach to incorporating phylogenetic
  information into the analysis of interspecific data. \emph{American
    Naturalist}, \bold{149}, 646--667.

  Royer-Carenzi, M. and Didier, G. (2016) A comparison of ancestral
  state reconstruction methods for quantitative
  characters. \emph{Journal of Theoretical Biology}, \bold{404},
  126--142.

  Schluter, D., Price, T., Mooers, A. O. and Ludwig, D. (1997)
  Likelihood of ancestor states in adaptive radiation. \emph{Evolution},
  \bold{51}, 1699--1711.

  Yang, Z. (2006) \emph{Computational Molecular Evolution}. Oxford:
  Oxford University Press.
}
\author{Manuela Royer-Carenzi, Gilles Didier}
\seealso{
  \code{\link{MPR}}, \code{\link{corBrownian}}, \code{\link{compar.ou}}

  Reconstruction of ancestral sequences can be done with the package
  \pkg{phangorn} (see function \code{?ancestral.pml}).
}
\note{
\code{GLS_ABM} should not be used on ultrametric tree.

\code{GLS_OU} may lead to aberrant reconstructions.
}
\examples{
### Some random data...
data(bird.orders)
x <- rnorm(23, m=100)
### Reconstruct ancestral quantitative characters:
reconstruct(x, bird.orders)
reconstruct(x, bird.orders, method = "GLS_OUS", alpha=NULL)
}
\keyword{models}
