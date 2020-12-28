\name{birthdeath}
\alias{birthdeath}
\alias{print.birthdeath}
\title{Estimation of Speciation and Extinction Rates With Birth-Death Models}
\usage{
birthdeath(phy)
\method{print}{birthdeath}(x, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{x}{an object of class \code{"birthdeath"}.}
  \item{\dots}{further arguments passed to the \code{print} function.}
}
\description{
  This function fits by maximum likelihood a birth-death model to the
  branching times computed from a phylogenetic tree using the method of
  Nee et al. (1994).
}
\details{
  Nee et al. (1994) used a re-parametrization of the birth-death model
  studied by Kendall (1948) so that the likelihood has to be maximized
  over \emph{d/b} and \emph{b - d}, where \emph{b} is the birth rate,
  and \emph{d} the death rate. This is the approach used by the present
  function.

  This function computes the standard-errors of the estimated parameters
  using a normal approximations of the maximum likelihood estimates:
  this is likely to be inaccurate because of asymmetries of the
  likelihood function (Nee et al. 1995). In addition, 95 % confidence
  intervals of both parameters are computed using profile likelihood:
  they are particularly useful if the estimate of \emph{d/b} is at the
  boundary of the parameter space (i.e. 0, which is often the case).

  Note that the function does not check that the tree is effectively
  ultrametric, so if it is not, the returned result may not be meaningful.
}
\value{
  An object of class \code{"birthdeath"} which is a list with the
  following components:
  \item{tree}{the name of the tree analysed.}
  \item{N}{the number of species.}
  \item{dev}{the deviance (= -2 log lik) at its minimum.}
  \item{para}{the estimated parameters.}
  \item{se}{the corresponding standard-errors.}
  \item{CI}{the 95\% profile-likelihood confidence intervals.}
}
\references{
  Kendall, D. G. (1948) On the generalized ``birth-and-death''
  process. \emph{Annals of Mathematical Statistics}, \bold{19}, 1--15.

  Nee, S., May, R. M. and Harvey, P. H. (1994) The reconstructed
  evolutionary process. \emph{Philosophical Transactions of the Royal
    Society of London. Series B. Biological Sciences}, \bold{344}, 305--311.

  Nee, S., Holmes, E. C., May, R. M. and Harvey, P. H. (1995) Estimating
  extinctions from molecular phylogenies. in \emph{Extinction Rates},
  eds. Lawton, J. H. and May, R. M., pp. 164--182, Oxford University Press.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{branching.times}}, \code{\link{diversi.gof}},
  \code{\link{diversi.time}}, \code{\link{ltt.plot}},
  \code{\link{yule}}, \code{\link{bd.ext}}, \code{\link{yule.cov}},
  \code{\link{bd.time}}
}
\keyword{models}
