\name{chronopl}
\alias{chronopl}
\title{Molecular Dating With Penalized Likelihood}
\description{
  This function estimates the node ages of a tree using a
  semi-parametric method based on penalized likelihood (Sanderson
  2002). The branch lengths of the input tree are interpreted as mean
  numbers of substitutions (i.e., per site).
}
\usage{
chronopl(phy, lambda, age.min = 1, age.max = NULL,
         node = "root", S = 1, tol = 1e-8,
         CV = FALSE, eval.max = 500, iter.max = 500, ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{lambda}{value of the smoothing parameter.}
  \item{age.min}{numeric values specifying the fixed node ages (if
    \code{age.max = NULL}) or the youngest bound of the nodes known to
    be within an interval.}
  \item{age.max}{numeric values specifying the oldest bound of the nodes
    known to be within an interval.}
  \item{node}{the numbers of the nodes whose ages are given by
    \code{age.min}; \code{"root"} is a short-cut for the root.}
  \item{S}{the number of sites in the sequences; leave the default if
    branch lengths are in mean number of substitutions.}
  \item{tol}{the value below which branch lengths are considered
    effectively zero.}
  \item{CV}{whether to perform cross-validation.}
  \item{eval.max}{the maximal number of evaluations of the penalized
    likelihood function.}
  \item{iter.max}{the maximal number of iterations of the optimization
    algorithm.}
  \item{\dots}{further arguments passed to control \code{nlminb}.}
}
\details{
  The idea of this method is to use a trade-off between a parametric
  formulation where each branch has its own rate, and a nonparametric
  term where changes in rates are minimized between contiguous
  branches. A smoothing parameter (lambda) controls this trade-off. If
  lambda = 0, then the parametric component dominates and rates vary as
  much as possible among branches, whereas for increasing values of
  lambda, the variation are smoother to tend to a clock-like model (same
  rate for all branches).

  \code{lambda} must be given. The known ages are given in
  \code{age.min}, and the correponding node numbers in \code{node}.
  These two arguments must obviously be of the same length. By default,
  an age of 1 is assumed for the root, and the ages of the other nodes
  are estimated.

  If \code{age.max = NULL} (the default), it is assumed that
  \code{age.min} gives exactly known ages. Otherwise, \code{age.max} and
  \code{age.min} must be of the same length and give the intervals for
  each node. Some node may be known exactly while the others are
  known within some bounds: the values will be identical in both
  arguments for the former (e.g., \code{age.min = c(10, 5), age.max =
    c(10, 6), node = c(15, 18)} means that the age of node 15 is 10
  units of time, and the age of node 18 is between 5 and 6).

  If two nodes are linked (i.e., one is the ancestor of the other) and
  have the same values of \code{age.min} and \code{age.max} (say, 10 and
  15) this will result in an error because the medians of these values
  are used as initial times (here 12.5) giving initial branch length(s)
  equal to zero. The easiest way to solve this is to change slightly the
  given values, for instance use \code{age.max = 14.9} for the youngest
  node, or \code{age.max = 15.1} for the oldest one (or similarly for
  \code{age.min}).

  The input tree may have multichotomies. If some internal branches are
  of zero-length, they are collapsed (with a warning), and the returned
  tree will have less nodes than the input one. The presence of
  zero-lengthed terminal branches of results in an error since it makes
  little sense to have zero-rate branches.

  The cross-validation used here is different from the one proposed by
  Sanderson (2002). Here, each tip is dropped successively and the
  analysis is repeated with the reduced tree: the estimated dates for
  the remaining nodes are compared with the estimates from the full
  data. For the \eqn{i}{i}th tip the following is calculated:

  \deqn{\sum_{j=1}^{n-2}{\frac{(t_j - t_j^{-i})^2}{t_j}}}{SUM[j = 1, ..., n-2] (tj - tj[-i])^2/tj},

  where \eqn{t_j}{tj} is the estimated date for the \eqn{j}{j}th node
  with the full phylogeny, \eqn{t_j^{-i}}{tj[-i]} is the estimated date
  for the \eqn{j}{j}th node after removing tip \eqn{i}{i} from the tree,
  and \eqn{n}{n} is the number of tips.

  The present version uses the \code{\link[stats]{nlminb}} to optimise
  the penalized likelihood function: see its help page for details on
  parameters controlling the optimisation procedure.
}
\value{
  an object of class \code{"phylo"} with branch lengths as estimated by
  the function. There are three or four further attributes:

  \item{ploglik}{the maximum penalized log-likelihood.}
  \item{rates}{the estimated rates for each branch.}
  \item{message}{the message returned by \code{nlminb} indicating
    whether the optimisation converged.}
  \item{D2}{the influence of each observation on overall date
    estimates (if \code{CV = TRUE}).}
}
\note{
  The new function \code{\link{chronos}} replaces the present one which
  is no more maintained.
}
\references{
  Sanderson, M. J. (2002) Estimating absolute rates of molecular
  evolution and divergence times: a penalized likelihood
  approach. \emph{Molecular Biology and Evolution}, \bold{19},
  101--109.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{chronos}}, \code{\link{chronoMPL}}
}
\keyword{models}
