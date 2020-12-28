\name{bd.time}
\alias{bd.time}
\title{Time-Dependent Birth-Death Models}
\description{
  This function fits a used-defined time-dependent birth-death
  model.
}
\usage{
bd.time(phy, birth, death, BIRTH = NULL, DEATH = NULL,
        ip, lower, upper, fast = FALSE, boot = 0, trace = 0)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{birth}{either a numeric (if speciation rate is assumed
    constant), or a (vectorized) function specifying how the birth
    (speciation) probability changes through time (see details).}
  \item{death}{id. for extinction probability.}
  \item{BIRTH}{(optional) a vectorized function giving the primitive
    of \code{birth}.}
  \item{DEATH}{id. for \code{death}.}
  \item{ip}{a numeric vector used as initial values for the estimation
    procedure. If missing, these values are guessed.}
  \item{lower, upper}{the lower and upper bounds of the parameters. If
    missing, these values are guessed too.}
  \item{fast}{a logical value specifying whether to use faster
    integration (see details).}
  \item{boot}{the number of bootstrap replicates to assess the
    confidence intervals of the parameters. Not run by default.}
  \item{trace}{an integer value. If non-zero, the fitting procedure is
    printed every \code{trace} steps. This can be helpful if convergence
    is particularly slow.}
}
\details{
  Details on how to specify the birth and death functions and their
  primitives can be found in the help page of \code{\link{yule.time}}.

  The model is fitted by minimizing the least squares deviation between
  the observed and the predicted distributions of branching times. These
  computations rely heavily on numerical integrations. If \code{fast =
  FALSE}, integrations are done with R's \code{\link[stats]{integrate}}
  function. If \code{fast = TRUE}, a faster but less accurate function
  provided in \pkg{ape} is used. If fitting a complex model to a large
  phylogeny, a strategy might be to first use the latter option, and
  then to use the estimates as starting values with \code{fast = FALSE}.
}
\value{
  A list with the following components:

\itemize{
  \item{par}{a vector of estimates with names taken from the parameters
    in the specified functions.}
  \item{SS}{the minimized sum of squares.}
  \item{convergence}{output convergence criterion from
    \code{\link[stats]{nlminb}}.}
  \item{message}{id.}
  \item{iterations}{id.}
  \item{evaluations}{id.}
}}
\references{
  Paradis, E. (2011) Time-dependent speciation and extinction from
  phylogenies: a least squares approach. \emph{Evolution}, \bold{65},
  661--672.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{ltt.plot}}, \code{\link{birthdeath}},
  \code{\link{yule.time}}, \code{\link{LTT}}
}
\examples{
set.seed(3)
tr <- rbdtree(0.1, 0.02)
bd.time(tr, 0, 0) # fits a simple BD model
bd.time(tr, 0, 0, ip = c(.1, .01)) # 'ip' is useful here
## the classic logistic:
birth.logis <- function(a, b) 1/(1 + exp(-a*t - b))
\dontrun{
bd.time(tr, birth.logis, 0, ip = c(0, -2, 0.01))
## slow to get:
## $par
##            a            b        death
## -0.003486961 -1.995983179  0.016496454
##
## $SS
## [1] 20.73023
}
}
\keyword{models}
