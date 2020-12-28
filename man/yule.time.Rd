\name{yule.time}
\alias{yule.time}
\title{Fits the Time-Dependent Yule Model}
\usage{
yule.time(phy, birth, BIRTH = NULL, root.time = 0, opti = "nlm", start = 0.01)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{birth}{a (vectorized) function specifying how the birth
    (speciation) probability changes through time (see details).}
  \item{BIRTH}{a (vectorized) function giving the primitive of
    \code{birth}.}
  \item{root.time}{a numeric value giving the time of the root node
    (time is measured from the past towards the present).}
  \item{opti}{a character string giving the function used for
    optimisation of the likelihood function. Three choices are possible:
    \code{"nlm"}, \code{"nlminb"}, or \code{"optim"}, or any unambiguous
    abbreviation of these.}
  \item{start}{the initial values used in the optimisation.}
}
\description{
  This function fits by maximum likelihood the time-dependent Yule
  model. The time is measured from the past (\code{root.time}) to the
  present.
}
\details{
  The model fitted is a straightforward extension of the Yule model with
  covariates (see \code{\link{yule.cov}}). Rather than having
  heterogeneity among lineages, the speciation probability is the same
  for all lineages at a given time, but can change through time.

  The function \code{birth} \emph{must} meet these two requirements: (i)
  the parameters to be estimated are the formal arguments; (ii) time is
  named \code{t} in the body of the function. However, this is the
  opposite for the primitive \code{BIRTH}: \code{t} is the formal
  argument, and the parameters are used in its body. See the examples.

  It is recommended to use \code{BIRTH} if possible, and required if
  speciation probability is constant on some time interval. If this
  primitive cannot be provided, a numerical integration is done with
  \code{\link[stats]{integrate}}.

  The standard-errors of the parameters are computed with the Hessian of
  the log-likelihood function.
}
\value{
  An object of class \code{"yule"} (see \code{\link{yule}}).
}
\author{Emmanuel Paradis}
\references{
  Hubert, N., Paradis, E., Bruggemann, H. and Planes, S. (2011) Community
  assembly and diversification in Indo-Pacific coral reef
  fishes. \emph{Ecology and Evolution}, \bold{1}, 229--277.
}
\seealso{
  \code{\link{branching.times}}, \code{\link{ltt.plot}},
  \code{\link{birthdeath}}, \code{\link{yule}}, \code{\link{yule.cov}},
  \code{\link{bd.time}}
}
\examples{
### define two models...
birth.logis <- function(a, b) 1/(1 + exp(-a*t - b)) # logistic
birth.step <- function(l1, l2, Tcl) { # 2 rates with one break-point
    ans <- rep(l1, length(t))
    ans[t > Tcl] <- l2
    ans
}
### ... and their primitives:
BIRTH.logis <- function(t) log(exp(-a*t) + exp(b))/a + t
BIRTH.step <- function(t)
{
    out <- numeric(length(t))
    sel <- t <= Tcl
    if (any(sel)) out[sel] <- t[sel] * l1
    if (any(!sel)) out[!sel] <- Tcl * l1 + (t[!sel] - Tcl) * l2
    out
}
data(bird.families)
### fit both models:
yule.time(bird.families, birth.logis)
yule.time(bird.families, birth.logis, BIRTH.logis) # same but faster
\dontrun{yule.time(bird.families, birth.step)}  # fails
yule.time(bird.families, birth.step, BIRTH.step,
          opti = "nlminb", start = c(.01, .01, 100))
}
\keyword{models}
