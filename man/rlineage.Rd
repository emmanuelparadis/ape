\name{rlineage}
\alias{rlineage}
\alias{rbdtree}
\alias{rphylo}
\alias{drop.fossil}
\title{Tree Simulation Under the Time-Dependent Birth--Death Models}
\description{
  These three functions simulate phylogenies under any time-dependent
  birth--death model: \code{rlineage} generates a complete tree including
  the species going extinct before present; \code{rbdtree} generates a
  tree with only the species living at present (thus the tree is
  ultrametric); \code{rphylo} generates a tree with a fixed number of
  species at present time. \code{drop.fossil} is a utility function to
  remove the extinct species.
}
\usage{
rlineage(birth, death, Tmax = 50, BIRTH = NULL,
         DEATH = NULL, eps = 1e-6)
rbdtree(birth, death, Tmax = 50, BIRTH = NULL,
        DEATH = NULL, eps = 1e-6)
rphylo(n, birth, death, BIRTH = NULL, DEATH = NULL,
       T0 = 50, fossils = FALSE, eps = 1e-06)
drop.fossil(phy, tol = 1e-8)
}
\arguments{
  \item{birth, death}{a numeric value or a (vectorized) function
    specifying how speciation and extinction rates vary through time.}
  \item{Tmax}{a numeric value giving the length of the simulation.}
  \item{BIRTH, DEATH}{a (vectorized) function which is the primitive
    of \code{birth} or \code{death}. This can be used to speed-up the
    computation. By default, a numerical integration is done.}
  \item{eps}{a numeric value giving the time resolution of the
    simulation; this may be increased (e.g., 0.001) to shorten
    computation times.}
  \item{n}{the number of species living at present time.}
  \item{T0}{the time at present (for the backward-in-time algorithm).}
  \item{fossils}{a logical value specifying whether to output the
    lineages going extinct.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{tol}{a numeric value giving the tolerance to consider a species
    as extinct.}
}
\details{
  These three functions use continuous-time algorithms: \code{rlineage}
  and \code{rbdtree} use the forward-in-time algorithms described in
  Paradis (2011), whereas \code{rphylo} uses a backward-in-time
  algorithm from Stadler (2011). The models are time-dependent
  birth--death models as  described in Kendall (1948). Speciation
  (birth) and extinction (death) rates may be constant or vary through
  time according to an \R function specified by the user. In the latter
  case, \code{BIRTH} and/or \code{DEATH} may be used if the primitives
  of \code{birth} and \code{death} are known. In these functions time is
  the formal argument and must be named \code{t}.

  Note that \code{rphylo} simulates trees in a way similar to what
  the package \pkg{TreeSim} does, the difference is in the
  parameterization of the time-dependent models which is here the same
  than used in the two other functions. In this parameterization scheme,
  time is measured from past to present (see details in Paradis 2015
  which includes a comparison of these algorithms).

  The difference between \code{rphylo} and \code{rphylo(... fossils
    = TRUE)} is the same than between \code{rbdtree} and \code{rlineage}.
}
\value{
  An object of class \code{"phylo"}.
}
\references{
  Kendall, D. G. (1948) On the generalized ``birth-and-death''
  process. \emph{Annals of Mathematical Statistics}, \bold{19}, 1--15.

  Paradis, E. (2011) Time-dependent speciation and extinction from
  phylogenies: a least squares approach. \emph{Evolution}, \bold{65},
  661--672.

  Paradis, E. (2015) Random phylogenies and the distribution of
  branching times. \emph{Journal of Theoretical Biology}, \bold{387},
  39--45.

  Stadler, T. (2011) Simulating trees with a fixed number of extant
  species. \emph{Systematic Biology}, \bold{60}, 676--684.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{yule}}, \code{\link{yule.time}}, \code{\link{birthdeath}},
  \code{\link{rtree}}, \code{\link{stree}}
}
\examples{
set.seed(10)
plot(rlineage(0.1, 0)) # Yule process with lambda = 0.1
plot(rlineage(0.1, 0.05)) # simple birth-death process
b <- function(t) 1/(1 + exp(0.2*t - 1)) # logistic
layout(matrix(0:3, 2, byrow = TRUE))
curve(b, 0, 50, xlab = "Time", ylab = "")
mu <- 0.07
segments(0, mu, 50, mu, lty = 2)
legend("topright", c(expression(lambda), expression(mu)),
       lty = 1:2, bty = "n")
plot(rlineage(b, mu), show.tip.label = FALSE)
title("Simulated with 'rlineage'")
plot(rbdtree(b, mu), show.tip.label = FALSE)
title("Simulated with 'rbdtree'")
}
\keyword{datagen}
