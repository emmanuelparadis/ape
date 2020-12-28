\name{skyline}
\alias{skyline}
\alias{skyline.phylo}
\alias{skyline.coalescentIntervals}
\alias{skyline.collapsedIntervals}
\alias{find.skyline.epsilon}
\title{Skyline Plot Estimate of Effective Population Size}
\usage{
skyline(x, \dots)
\method{skyline}{phylo}(x, \dots)
\method{skyline}{coalescentIntervals}(x, epsilon=0, \dots)
\method{skyline}{collapsedIntervals}(x, old.style=FALSE, \dots)
find.skyline.epsilon(ci, GRID=1000, MINEPS=1e-6, \dots)
}
\arguments{
  \item{x}{Either an ultrametric tree (i.e. an object of class
    \code{"phylo"}), or coalescent intervals (i.e. an object of class
    \code{"coalescentIntervals"}), or collapsed coalescent intervals
    (i.e. an object of class \code{"collapsedIntervals"}).}
  \item{epsilon}{collapsing parameter that controls the amount of smoothing
    (allowed range: from \code{0} to \code{ci$total.depth}, default value: 0). This is the same parameter as in
    \link{collapsed.intervals}.}
  \item{old.style}{Parameter to choose between two slightly different variants of the
     generalized skyline plot (Strimmer and Pybus, pers. comm.). The default value \code{FALSE} is
     recommended.}
  \item{ci}{coalescent intervals (i.e. an object of class \code{"coalescentIntervals"})}
  \item{GRID}{Parameter for the grid search for \code{epsilon} in \code{find.skyline.epsilon}.}
  \item{MINEPS}{Parameter for the grid search for \code{epsilon} in \code{find.skyline.epsilon}.}
  \item{\dots}{Any of the above parameters.}
}
\description{
 \code{skyline} computes the \emph{generalized skyline plot} estimate of effective population size
 from an estimated phylogeny.  The demographic history is approximated by
 a step-function.  The number of parameters of the skyline plot (i.e. its smoothness)
 is controlled by a parameter \code{epsilon}.

 \code{find.skyline.epsilon} searches for an optimal value of the \code{epsilon} parameter,
 i.e. the value that maximizes the AICc-corrected log-likelihood (\code{logL.AICc}).
}
\details{
\code{skyline} implements the \emph{generalized skyline plot}  introduced in
Strimmer and Pybus (2001).  For \code{epsilon = 0} the
generalized skyline plot degenerates to the
\emph{classic skyline plot} described in
Pybus et al. (2000).  The latter is in turn directly related to lineage-through-time plots
(Nee et al., 1995).
}

\value{
\code{skyline} returns an object of class \code{"skyline"} with the following entries:

  \item{time}{ A vector with the time at the end of each coalescent
    interval (i.e. the accumulated interval lengths from the beginning of the first interval
    to the end of an interval)}
  \item{interval.length}{ A vector with the length of each interval.}
  \item{population.size}{A vector with the effective population size of each interval.}
  \item{parameter.count}{ Number of free parameters in the skyline plot.}
  \item{epsilon}{The value of the underlying smoothing parameter.}
  \item{logL}{Log-likelihood of skyline plot (see Strimmer and Pybus, 2001).}
  \item{logL.AICc}{AICc corrected log-likelihood (see Strimmer and Pybus, 2001).}

  \code{find.skyline.epsilon} returns the value of the \code{epsilon} parameter
  that maximizes \code{logL.AICc}.
}
\author{Korbinian Strimmer}
\seealso{
\code{\link{coalescent.intervals}}, \code{\link{collapsed.intervals}},
\code{\link{skylineplot}}, \code{\link{ltt.plot}}.
}
\references{
  Strimmer, K. and Pybus, O. G. (2001) Exploring the demographic history
  of DNA sequences using the generalized skyline plot. \emph{Molecular
    Biology and Evolution}, \bold{18}, 2298--2305.

  Pybus, O. G, Rambaut, A. and Harvey, P. H. (2000) An integrated
  framework for the inference of viral population history from
  reconstructed genealogies. \emph{Genetics}, \bold{155}, 1429--1437.

  Nee, S., Holmes, E. C., Rambaut, A. and Harvey, P. H. (1995) Inferring
  population history from molecular phylogenies. \emph{Philosophical
    Transactions of the Royal Society of London. Series B. Biological
    Sciences}, \bold{349}, 25--31.
}
\examples{
# get tree
data("hivtree.newick") # example tree in NH format
tree.hiv <- read.tree(text = hivtree.newick) # load tree

# corresponding coalescent intervals
ci <- coalescent.intervals(tree.hiv) # from tree

# collapsed intervals
cl1 <- collapsed.intervals(ci,0)
cl2 <- collapsed.intervals(ci,0.0119)

#### classic skyline plot ####
sk1 <- skyline(cl1)        # from collapsed intervals 
sk1 <- skyline(ci)         # from coalescent intervals
sk1 <- skyline(tree.hiv)   # from tree
sk1

plot(skyline(tree.hiv))
skylineplot(tree.hiv) # shortcut

plot(sk1, show.years=TRUE, subst.rate=0.0023, present.year = 1997)

#### generalized skyline plot ####

sk2 <- skyline(cl2)              # from collapsed intervals
sk2 <- skyline(ci, 0.0119)       # from coalescent intervals
sk2 <- skyline(tree.hiv, 0.0119) # from tree
sk2

plot(sk2)

# classic and generalized skyline plot together in one plot
plot(sk1, show.years=TRUE, subst.rate=0.0023, present.year = 1997, col=c(grey(.8),1))
lines(sk2,  show.years=TRUE, subst.rate=0.0023, present.year = 1997)
legend(.15,500, c("classic", "generalized"), col=c(grey(.8),1),lty=1)


# find optimal epsilon parameter using AICc criterion
find.skyline.epsilon(ci)

sk3 <- skyline(ci, -1) # negative epsilon also triggers estimation of epsilon
sk3$epsilon
}
\keyword{manip}
