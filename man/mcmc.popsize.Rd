\name{mcmc.popsize}
\alias{mcmc.popsize}
\alias{extract.popsize}
\alias{plot.popsize}
\alias{lines.popsize}
\title{Reversible Jump MCMC to Infer Demographic History}
\usage{
mcmc.popsize(tree,nstep, thinning=1, burn.in=0,progress.bar=TRUE,
    method.prior.changepoints=c("hierarchical", "fixed.lambda"), max.nodes=30,
   lambda=0.5, gamma.shape=0.5, gamma.scale=2,
    method.prior.heights=c("skyline", "constant", "custom"),
    prior.height.mean,
    prior.height.var)
extract.popsize(mcmc.out, credible.interval=0.95, time.points=200, thinning=1, burn.in=0)
\method{plot}{popsize}(x, show.median=TRUE, show.years=FALSE,
             subst.rate, present.year, xlab = NULL,
             ylab = "Effective population size", log = "y", ...)
\method{lines}{popsize}(x, show.median=TRUE,show.years=FALSE, subst.rate, present.year, ...)

}
\arguments{
  \item{tree}{Either an ultrametric tree (i.e. an object of class \code{"phylo"}),
           or coalescent intervals (i.e. an object of class \code{"coalescentIntervals"}). }
  \item{nstep}{Number of MCMC steps, i.e. length of the Markov chain (suggested value: 10,000-50,000).}
  \item{thinning}{Thinning factor (suggest value: 10-100).}
  \item{burn.in}{Number of steps dropped from the chain to allow for a burn-in phase (suggest value: 1000).}

  \item{progress.bar}{Show progress bar during the MCMC run.}

  \item{method.prior.changepoints}{If \code{hierarchical}is chosen (the default) then the smoothing parameter lambda is drawn from
     a gamma distribution with some specified shape and scale parameters.
     Alternatively, for \code{fixed.lambda} the value of lambda is   a given constant.
  }

  \item{max.nodes}{Upper limit for the number of internal nodes of the approximating spline (default: 30).}
  \item{lambda}{Smoothing parameter. For \code{method="fixed.lambda"} the specifed value of lambda determines
      the mean of the prior distribution   for the number of internal nodes of the approximating
      spline for the demographic function (suggested value: 0.1-1.0).}
  \item{gamma.shape}{Shape parameter of the gamma function from which \code{lambda} is drawn for
    \code{method="hierarchical"}.}
   \item{gamma.scale}{Scale parameter of the gamma function from which \code{lambda} is drawn for
    \code{method="hierarchical"}.}
  \item{method.prior.heights}{Determines the prior for the heights of the change points.
          If \code{custom} is chosen then two functions describing the mean and variance
	  of the heigths in depence of time have to be specified (via \code{prior.height.mean}
	  and \code{prior.height.var} options).  Alternatively, two built-in priors are available:
	  \code{constant} assumes constant population size and variance determined by Felsenstein
	  (1992), and \code{skyline} assumes a skyline plot (see Opgen-Rhein et al. 2004 for
	  more details).}
  \item{prior.height.mean}{Function describing the mean of the prior distribution for the heights
                           (only used if \code{method.prior.heights = custom}).}

  \item{prior.height.var}{Function describing the variance of the prior distribution for the heights
                           (only used if \code{method.prior.heights = custom}).}
  \item{mcmc.out}{Output from \code{mcmc.popsize} - this is needed as input for \code{extract.popsize}.}
 \item{credible.interval}{Probability mass of the confidence band (default: 0.95).}

 \item{time.points}{Number of discrete time points in the table output by \code{extract.popsize}.}

 \item{x}{Table with population size versus time, as computed by \code{extract.popsize}. }

 \item{show.median}{Plot median rather than mean as point estimate for demographic function (default: TRUE).}

  \item{show.years}{Option that determines whether the time is plotted in units of
    of substitutions (default) or in years (requires specification of substution rate
    and year of present).}
  \item{subst.rate}{Substitution rate (see option show.years).}
  \item{present.year}{Present year (see option show.years).}
  \item{xlab}{label on the x-axis (depends on the value of
    \code{show.years}).}
  \item{ylab}{label on the y-axis.}
  \item{log}{log-transformation of axes; by default, the y-axis is
    log-transformed.}
  \item{\dots}{Further arguments to be passed on  to \code{plot} or \code{lines}.}
}
\description{
 These functions implement a reversible jump MCMC framework to infer the demographic history,
 as well as corresponding confidence bands,
 from a genealogical tree. The computed demographic history is a continous
 and smooth function in time.
 \code{mcmc.popsize} runs the actual MCMC chain and outputs information about the
 sampling steps, \code{extract.popsize} generates from this MCMC
 output a table of population size in time, and  \code{plot.popsize} and \code{lines.popsize}
 provide utility functions to plot the corresponding demographic functions.
}

\details{
 Please refer to Opgen-Rhein et al. (2005) for methodological details, and the help page of
 \code{\link{skyline}} for information on a related approach.
}
\author{
  Rainer Opgen-Rhein and Korbinian Strimmer. Parts of the rjMCMC
  sampling procedure are adapted from \R code by Karl Broman.
}
\seealso{
\code{\link{skyline}} and \code{\link{skylineplot}}. }
\references{
  Opgen-Rhein, R., Fahrmeir, L. and Strimmer, K. 2005. Inference of
  demographic history from genealogical trees using reversible jump
  Markov chain Monte Carlo. \emph{BMC Evolutionary Biology}, \bold{5},
  6.
}
\examples{
# get tree
data("hivtree.newick") # example tree in NH format
tree.hiv <- read.tree(text = hivtree.newick) # load tree

# run mcmc chain
mcmc.out <- mcmc.popsize(tree.hiv, nstep=100, thinning=1, burn.in=0,progress.bar=FALSE) # toy run
#mcmc.out <- mcmc.popsize(tree.hiv, nstep=10000, thinning=5, burn.in=500) # remove comments!!

# make list of population size versus time
popsize  <- extract.popsize(mcmc.out)

# plot and compare with skyline plot
sk <- skyline(tree.hiv)
plot(sk, lwd=1, lty=3, show.years=TRUE, subst.rate=0.0023, present.year = 1997)
lines(popsize, show.years=TRUE, subst.rate=0.0023, present.year = 1997)
}
\keyword{manip}
