\name{rtt}
\alias{rtt}
\title{Root a Tree by Root-to-Tip Regression}
\description{
  This function roots a phylogenetic tree with dated tips in the location
  most compatible with the assumption of a strict molecular clock.
}
\usage{
rtt(t, tip.dates, ncpu = 1, objective = correlation,
    opt.tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{t}{an object of class \code{"phylo"}.}
  \item{tip.dates}{a vector of sampling times associated to the tips of
    \code{t}, in the same order as \code{t$tip.label}.}
  \item{ncpu}{number of cores to use.}
  \item{objective}{one of \code{"correlation"}, \code{"rms"}, or
    \code{"rsquared"}.}
  \item{opt.tol}{tolerance for optimization precision.}
}
\details{
  This function duplicates one part the functionality of the program
  Path-O-Gen (see references). The root position is chosen to produce
  the best linear regression of root-to-tip distances against sampling
  times.

  \code{t} must have branch lengths in units of expected substitutions
  per site.

  \code{tip.dates} should be a vector of sampling times, in any time
  unit, with time increasing toward the present. For example, this may
  be in units of ``days since study start'' or ``years since 10,000
  BCE'', but not ``millions of yearsago''.

  Setting \code{ncpu} to a value larger than 1 requires the \code{parallel}
  library.

  \code{objective} is the measure which will be used to define the
  ``goodness'' of a regression fit. It may be one of \code{"correlation"}
  (strongest correlation between tip date and distance from root),
  \code{"rms"} (lowest root-mean-squared error), or \code{"rsquared"}
  (highest R-squared value).

  \code{opt.tol} is used to optimize the location of the root along the best
  branch. By default, R's \code{optimize} function uses a precision of
  \code{.Machine$double.eps^0.25}, which is about 0.0001 on a 64-bit system.
  This should be set to a smaller value if the branch lengths of \code{t} are
  very short.
}
\value{
  an object of class \code{"phylo"}.
}
\note{
  This function only chooses the best root. It does not rescale the branch
  lengths to time, or perform a statistical test of the molecular clock
  hypothesis.
}
\author{
  Rosemary McCloskey\email{rmccloskey@cfenet.ubc.ca},
  Emmanuel Paradis
}
\references{
  Rambaut, A. (2009). Path-O-Gen: temporal signal investigation
  tool.

  Rambaut, A. (2000). Estimating the rate of molecular evolution:
  incorporating non-contemporaneous sequences into maximum likelihood
  phylogenies. \emph{Bioinformatics}, \bold{16}, 395-399.
}
\examples{
t <- rtree(100)
tip.date <- rnorm(t$tip.label)^2
rtt(t, tip.date)
}
