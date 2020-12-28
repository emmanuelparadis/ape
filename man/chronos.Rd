\name{chronos}
\alias{chronos}
\alias{makeChronosCalib}
\alias{chronos.control}
\alias{print.chronos}
\title{Molecular Dating by Penalised Likelihood and Maximum Likelihood}
\description{
  \code{chronos} is the main function fitting a chronogram to a
  phylogenetic tree whose branch lengths are in number of substitution
  per sites.

  \code{makeChronosCalib} is a tool to prepare data frames with the
  calibration points of the phylogenetic tree.

  \code{chronos.control} creates a list of parameters to be passed
  to \code{chronos}.
}
\usage{
chronos(phy, lambda = 1, model = "correlated", quiet = FALSE,
        calibration = makeChronosCalib(phy),
        control = chronos.control())
\method{print}{chronos}(x, ...)
makeChronosCalib(phy, node = "root", age.min = 1,
   age.max = age.min, interactive = FALSE, soft.bounds = FALSE)
chronos.control(...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{lambda}{value of the smoothing parameter.}
  \item{model}{a character string specifying the model of substitution
    rate variation among branches. The possible choices are:
    ``correlated'', ``relaxed'', ``discrete'', or an unambiguous
    abbreviation of these.}
  \item{quiet}{a logical value; by default the calculation progress are
    displayed.}
  \item{calibration}{a data frame (see details).}
  \item{control}{a list of parameters controlling the optimisation
    procedure (see details).}
  \item{x}{an object of class \code{c("chronos", "phylo")}.}
  \item{node}{a vector of integers giving the node numbers for which a
    calibration point is given. The default is a short-cut for the
    root.}
  \item{age.min, age.max}{vectors of numerical values giving the minimum
    and maximum ages of the nodes specified in \code{node}.}
  \item{interactive}{a logical value. If \code{TRUE}, then \code{phy} is
    plotted and the user is asked to click close to a node and enter the
    ages on the keyboard.}
  \item{soft.bounds}{(currently unused)}
  \item{\dots}{in the case of \code{chronos.control}: one of the five
    parameters controlling optimisation (unused in the case of
    \code{print.chronos}).}
}
\details{
  \code{chronos} replaces \code{chronopl} but with a different interface
  and some extensions (see References).

  The known dates (argument \code{calibration}) must be given in a data
  frame with the following column names: node, age.min, age.max, and
  soft.bounds (the last one is yet unused). For each row, these are,
  respectively: the number of the node in the ``phylo'' coding standard,
  the minimum age for this node, the maximum age, and a logical value
  specifying whether the bounds are soft. If age.min = age.max, this
  means that the age is exactly known. This data frame can be built with
  \code{makeChronosCalib} which returns by default a data frame with a
  single row giving age = 1 for the root. The data frame can be built
  interactively by clicking on the plotted tree.

  The argument \code{control} allows one to change some parameters of
  the optimisation procedure. This must be a list with names. The
  available options with their default values are:

  \itemize{
    \item{tol = 1e-8: }{tolerance for the estimation of the substitution
      rates.}
    \item{iter.max = 1e4: }{the maximum number of iterations at each
      optimization step.}
    \item{eval.max = 1e4: }{the maximum number of function evaluations at
      each optimization step.}
    \item{nb.rate.cat = 10: }{the number of rate categories if \code{model
	= "discrete"} (set this parameter to 1 to fit a strict clock
      model).}
    \item{dual.iter.max = 20: }{the maximum number of alternative
      iterations between rates and dates.}
    \item{epsilon = 1e-6: }{the convergence diagnostic criterion.}
  }

  The command \code{chronos.control()} returns a list with the default
  values of these parameters. They may be modified by passing them to
  this function, or directly in the list.
}
\value{
  \code{chronos} returns an object of class \code{c("chronos",
  "phylo")}. There is a print method for it. There are additional
  attributes which can be visualised with \code{str} or extracted with
  \code{attr}.

  \code{makeChronosCalib} returns a data frame.

  \code{chronos.control} returns a list.
}
\references{
  Kim, J. and Sanderson, M. J. (2008) Penalized likelihood phylogenetic
  inference: bridging the parsimony-likelihood gap. \emph{Systematic
    Biology}, \bold{57}, 665--674.

  Paradis, E. (2013) Molecular dating of phylogenies by likelihood
  methods: a comparison of models and a new information
  criterion. \emph{Molecular Phylogenetics and Evolution}, \bold{67},
  436--444.

  Sanderson, M. J. (2002) Estimating absolute rates of molecular
  evolution and divergence times: a penalized likelihood
  approach. \emph{Molecular Biology and Evolution}, \bold{19},
  101--109.
}
\author{Emmanuel Paradis, Santiago Claramunt, Guillaume Louvel}
\seealso{\code{\link{chronoMPL}}}
\examples{
tr <- rtree(10)
### the default is the correlated rate model:
chr <- chronos(tr)
### strict clock model:
ctrl <- chronos.control(nb.rate.cat = 1)
chr.clock <- chronos(tr, model = "discrete", control = ctrl)
### How different are the rates?
attr(chr, "rates")
attr(chr.clock, "rates")
\dontrun{
cal <- makeChronosCalib(tr, interactive = TRUE)
cal
### if you made mistakes, you can edit the data frame with:
### fix(cal)
chr <- chronos(tr, calibration = cal)
}
}
\keyword{models}
