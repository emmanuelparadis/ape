\name{diversi.time}
\alias{diversi.time}
\title{Analysis of Diversification with Survival Models}
\usage{
diversi.time(x, census = NULL, censoring.codes = c(1, 0), Tc = NULL)
}
\arguments{
  \item{x}{a numeric vector with the branching times.}
  \item{census}{a vector of the same length than `x' used as an
    indicator variable; thus, it must have only two values, one coding
    for accurately known branching times, and the other for censored
    branching times. This argument can be of any mode (numeric, character,
    logical), or can even be a factor.}
  \item{censoring.codes}{a vector of length two giving the codes used
    for \code{census}: by default 1 (accurately known times) and 0 (censored
    times). The mode must be the same than the one of \code{census}.}
  \item{Tc}{a single numeric value specifying the break-point time to
    fit Model C. If none is provided, then it is set arbitrarily to the
    mean of the analysed branching times.}
}
\description{
  This functions fits survival models to a set of branching times, some
  of them may be known approximately (censored). Three models are
  fitted, Model A assuming constant diversification, Model B assuming
  that diversification follows a Weibull law, and Model C assuming that
  diversification changes with a breakpoint at time `Tc'. The models are
  fitted by maximum likelihood.
}
\details{
  The principle of the method is to consider each branching time as an
  event: if the branching time is accurately known, then it is a failure
  event; if it is approximately knwon then it is a censoring event. An
  analogy is thus made between the failure (or hazard) rate estimated by
  the survival models and the diversification rate of the lineage. Time
  is here considered from present to past.

  Model B assumes a monotonically changing diversification rate. The
  parameter that controls the change of this rate is called beta. If
  beta is greater than one, then the diversification rate decreases
  through time; if it is lesser than one, the the rate increases through
  time. If beta is equal to one, then Model B reduces to Model A.
}
\value{
  A NULL value is returned, the results are simply printed.
}
\references{
  Paradis, E. (1997) Assessing temporal variations in diversification
  rates from phylogenies: estimation and hypothesis
  testing. \emph{Proceedings of the Royal Society of London. Series
    B. Biological Sciences}, \bold{264}, 1141--1147.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{branching.times}}, \code{\link{diversi.gof}}
  \code{\link{ltt.plot}}, \code{\link{birthdeath}},
  \code{\link{bd.ext}}, \code{\link{yule}}, \code{\link{yule.cov}}
}
\keyword{models}
