\name{LTT}
\alias{LTT}
\title{Theoretical Lineage-Through Time Plots}
\description{
  This function draws the lineage-through time (LTT) plots predicted
  under a speciation-extinction model (aka birth-death model) with
  specified values of speciation and extinction rates (which may vary
  with time).

  A prediction interval is plotted by default which requires to define a
  sample size (100 by default), and different curves can be combined.
}
\usage{
LTT(birth = 0.1, death = 0, N = 100, Tmax = 50, PI = 95,
    scaled = TRUE, eps = 0.1, add = FALSE, backward = TRUE,
    ltt.style = list("black", 1, 1), pi.style = list("blue", 1, 2), ...)
}
\arguments{
  \item{birth}{the speciation rate, this may be either a numeric value
    or a funtion of time (named \code{t} in the code of the function).}
  \item{death}{id. for the extinction rate.}
  \item{N}{the size of the tree.}
  \item{Tmax}{the age of the root of the tree.}
  \item{PI}{the percentage value of the prediction interval; set this
    value to 0 to not draw this interval.}
  \item{scaled}{a logical values specifying whether to scale the
    \eqn{y}-axis between 0 and 1.}
  \item{eps}{a numerical value giving the resolution of the time axis.}
  \item{add}{a logical values specifying whether to make a new plot (the
    default).}
  \item{backward}{a logical value: should the time axis be traced from
    the present (the default), or from the root of the tree?}
  \item{ltt.style}{a list with three elements giving the style of the
    LTT curve with, respectively, the colour (\code{"col"}), the line
    thickness (\code{"lwd"}), and the line type (\code{"lty"}).}
  \item{pi.style}{id. for the prediction interval.}
  \item{\dots}{arguments passed to \code{plot} (e.g., \code{log="y"}).}
}
\details{
  For the moment, this works well when \code{birth} and \code{death} are
  constant. Some improvements are under progress for time-dependent
  rates (but see below for an example).
}
\references{
  Hallinan, N. (2012) The generalized time variable reconstructed
  birth--death process. \emph{Journal of Theoretical Biology},
  \bold{300}, 265--276.

  Paradis, E. (2011) Time-dependent speciation and extinction from
  phylogenies: a least squares approach. \emph{Evolution}, \bold{65},
  661--672.

  Paradis, E. (2015) Random phylogenies and the distribution of
  branching times. \emph{Journal of Theoretical Biology}, \bold{387},
  39--45.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{ltt.plot}}
}
\examples{
### predicted LTT plot under a Yule model with lambda = 0.1
### and 50 species after 50 units of time...
LTT(N = 50)
### ... and with a birth-death model with the same rate of
### diversification (try with N = 500):
LTT(0.2, 0.1, N = 50, PI = 0, add = TRUE, ltt.style = list("red", 2, 1))
### predictions under different tree sizes:
layout(matrix(1:4, 2, 2, byrow = TRUE))
for (N in c(50, 100, 500, 1000)) {
    LTT(0.2, 0.1, N = N)
    title(paste("N =", N))
}
layout(1)
\dontrun{
### speciation rate decreasing with time
birth.logis <- function(t) 1/(1 + exp(0.02 * t + 4))
LTT(birth.logis)
LTT(birth.logis, 0.05)
LTT(birth.logis, 0.1)
}
}
\keyword{hplot}
