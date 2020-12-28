\name{ltt.plot}
\alias{ltt.plot}
\alias{ltt.lines}
\alias{mltt.plot}
\alias{ltt.coplot}
\alias{ltt.plot.coords}
\title{Lineages Through Time Plot}
\description{
  These functions provide tools for plotting the numbers of lineages
  through time from phylogenetic trees.
}
\usage{
ltt.plot(phy, xlab = "Time", ylab = "N",
         backward = TRUE, tol = 1e-6, ...)
ltt.lines(phy, backward = TRUE, tol = 1e-6, ...)
mltt.plot(phy, ..., dcol = TRUE, dlty = FALSE, legend = TRUE,
          xlab = "Time", ylab = "N", log = "", backward = TRUE,
          tol = 1e-6)
ltt.coplot(phy, backward = TRUE, ...)
ltt.plot.coords(phy, backward = TRUE, tol = 1e-6)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}; this could be an object
    of class \code{"multiPhylo"} in the case of \code{mltt.plot}.}
  \item{xlab}{a character string (or a variable of mode character)
    giving the label for the \eqn{x}-axis (default is "Time").}
  \item{ylab}{idem for the \eqn{y}-axis (default is "N").}
  \item{backward}{a logical value: should the time axis be traced from
    the present (the default), or from the root of the tree?}
  \item{tol}{a numeric value (see details).}
  \item{\dots}{in the cases of \code{ltt.plot()}, \code{ltt.lines()},
    or \code{ltt.coplot()} these are further (graphical) arguments to be
    passed to \code{plot()}, \code{lines()}, or \code{plot.phylo()},
    respectively (see details on how to transform the axes); in
    the case of \code{mltt.plot()} these are additional trees to be plotted
    (see details).}
  \item{dcol}{a logical specifying whether the different curves should
    be differentiated with colors (default is \code{TRUE}).}
  \item{dlty}{a logical specifying whether the different curves should
    be differentiated with patterns of dots and dashes (default is
    \code{FALSE}).}
  \item{legend}{a logical specifying whether a legend should be
    plotted.}
  \item{log}{a character string specifying which axis(es) to be
    log-transformed; must be one of the followings: \code{""},
    \code{"x"}, \code{"y"}, or \code{"xy"}.}
}
\details{
  \code{ltt.plot} does a simple lineages through time (LTT)
  plot. Additional arguments (\code{\dots}) may be used to change, for
  instance, the limits on the axes (with \code{xlim} and/or
  \code{ylim}) or other graphical settings (\code{col} for the color,
  \code{lwd} for the line thickness, \code{lty} for the line type may be
  useful; see \code{\link[graphics]{par}} for an exhaustive listing of
  graphical parameters). The \eqn{y}-axis can be log-transformed by
  adding the following option: \code{log = "y"}.

  The option \code{tol} is used as follows: first the most distant tip
  from the root is found, then all tips whose distance to the root is
  not different from the previous one more than \code{tol} are
  considered to be contemporaneous with it.

  If the tree is not ultrametric, the plot is done assuming the tips,
  except the most distant from the root, represent extinction events. If
  a root edge is present, it is taken into account.

  \code{ltt.lines} adds a LTT curve to an existing plot. Additional
  arguments (\code{\dots}) may be used to change the settings of the added
  line.

  \code{mltt.plot} does a multiple LTT plot taking as arguments one or
  several trees. These trees may be given as objects of class
  \code{"phylo"} (single trees) and/or \code{"multiPhylo"} (multiple
  trees). Any number of objects may be given. This function is mainly
  for exploratory analyses with the advantages that the axes are set
  properly to view all lines, and the legend is plotted by default. The
  plot will certainly make sense if all trees have their
  most-distant-from-the-root tips contemporaneous (i.e., trees with only
  extinct lineages will not be represented properly). For more flexible
  settings of line drawings, it may be better to combine
  \code{ltt.plot()} with successive calls of \code{ltt.lines()} (see
  examples).

  \code{ltt.coplot} is meant to show how to set a tree and a LTT plots
  on the same scales. All extra arguments modify only the appearance of
  the tree. The code can be easily edited and tailored.
}
\value{
  \code{ltt.plot.coords} returns a two-column matrix with the time
  points and the number of lineages, respectively. The \eqn{i}th value
  of the second column is the number of lineages for the interval
  defined by the \eqn{(i - 1)}th and the \eqn{i}th values of the first
  column. These are then plotted with the option \code{type = "S"} by
  the other functions.
}
\references{
  Harvey, P. H., May, R. M. and Nee, S. (1994) Phylogenies without
  fossils. \emph{Evolution}, \bold{48}, 523--529.

  Nee, S., Holmes, E. C., Rambaut, A. and Harvey, P. H. (1995) Inferring
  population history from molecular phylogenies. \emph{Philosophical
    Transactions of the Royal Society of London. Series B. Biological
    Sciences}, \bold{349}, 25--31.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{kronoviz}}, \code{\link{skyline}}, \code{\link{LTT}},
  \code{\link{branching.times}}, \code{\link{birthdeath}},
  \code{\link{bd.ext}}, \code{\link{yule.cov}}, \code{\link{bd.time}};
  \code{\link[graphics]{plot}} for the basic plotting function in R
}
\examples{
data(bird.families)
opar <- par(mfrow = c(2, 1))
ltt.plot(bird.families)
title("Lineages Through Time Plot of the Bird Families")
ltt.plot(bird.families, log = "y")
title(main = "Lineages Through Time Plot of the Bird Families",
      sub = "(with logarithmic transformation of the y-axis)")
par(opar)

### to plot the tree and the LTT plot together
data(bird.orders)
layout(matrix(1:4, 2, 2))
plot(bird.families, show.tip.label = FALSE)
ltt.plot(bird.families, main = "Bird families")
plot(bird.orders, show.tip.label = FALSE)
ltt.plot(bird.orders, main = "Bird orders")
layout(1)

### better with ltt.coplot():
ltt.coplot(bird.families, show.tip.label = FALSE, x.lim = 27.5)
data(chiroptera)
chiroptera <- compute.brlen(chiroptera)
ltt.coplot(chiroptera, show.tip.label = FALSE, type = "c")

### with extinct lineages and a root edge:
omar <- par("mar")
set.seed(31)
tr <- rlineage(0.2, 0.15)
tr$root.edge <- 5
ltt.coplot(tr, show.tip.label = FALSE, x.lim = 55)
## compare with:
## ltt.coplot(drop.fossil(tr), show.tip.label = FALSE)
layout(1)
par(mar = omar)

mltt.plot(bird.families, bird.orders)
### Generates 10 random trees with 23 tips:
TR <- replicate(10, rcoal(23), FALSE)
### Give names to each tree:
names(TR) <- paste("random tree", 1:10)
### And specify the class of the list so that mltt.plot()
### does not trash it!
class(TR) <- "multiPhylo"
mltt.plot(TR, bird.orders)
### And now for something (not so) completely different:
ltt.plot(bird.orders, lwd = 2)
for (i in 1:10) ltt.lines(TR[[i]], lty = 2)
legend(-20, 10, lwd = c(2, 1), lty = c(1, 2), bty = "n",
       legend = c("Bird orders", "Random (coalescent) trees"))
}
\keyword{hplot}
\keyword{aplot}
