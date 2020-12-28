\name{phydataplot}
\alias{phydataplot}
\alias{ring}
\title{Tree Annotation}
\description{
  \code{phydataplot} plots data on a tree in a way that adapts to the
  type of tree. \code{ring} does the same for circular trees.

  Both functions match the data with the labels of the tree.
}
\usage{
phydataplot(x, phy, style = "bars", offset = 1, scaling = 1,
            continuous = FALSE, width = NULL, legend = "below",
            funcol = rainbow, ...)
ring(x, phy, style = "ring", offset = 1, ...)
}
\arguments{
  \item{x}{a vector, a factor, a matrix, or a data frame.}
  \item{phy}{the tree (which must be already plotted).}
  \item{style}{a character string specifying the type of graphics; can
    be abbreviated (see details).}
  \item{offset}{the space between the tips of the tree and the plot.}
  \item{scaling}{the scaling factor to apply to the data.}
  \item{continuous}{(used if style="mosaic") a logical specifying
    whether to treat the values in \code{x} as continuous or not; can be
    an integer value giving the number of categories.}
  \item{width}{(used if style = "mosaic") the width of the cells; by
    default, all the available space is used.}
  \item{legend}{(used if style = "mosaic") the place where to draw the
    legend; one of \code{"below"} (the default), \code{"side"}, or
    \code{"none"}, or an unambiguous abbreviation of these.}
  \item{funcol}{(used if style = "mosaic") the function used to generate
    the colours (see details and examples).}
  \item{\dots}{further arguments passed to the graphical functions.}
}
\details{
  The possible values for \code{style} are ``bars'', ``segments'',
  ``image'', ``arrows'', ``boxplot'', ``dotchart'', or ``mosaic'' for
  \code{phydataplot}, and ``ring'', ``segments'', or ``arrows'' for
  \code{ring}.

  \code{style = "image"} works only with square matrices (e.g.,
  similarities). If you want to plot a DNA alignment in the same way
  than \code{\link{image.DNAbin}}, try \code{style = "mosaic"}.

  \code{style = "mosaic"} can plot any kind of matrices, possibly after
  discretizing its values (using \code{continuous}). The default colour
  palette is taken from the function \code{\link[grDevices]{rainbow}}.
  If you want to use specified colours, a function simply returning the
  vector of colours must be used, possibly with names if you want to
  assign a specific colour to each value (see examples).
}
\note{
  For the moment, only rightwards trees are supported (does not apply to
  circular trees).
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{nodelabels}},
  \code{\link{fancyarrows}}
}
\examples{
## demonstrates matching with names:
tr <- rcoal(n <- 10)
x <- 1:n
names(x) <- tr$tip.label
plot(tr, x.lim = 11)
phydataplot(x, tr)
## shuffle x but matching names with tip labels reorders them:
phydataplot(sample(x), tr, "s", lwd = 3, lty = 3)

## adapts to the tree:
plot(tr, "f", x.l = c(-11, 11), y.l = c(-11, 11))
phydataplot(x, tr, "s")

## leave more space with x.lim to show a barplot and a dotchart:
plot(tr, x.lim = 22)
phydataplot(x, tr, col = "yellow")
phydataplot(x, tr, "d", offset = 13)

ts <- rcoal(N <- 100)
X <- rTraitCont(ts) # names are set
dd <- dist(X)
op <- par(mar = rep(0, 4))
plot(ts, x.lim = 10, cex = 0.4, font = 1)
phydataplot(as.matrix(dd), ts, "i", offset = 0.2)

par(xpd = TRUE, mar = op$mar)
co <- c("blue", "red"); l <- c(-2, 2)
X <- X + abs(min(X)) # move scale so X >= 0
plot(ts, "f", show.tip.label = FALSE, x.lim = l, y.lim = l, open.angle = 30)
phydataplot(X, ts, "s", col = co, offset = 0.05)
ring(X, ts, "ring", col = co, offset = max(X) + 0.1) # the same info as a ring

## as many rings as you want...
co <- c("blue", "yellow")
plot(ts, "r", show.tip.label = FALSE, x.l = c(-1, 1), y.l = c(-1, 1))
for (o in seq(0, 0.4, 0.2)) {
    co <- rev(co)
    ring(0.2, ts, "r", col = rep(co, each = 5), offset = o)
}

lim <- c(-5, 5)
co <- rgb(0, 0.4, 1, alpha = 0.1)
y <- seq(0.01, 1, 0.01)
plot(ts, "f", x.lim = lim, y.lim = lim, show.tip.label = FALSE)
ring(y, ts, offset = 0, col = co, lwd = 0.1)
for (i in 1:3) {
    y <- y + 1
    ring(y, ts, offset = 0, col = co, lwd = 0.1)
}

## rings can be in the background
plot(ts, "r", plot = FALSE)
ring(1, ts, "r", col = rainbow(100), offset = -1)
par(new = TRUE)
plot(ts, "r", font = 1, edge.color = "white")

## might be more useful:
co <- c("lightblue", "yellow")
plot(ts, "r", plot = FALSE)
ring(0.1, ts, "r", col = sample(co, size = N, rep = TRUE), offset = -.1)
par(new = TRUE)
plot(ts, "r", font = 1)

## if x is matrix:
tx <- rcoal(m <- 20)
X <- runif(m, 0, 0.5); Y <- runif(m, 0, 0.5)
X <- cbind(X, Y, 1 - X - Y)
rownames(X) <- tx$tip.label
plot(tx, x.lim = 6)
co <- rgb(diag(3))
phydataplot(X, tx, col = co)
## a variation:
plot(tx, show.tip.label = FALSE, x.lim = 5)
phydataplot(X, tx, col = co, offset = 0.05, border = NA)

plot(tx, "f", show.tip.label = FALSE, open.angle = 180)
ring(X, tx, col = co, offset = 0.05)

Z <- matrix(rnorm(m * 5), m)
rownames(Z) <- rownames(X)
plot(tx, x.lim = 5)
phydataplot(Z, tx, "bo", scaling = .5, offset = 0.5,
            boxfill = c("gold", "skyblue"))

## plot an alignment with a NJ tree:
data(woodmouse)
trw <- nj(dist.dna(woodmouse))
plot(trw, x.lim = 0.1, align.tip = TRUE, font = 1)
phydataplot(woodmouse[, 1:50], trw, "m", 0.02, border = NA)

## use type = "mosaic" on a 30x5 matrix:
tr <- rtree(n <- 30)
p <- 5
x <- matrix(sample(3, size = n*p, replace = TRUE), n, p)
dimnames(x) <- list(paste0("t", 1:n), LETTERS[1:p])
plot(tr, x.lim = 35, align.tip = TRUE, adj = 1)
phydataplot(x, tr, "m", 2)
## change the aspect:
plot(tr, x.lim = 35, align.tip = TRUE, adj = 1)
phydataplot(x, tr, "m", 2, width = 2, border = "white", lwd = 3, legend = "side")
## user-defined colour:
f <- function(n) c("yellow", "blue", "red")
phydataplot(x, tr, "m", 18, width = 2, border = "white", lwd = 3,
            legend = "side", funcol = f)

## alternative colour function...:
## fb <- function(n) c("3" = "red", "2" = "blue", "1" = "yellow")
## ... but since the values are sorted alphabetically,
## both f and fb will produce the same plot.

## use continuous = TRUE with two different scales:
x[] <- 1:(n*p)
plot(tr, x.lim = 35, align.tip = TRUE, adj = 1)
phydataplot(x, tr, "m", 2, width = 1.5, continuous = TRUE, legend = "side",
            funcol = colorRampPalette(c("white", "darkgreen")))
phydataplot(x, tr, "m", 18, width = 1.5, continuous = 5, legend = "side",
            funcol = topo.colors)
}
\keyword{aplot}
