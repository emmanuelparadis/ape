\name{skylineplot}
\alias{skylineplot}
\alias{plot.skyline}
\alias{lines.skyline}
\alias{skylineplot.deluxe}
\title{Drawing Skyline Plot Graphs}
\usage{
\method{plot}{skyline}(x, show.years=FALSE, subst.rate, present.year, \dots)
\method{lines}{skyline}(x, show.years=FALSE, subst.rate, present.year, \dots)
skylineplot(z, \dots)
skylineplot.deluxe(tree, \dots)
}
\arguments{
  \item{x}{skyline plot data (i.e. an object of class \code{"skyline"}).}
  \item{z}{Either an ultrametric tree (i.e. an object of class \code{"phylo"}),
           or coalescent intervals (i.e. an object of class \code{"coalescentIntervals"}), or
	   collapsed coalescent intervals (i.e. an object of class \code{"collapsedIntervals"}).}
  \item{tree}{ultrametric tree (i.e. an object of class \code{"phylo"}).}
  \item{show.years}{option that determines whether the time is plotted in units of
        of substitutions (default) or in years (requires specification of substution rate
	and year of present).}
 \item{subst.rate}{substitution rate (see option show.years).}
 \item{present.year}{present year (see option show.years).}

  \item{\dots}{further arguments to be passed on to \code{skyline()} and \code{plot()}.}
}
\description{
 These functions provide various ways to draw \emph{skyline plot} graphs
 on the current graphical device. Note that \code{skylineplot(z, \dots)} is simply
 a shortcut for \code{plot(skyline(z, \dots))}.
 The skyline plot itself is an estimate of effective population size through time,
 and is computed using the function \code{\link{skyline}}.
}

\details{
 See \code{\link{skyline}} for more details (incl. references) about the skyline plot method.
}
\author{Korbinian Strimmer}
\seealso{
\code{\link[graphics]{plot}} and \code{\link[graphics]{lines}} for the basic plotting
function in R, \code{\link{coalescent.intervals}}, \code{\link{skyline}}
}
\examples{
# get tree
data("hivtree.newick") # example tree in NH format
tree.hiv <- read.tree(text = hivtree.newick) # load tree

#### classic skyline plot
skylineplot(tree.hiv) # shortcut

#### plot classic and generalized skyline plots and estimate epsilon
sk.opt <- skylineplot.deluxe(tree.hiv)
sk.opt$epsilon

#### classic and generalized skyline plot ####
sk1 <- skyline(tree.hiv)
sk2 <- skyline(tree.hiv, 0.0119)

# use years rather than substitutions as unit for the time axis
plot(sk1, show.years=TRUE, subst.rate=0.0023, present.year = 1997, col=c(grey(.8),1))
lines(sk2,  show.years=TRUE, subst.rate=0.0023, present.year = 1997)
legend(.15,500, c("classic", "generalized"), col=c(grey(.8),1),lty=1)

#### various skyline plots for different epsilons
layout(mat= matrix(1:6,2,3,byrow=TRUE))
ci <- coalescent.intervals(tree.hiv)
plot(skyline(ci, 0.0));title(main="0.0")
plot(skyline(ci, 0.007));title(main="0.007")
plot(skyline(ci, 0.0119),col=4);title(main="0.0119")
plot(skyline(ci, 0.02));title(main="0.02")
plot(skyline(ci, 0.05));title(main="0.05")
plot(skyline(ci, 0.1));title(main="0.1")
layout(mat= matrix(1:1,1,1,byrow=TRUE))
}
\keyword{hplot}
