\name{axisPhylo}
\alias{axisPhylo}
\title{Axis on Side of Phylogeny}
\usage{
axisPhylo(side = 1, root.time = NULL, backward = TRUE, ...)
}
\arguments{
  \item{side}{a numeric value specifying the side where the axis is
    plotted: 1: below, 2: left, 3: above, 4: right.}
  \item{root.time}{the time assigned to the root node of the tree. By
    default, this is taken from the \code{root.time} element of the
    tree. If it is absent, this is determined from the next option.}
  \item{backward}{a logical value; if TRUE, the most distant tip from
    the root is considered as the origin of the time scale; if FALSE,
    this is the root node.}
  \item{\dots}{further arguments to be passed to \code{axis}.}
}
\description{
  This function adds a scaled axis on the side of a phylogeny plot.
}
\details{
  The further arguments (\code{...}) are used to format the axis. They
  may be \code{font}, \code{cex}, \code{col}, \code{las}, and so on (see
  the help pages on \code{\link[graphics]{axis}} and
  \code{\link[graphics]{par}}).
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{add.scale.bar}},
  \code{\link[graphics]{axis}}, \code{\link[graphics]{par}}
}
\examples{
tr <- rtree(30)
ch <- rcoal(30)
plot(ch)
axisPhylo()
plot(tr, "c", FALSE, direction = "u")
axisPhylo(2, las = 1)
}
\keyword{aplot}
