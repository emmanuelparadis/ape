\name{collapsed.intervals}
\alias{collapsed.intervals}
\title{Collapsed Coalescent Intervals}
\usage{
collapsed.intervals(ci, epsilon=0)
}
\arguments{
  \item{ci}{coalescent intervals (i.e. an object of class \code{"coalescentIntervals"}).}
  \item{epsilon}{collapsing parameter that controls the amount of smoothing
  (allowed range: from \code{0} to \code{ci$total.depth})}
}
\description{
 This function takes a \code{"coalescentIntervals"} objects and collapses neighbouring
 coalescent intervals into a single combined interval so that every collapsed interval is
 larger than \code{epsilon}. Collapsed coalescent intervals are used, e.g., to obtain the
 generalized skyline plot (\code{\link{skyline}}). For \code{epsilon = 0} no interval
 is collapsed.
}
\details{
Proceeding from the tips to the root of the tree each small
interval is pooled with the neighboring interval closer to the root. If the
neighboring interval is also small, then pooling continues until the composite
interval is larger than \code{epsilon}. Note that this approach prevents the
occurrence of zero-length intervals at the present.
For more details see Strimmer and Pybus (2001).
}
\value{
An object of class \code{"collapsedIntervals"} with the following entries:

  \item{lineages}{ A vector with the number of lineages at the start of each coalescent
    interval.}
  \item{interval.length}{ A vector with the length of each coalescent
    interval.}
   \item{collapsed.interval}{A vector indicating for each coalescent interval to which
     collapsed interval it belongs.}
  \item{interval.count}{ The total number of coalescent
    intervals.}
   \item{collapsed.interval.count}{The number of collapsed intervals.}
  \item{total.depth}{ The sum of the lengths of all coalescent
    intervals.}
  \item{epsilon}{The value of the underlying smoothing parameter.}
}
\author{Korbinian Strimmer}
\seealso{
\code{\link{coalescent.intervals}},\code{\link{skyline}}.
}
\references{
  Strimmer, K. and Pybus, O. G. (2001) Exploring the demographic history
  of DNA sequences using the generalized skyline plot. \emph{Molecular
    Biology and Evolution}, \bold{18}, 2298--2305.
}
\examples{
data("hivtree.table") # example tree
# colescent intervals from vector of interval lengths
ci <- coalescent.intervals(hivtree.table$size)
ci
# collapsed intervals
cl1 <- collapsed.intervals(ci,0)
cl2 <- collapsed.intervals(ci,0.0119)
cl1
cl2
}
\keyword{manip}
