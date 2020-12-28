## coalescent.intervals.R (2002-09-12)

##   Constructs objects with information on coalescent intervals

## Copyright 2002 Korbinian Strimmer

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

coalescent.intervals <- function(x) UseMethod("coalescent.intervals")

# set up coalescent interval object (from NH tree)
coalescent.intervals.phylo <- function(x)
{
    if (class(x) != "phylo") stop("object \"x\" is not of class \"phylo\"")

    # ensure we have a BINARY tree
    if (!is.binary.phylo(x)) stop("object \"x\" is not a binary tree")
    # ordered branching times
    t <- sort(branching.times(x))
    lt <- length(t)

    # interval widths
    w <- numeric(lt)
    w[1] <- t[1]
    for (i in 2:lt) w[i] <- t[i] - t[i - 1]

    l <- (lt+1):2       # number of lineages

    obj <- list(
     lineages=l,
     interval.length=w,
     interval.count=lt,
     total.depth =sum(w))
    class(obj) <- "coalescentIntervals"
    return(obj)
}


# set up coalescent interval object from vector of interval length
coalescent.intervals.default <- function(x)
{
  if (!is.vector(x)) stop("argument \"x\" is not a vector of interval lengths")

  # x = list of the widths of each interval
  lt <- length(x)
  l <- (lt+1):2           # number of lineages at the beginning of each interval

  obj <- list(
     lineages=l,
     interval.length=x,
     interval.count=lt,
     total.depth =sum(x))
    class(obj) <- "coalescentIntervals"
    return(obj)
}
