## collapsed.intervals.R (2002-09-12)

##   Collapsed coalescent intervals (e.g. for the skyline plot)

## Copyright 2002 Korbinian Strimmer

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

# construct collapsed intervals from coalescent intervals
collapsed.intervals <- function(ci, epsilon=0.0)
{
  if (class(ci) != "coalescentIntervals")
    stop("object \"ci\" is not of class \"coalescentIntervals\"")

  sz <- ci$interval.length
  lsz <- length(sz)
  idx <- c <- 1:lsz

  p <- 1
  w <- 0

  # starting from tips collapes intervals
  # until total size is >= epsilon
  for (i in 1:lsz)
  {
    idx[[i]] <- p
    w <- w + sz[[i]]
    if (w >= epsilon)
    {
      p <- p+1
      w <- 0
    }
  }

  # if last interval is smaller than epsilon merge
  # with second last interval
  lastInterval <- idx==p
  if ( sum(sz[lastInterval]) < epsilon )
  {
    p <- p-1
    idx[lastInterval] <- p
  }

  obj <- list(
     lineages=ci$lineages,
     interval.length=ci$interval.length,
     collapsed.interval=idx, # collapsed intervals (via reference)
     interval.count=ci$interval.count,
     collapsed.interval.count = idx[[ci$interval.count]],
     total.depth =ci$total.depth,
     epsilon = epsilon
    )
  class(obj) <- "collapsedIntervals"

  return(obj)
}
