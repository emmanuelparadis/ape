## skyline.R (2002-09-12)

##   Methods to construct skyline objects (data underlying skyline plot)

## Copyright 2002 Korbinian Strimmer

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

skyline <- function(x, ...) UseMethod("skyline")

# input: phylogenetic tree
skyline.phylo <- function(x, ...)
{
  if (class(x) != "phylo")
    stop("object \"x\" is not of class \"phylo\"")

  skyline(coalescent.intervals(x), ...)
}

# input: coalescent intervals and epsilon
skyline.coalescentIntervals <- function(x, epsilon=0, ...)
{
  if (class(x) != "coalescentIntervals")
    stop("object \"x\" is not of class \"coalescentIntervals\"")

  if (epsilon < 0)
  {
    eps <- find.skyline.epsilon(x, ...)
  }
  else
    eps <- epsilon

  skyline(collapsed.intervals(x, epsilon=eps), ...)
}


# input: collapsed intervals
skyline.collapsedIntervals <- function(x, old.style=FALSE, ...)
{
  if (class(x) != "collapsedIntervals")
    stop("object \"x\" is not of class \"collapsedIntervals\"")

  link <- x$collapsed.interval
  params <- x$collapsed.interval.count
  l <- x$lineages
  w <- x$interval.length

  b <- choose(l,2) # binomial coefficients

  sg <- rep(0,params)   # sizes of collapsed intervals
  cg <- rep(0,params)   # coalescent events in interval

  if(old.style)
    ng <- rep(0,params) # lineages at beginning of an in interval
  else
  {
    ng <- rep(0,params) # sum of classic skp estimates in an interval
    m.classic <- w*b
  }

  for (i in 1:params)
  {
    group <- link==i
    sgr <- w[group]
    sg[[i]] <- sum(sgr)
    cg[[i]] <- length(sgr)

    if(old.style)
      ng[[i]] <- l[group][[1]]
    else
      ng[[i]] <- sum(m.classic[group])
  }

  # generalized skp estimate
  t <- cumsum(sg)
  if (old.style)
    m <- sg*(ng*(ng-cg)/(2.0*cg) )
  else
    m <- ng/cg

  # log-likelihood
  logL <- sum(log(b/m[link]) - b/m[link]*w)

  # AICc corrected log-likelihood
  K <- x$collapsed.interval.count
  S <- x$interval.count
  if (S-K > 1)
    logL.AICc <- logL - K- K*(K+1)/(S-K-1)
  else
    logL.AICc <- NA

  obj <- list(
    time=t,
    interval.length=sg,
    population.size=m,
    parameter.count=length(t),
    epsilon = x$epsilon,
    logL = logL,
    logL.AICc = logL.AICc
  )
  class(obj) <- "skyline"
  return(obj)
}

# grid search for finding optimal epsilon parameter
find.skyline.epsilon <- function(ci, GRID=1000, MINEPS=1e-6, ...)
{
  # Why MINEPS?
  # Because most "clock-like" trees are not properly
  # clock-like for a variety of reasons, i.e. the heights
  # of the tips are not exactly zero.

  cat("Searching for the optimal epsilon... ")

  # a grid search is a naive way but still effective of doing this ...

  size <- ci$interval.count
  besteps <- ci$total.depth
  eps <- besteps

  cli <- collapsed.intervals(ci,eps)
  skpk <- skyline(cli, ...)
  bestaicc <- skpk$ logL.AICc
  params <- skpk$parameter.count

  delta <- besteps/GRID

  eps <- eps-delta
  while(eps > MINEPS)
  {
    cli <- collapsed.intervals(ci,eps)
    skpk <- skyline(cli, ...)
    aicc <- skpk$ logL.AICc
    params <- skpk$parameter.count

    if (aicc > bestaicc && params < size-1)
    {
      besteps <- eps
      bestaicc <- aicc
    }
    eps <- eps-delta
  }

   cat("epsilon =", besteps, "\n")

  besteps
}

