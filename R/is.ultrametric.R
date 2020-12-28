## is.ultrametric.R (2016-10-04)

##   Test if a Tree is Ultrametric

## Copyright 2003-2016 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.ultrametric <- function(phy, ...) UseMethod("is.ultrametric")

## the main driver code (n = number of tips):
.is.ultrametric_ape <- function(phy, tol, option, n)
{
    if (is.null(phy$edge.length))
        stop("the tree has no branch lengths")
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    EL <- phy$edge.length

    ## xx: distance from a node or a tip to the root
    xx <- numeric(n + phy$Nnode)

    ## the following must start at the root and follow the
    ## edges contiguously; so the tree must be either in cladewise
    ## order (or in pruningwise but the for loop must start from
    ## the bottom of the edge matrix)

    for (i in seq_len(length(e1)))
        xx[e2[i]] <- xx[e1[i]] + EL[i]

    xx.tip <- xx[1:n]

    crit <- switch(option, {
        mn <- min(xx.tip)
        mx <- max(xx.tip)
        (mx - mn)/mx
    }, var(xx.tip))

    isTRUE(all.equal.numeric(crit, 0, tolerance = tol))
}

is.ultrametric.phylo <- function(phy, tol = .Machine$double.eps^0.5,
                                 option = 1, ...)
{
    phy <- reorder.phylo(phy)
    .is.ultrametric_ape(phy, tol, option, length(phy$tip.label))
}

is.ultrametric.multiPhylo <- function(phy, tol = .Machine$double.eps^0.5,
                                      option = 1, ...)
{
    phy <- reorder.multiPhylo(phy)
    labs <- attr(phy, "TipLabel")
    if (is.null(labs))
        sapply(phy, is.ultrametric.phylo, tol = tol, option = option)
    else
        sapply(phy, .is.ultrametric_ape, tol = tol, option = option, n = length(labs))
}
