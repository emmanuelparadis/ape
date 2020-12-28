## is.binary.tree.R (2016-11-03)

##    Test for Binary Tree

## Copyright 2016 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.binary <- function(phy) UseMethod("is.binary")

is.binary.phylo <- function(phy)
    length(phy$tip.label) - phy$Nnode + is.rooted.phylo(phy) == 2

is.binary.tree <- function(phy)
{
    ##warning("is.binary.tree() is deprecated; using is.binary() instead.\n\nis.binary.tree() will be removed soon: see ?is.binary and update your code.")
    is.binary(phy)
}

is.binary.multiPhylo <- function(phy)
{
    phy <- unclass(phy)
    n <- length(attr(phy, "TipLabel"))
    if (n)
        n - sapply(phy, "[[", "Nnode") + is.rooted.multiPhylo(phy) == 2
    else
        sapply(phy, is.binary.phylo)
}
