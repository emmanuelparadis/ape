## is.binary.tree.R (2023-02-07)

##    Test for Binary Tree

## Copyright 2016-2023 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.binary <- function(phy) UseMethod("is.binary")

is.binary.phylo <- function(phy)
{
    n <- length(phy$tip.label)
    m <- phy$Nnode
    dgr <- tabulate(phy$edge, n + m)
    ref <- c(rep.int(1L, n), rep.int(3L, m))
    ## the root is assumed to be numbered n+1
    if (.is.rooted_ape(phy, n)) ref[n + 1L] <- 2L
    ## can use identical() as long as tabulate() returns integers
    identical(dgr, ref)
}

is.binary.tree <- function(phy)
{
    message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nis.binary.tree() is deprecated; using is.binary() instead.\n\nis.binary.tree() will be removed soon: see ?is.binary and update your code.\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
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
