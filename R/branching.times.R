## branching.times.R (2018-01-16)

##    Branching Times of a Phylogenetic Tree

## Copyright 2002-2018 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

branching.times <- function(phy)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    phy <- reorder(phy)
    n <- length(phy$tip.label)
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    EL <- phy$edge.length
    if (is.null(EL)) {
        warning("no branch length in tree")
        return(numeric())
    }
    N <- length(e1)
    xx <- numeric(phy$Nnode)
    interns <- which(e2 > n)
    ## we loop only on the internal edges, this assumes
    ## that `xx' is already set with 0
    for (i in interns) xx[e2[i] - n] <- xx[e1[i] - n] + EL[i]
    depth <- xx[e1[N] - n] + EL[N]
    xx <- depth - xx
    names(xx) <-
        if (is.null(phy$node.label)) (n + 1):(n + phy$Nnode)
        else phy$node.label
    xx
}
