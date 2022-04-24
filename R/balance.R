## balance.R (2015-07-08)

##   Balance of a Dichotomous Phylogenetic Tree

## Copyright 2002-2015 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

balance <- function(phy)
{
    if (!inherits(phy, "phylo"))
      stop('object "phy" is not of class "phylo"')
    phy <- reorder(phy, "postorder")
    N <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != N - 1)
      stop('"phy" is not rooted and fully dichotomous')
    ans <- matrix(NA, nb.node, 2)
    nd <- node.depth(phy)
    i <- 1
    while(i < nrow(phy$edge)){
        node <- phy$edge[i, 1] - N
        ans[node, 1] <- nd[phy$edge[i,2]]
        ans[node, 2] <- nd[phy$edge[i+1,2]]
        i <- i+2
    }
    rownames(ans) <-
      if (is.null(phy$node.label)) N + 1:nb.node else phy$node.label
    ans
}
