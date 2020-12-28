## balance.R (2015-07-08)

##   Balance of a Dichotomous Phylogenetic Tree

## Copyright 2002-2015 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

balance <- function(phy)
{
    if (!inherits(phy, "phylo"))
      stop('object "phy" is not of class "phylo"')
    phy <- reorder(phy) # fix a bug reported by G. Valiente in January 2009 (2015-07-08)
    N <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != N - 1)
      stop('"phy" is not rooted and fully dichotomous')
    ans <- matrix(NA, nb.node, 2)
    foo <- function(node, n) {
        s <- which(phy$edge[, 1] == node)
        desc <- phy$edge[s, 2]
        ans[node - N, 1] <<- n1 <- (s[2] - s[1] + 1)/2
        ans[node - N, 2] <<- n2 <- n - n1
        if (desc[1] > N) foo(desc[1], n1)
        if (desc[2] > N) foo(desc[2], n2)
    }
    foo(N + 1, N)
    rownames(ans) <-
      if (is.null(phy$node.label)) N + 1:nb.node else phy$node.label
    ans
}
