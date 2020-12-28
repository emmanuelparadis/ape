## ladderize.R (2017-04-25)

##   Ladderize a Tree

## Copyright 2007-2017 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

ladderize <- function(phy, right = TRUE)
{
    foo <- function(node, END, where) {
        start <- which(phy$edge[, 1] == node)
        end <- c(start[-1] - 1, END)
        size <- end - start + 1
        desc <- phy$edge[start, 2]
        Nclade <- length(desc)
        n <- N[desc]
        o <- order(n, decreasing = right)
        newpos <- c(0, cumsum(size[o][-Nclade])) + where
        desc <- desc[o]
        end <- end[o]
        start <- start[o]
        neworder[newpos] <<- start
        for (i in 1:Nclade)
            if (desc[i] > nb.tip) foo(desc[i], end[i], newpos[i] + 1)
    }
    phy <- reorder(phy) # fix by Klaus (2015-10-04)
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    nb.edge <- dim(phy$edge)[1]
    tmp <- reorder(phy, "postorder")
    N <- .C(node_depth, as.integer(nb.tip),
            as.integer(tmp$edge[, 1]), as.integer(tmp$edge[, 2]),
            as.integer(nb.edge), double(nb.tip + nb.node), 1L)[[5]]
    neworder <- integer(nb.edge)
    foo(nb.tip + 1, nb.edge, 1)
    phy$edge <- phy$edge[neworder, ]
    if (!is.null(phy$edge.length))
        phy$edge.length <- phy$edge.length[neworder]
    phy
}
